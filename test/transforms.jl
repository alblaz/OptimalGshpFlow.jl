using ForwardDiff

# analytic dx/du
σ(u) = 1 / (1 + exp(-u))
dfrom_R(::Free, u) = 1.0
dfrom_R(t::Linear, u) = t.scale
dfrom_R(::Unit, u) = σ(u) * σ(-u)
dfrom_R(t::Box, u) = (t.hi - t.lo) * σ(u) * σ(-u)
dfrom_R(t::LogBox, u) = from_R(t, u) * (log(t.hi) - log(t.lo)) * σ(u) * σ(-u)
dfrom_R(t::Above, u) = t.scale * exp(u)

@testset "variable transforms" begin
    us = range(-30., 30., 61)
    ts = (Free(), Linear(273.15, 5.), Unit(), Box(250., 400.), LogBox(1e3, 4.25e6), LogBox(0.000171, 4.25e6), Above(6e-5), Above(-1., 2.))
    for t in ts
        lo, hi = bounds(t)
        for u in us
            x = from_R(t, u)
            @test lo ≤ x ≤ hi                       # closed: the sigmoid saturates in Float64
            # Duals pass through; derivative against the analytic one
            @test ForwardDiff.derivative(v -> from_R(t, v), u) ≈ dfrom_R(t, u) rtol = 1e-12 atol = 1e-300
        end
        # round trip x -> u -> x on interior points
        for u in range(-10., 10., 21)
            x = from_R(t, u)
            @test from_R(t, to_R(t, x)) ≈ x rtol = 1e-12
        end
    end

    # values
    @test from_R(Unit(), 0.) == 0.5
    @test from_R(Box(250., 400.), 0.) == 325.
    @test from_R(LogBox(1e3, 1e5), 0.) ≈ 1e4 rtol = 1e-14        # geometric mean
    @test from_R(Above(2., 3.), 0.) == 5.
    @test from_R(Linear(1e3, 2e3), 1.5) == 4e3 && to_R(Linear(1e3, 2e3), 4e3) == 1.5
    @test to_R(Free(), 3.) == 3.
    # LogBox: equal steps in u are equal ratios of x around the midpoint
    t = LogBox(1e3, 1e5)
    @test from_R(t, 0.1) / from_R(t, 0.) ≈ from_R(t, 0.) / from_R(t, -0.1) rtol = 1e-12
    # same mapping as the sigmoid helpers
    @test from_R(Box(1., 3.), 0.7) == sigmoid_R_to_ab(0.7, 1., 3.)
    @test from_R(Unit(), 0.7) == sigmoid_R_to_01(0.7)

    # outside the bounds there is no real unknown
    @test_throws DomainError to_R(Box(0., 1.), 2.)
    @test_throws DomainError to_R(LogBox(1., 10.), -1.)

    # invalid intervals
    @test_throws ArgumentError Box(2., 1.)
    @test_throws ArgumentError Box(0., Inf)
    @test_throws ArgumentError LogBox(0., 1.)
    @test_throws ArgumentError Above(1., 0.)
    @test_throws ArgumentError Linear(0., -1.)
    @test Box(1, 2) isa Box{Float64}
end
