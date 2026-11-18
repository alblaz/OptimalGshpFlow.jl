# Infinite line source and temporal superposition (ROADMAP step 8)
@testset "infinite line source" begin
    # large Fo: ½ E₁(x) → ½(−γ − ln x), x = 1/(4Fo)
    for Fo in (1e3, 1e5)
        @test ILS(Fo) ≈ 0.5 * (log(4Fo) - Base.MathConstants.eulergamma) rtol = 1e-3
    end
    @test ILS(3600., 0.06, 1e-6) == ILS(1e-6 * 3600. / 0.06^2)
    @test ILS(3600, 0.06, 1e-6) == ILS(3600., 0.06, 1e-6)          # mixed number types

    g = InfiniteLineSource(λ = 3, α = 1e-6, rb = 0.06, Tg = 281.15)
    day = 86400.
    @test borehole_wall_temperature(g, Float64[], Float64[], 10day) == g.Tg     # no load
    # one step of 30 W/m from t = 0
    Tb1 = borehole_wall_temperature(g, [0.], [30.], 10day)
    @test Tb1 ≈ g.Tg - 30/(2π*3) * ILS(10day, 0.06, 1e-6)
    # loads not yet started have no effect
    @test borehole_wall_temperature(g, [0., 20day], [30., 0.], 10day) == Tb1
    # on for 10 days then off = difference of two responses (recovery)
    Tb2 = borehole_wall_temperature(g, [0., 10day], [30., 0.], 15day)
    @test Tb2 ≈ g.Tg - 30/(2π*3) * (ILS(15day, 0.06, 1e-6) - ILS(5day, 0.06, 1e-6))
    @test Tb1 < Tb2 < g.Tg
    @test (@inferred borehole_wall_temperature(g, [0., 10day], [30., 20.], 15day)) isa Float64

    # tabulated response for a uniform step = general superposition
    H = UniformLoadHistory(g, 3600.)
    @test borehole_wall_temperature(H) == g.Tg
    qs = [30., 28., 0., 25., 25., 31.]
    for (n, qn) in enumerate(qs)
        push_load!(H, qn)
        @test borehole_wall_temperature(H) ≈ borehole_wall_temperature(g, 3600. .* (0:n-1), qs[1:n], n * 3600.) rtol = 1e-13
    end
    @test (@inferred borehole_wall_temperature(H)) isa Float64
end
