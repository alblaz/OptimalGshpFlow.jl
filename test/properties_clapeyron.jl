using ForwardDiff

# Clapeyron-based refrigerant properties checked against CoolProp (same reference EoS).
@testset "ClapeyronRefrigerant thermodynamics" begin
    for f in ("R290", "R134a", "R32")
        r  = ClapeyronRefrigerant(f)
        pc = PropsSI("pcrit", f)
        @test r.pcrit ≈ pc rtol = 1e-10
        @test r.Tcrit ≈ PropsSI("Tcrit", f) rtol = 1e-10
        for pr in (0.02, 0.05, 0.1, 0.2, 0.4, 0.7)
            p = pr * pc
            @test T_dew(r, p)    ≈ PropsSI("T", "P", p, "Q", 1, f) rtol = 1e-8
            @test T_bubble(r, p) ≈ PropsSI("T", "P", p, "Q", 0, f) rtol = 1e-8
            @test h_dew(r, p)    ≈ PropsSI("H", "P", p, "Q", 1, f) rtol = 1e-8
            @test h_bubble(r, p) ≈ PropsSI("H", "P", p, "Q", 0, f) rtol = 1e-8
            @test ρ_dew(r, p)    ≈ PropsSI("D", "P", p, "Q", 1, f) rtol = 1e-8
            @test ρ_bubble(r, p) ≈ PropsSI("D", "P", p, "Q", 0, f) rtol = 1e-8
            @test s_dew(r, p)    ≈ PropsSI("S", "P", p, "Q", 1, f) rtol = 1e-8
            @test cp_dew(r, p)   ≈ PropsSI("C", "P", p, "Q", 1, f) rtol = 1e-6
            Ts = T_dew(r, p)
            for ΔT in (2., 5., 20.)
                T = Ts + ΔT
                @test h_pT(r, p, T)  ≈ PropsSI("H", "P", p, "T", T, f) rtol = 1e-8
                @test s_pT(r, p, T)  ≈ PropsSI("S", "P", p, "T", T, f) rtol = 1e-8
                @test ρ_pT(r, p, T)  ≈ PropsSI("D", "P", p, "T", T, f) rtol = 1e-8
                @test cp_pT(r, p, T) ≈ PropsSI("C", "P", p, "T", T, f) rtol = 1e-6
            end
            # isentropic compression from 5 K superheat to 3 × p
            s = s_pT(r, p, Ts + 5)
            p2 = min(3p, 0.9pc)
            @test h_ps(r, p2, s) ≈ PropsSI("H", "P", p2, "S", s, f) rtol = 1e-8
        end
    end
end

@testset "ClapeyronRefrigerant transport" begin
    # defaults: CoolProp for pure fluids, EntropyScaling for mixtures
    @test ClapeyronRefrigerant("R290").transport isa CoolPropTransport
    @test ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.5, 0.5]).transport isa EntropyScalingTransport
    # CoolProp backend: identical to CoolProp
    for f in ("R290", "R134a", "R32")
        r = ClapeyronRefrigerant(f)
        p = 0.1 * r.pcrit; T = T_dew(r, p) + 5
        @test μ_pT(r, p, T) ≈ PropsSI("V", "P", p, "T", T, f) rtol = 1e-12
        @test λ_pT(r, p, T) ≈ PropsSI("CONDUCTIVITY", "P", p, "T", T, f) rtol = 1e-12
        @test μ_dew(r, p)   ≈ PropsSI("V", "P", p, "Q", 1, f) rtol = 1e-12
        @test λ_dew(r, p)   ≈ PropsSI("CONDUCTIVITY", "P", p, "Q", 1, f) rtol = 1e-12
    end
    # EntropyScaling backend (RefpropRES): sanity bounds only. Measured deviations vs CoolProp in
    # the vapour: μ within ≈ 9 %, λ up to ≈ 16 % (R290); see examples/properties_validation.jl
    for f in ("R290", "R134a", "R32")
        r = ClapeyronRefrigerant(f; transport = EntropyScalingTransport())
        for pr in (0.05, 0.15, 0.4)
            p = pr * r.pcrit; T = T_dew(r, p) + 5
            @test μ_pT(r, p, T; phase = :vapor) ≈ PropsSI("V", "P", p, "T", T, f) rtol = 0.10
            @test λ_pT(r, p, T; phase = :vapor) ≈ PropsSI("CONDUCTIVITY", "P", p, "T", T, f) rtol = 0.20
        end
    end
    @test_throws ArgumentError ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.5, 0.5], transport = CoolPropTransport())
end

@testset "ClapeyronRefrigerant automatic differentiation" begin
    r  = ClapeyronRefrigerant("R290"; transport = EntropyScalingTransport())
    rc = ClapeyronRefrigerant("R290")       # CoolProp transport
    p  = 4.7e5
    T  = T_dew(r, p) + 5
    s  = s_pT(r, p, T)
    fd(f, x, h) = (f(x + h) - f(x - h)) / 2h
    for (f, x, h) in ((q -> T_dew(r, q), p, 1.),
                      (q -> h_dew(r, q), p, 1.),
                      (q -> ρ_dew(r, q), p, 1.),
                      (t -> h_pT(r, p, t), T, 1e-3),
                      (t -> cp_pT(r, p, t), T, 1e-3),
                      (x -> h_ps(r, 15e5, x), s, 1e-3),
                      (t -> μ_pT(r, p, t), T, 1e-3),
                      (t -> λ_pT(r, p, t), T, 1e-3),
                      (t -> λ_pT(rc, p, t), T, 1e-3),
                      (q -> λ_dew(rc, q), p, 1.))
        @test ForwardDiff.derivative(f, x) ≈ fd(f, x, h) rtol = 1e-4
    end
    # dh/ds at constant p is the temperature
    @test ForwardDiff.derivative(x -> h_ps(r, 15e5, x), s) ≈ T_ps(r, 15e5, s) rtol = 1e-8
    # phase hint and initial guess give the same result
    @test h_ps(r, 15e5, s; phase = :vapor, T0 = 330.) ≈ h_ps(r, 15e5, s) rtol = 1e-10
    @test h_pT(r, p, T; phase = :vapor) ≈ h_pT(r, p, T) rtol = 1e-12
end

@testset "ClapeyronRefrigerant mixture (R454C)" begin
    w = [0.215, 0.785]
    mx = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = w, name = "R454C")
    # mole fractions reproduce the mass fractions
    Mi = [PropsSI("M", "R32"), PropsSI("M", "R1234yf")]
    @test (mx.z .* Mi) ./ sum(mx.z .* Mi) ≈ w rtol = 1e-4
    @test mx.M ≈ sum(mx.z .* Mi) rtol = 1e-4
    p = 3e5
    @test T_bubble(mx, p) < T_dew(mx, p)                  # zeotropic: positive glide
    @test 0 < cp_twophase(mx, p) < Inf
    @test h_bubble(mx, p) < h_dew(mx, p)
    # dew and bubble states are consistent with single-phase evaluations just outside the dome
    @test h_pT(mx, p, T_dew(mx, p) + 1e-3; phase = :vapor) ≈ h_dew(mx, p) rtol = 1e-5
    @test h_pT(mx, p, T_bubble(mx, p) - 1e-3; phase = :liquid) ≈ h_bubble(mx, p) rtol = 1e-5
    @test ForwardDiff.derivative(q -> T_dew(mx, q), p) ≈ (T_dew(mx, p + 1) - T_dew(mx, p - 1)) / 2 rtol = 1e-4
    @test cp_twophase(ClapeyronRefrigerant("R290"), p) == Inf
    @test_throws ArgumentError ClapeyronRefrigerant(["R32", "R1234yf"])
end
