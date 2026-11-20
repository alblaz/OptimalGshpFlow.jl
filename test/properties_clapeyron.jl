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
            @test cp_bubble(r, p) ≈ PropsSI("C", "P", p, "Q", 0, f) rtol = 1e-6
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
        @test μ_bubble(r, p) ≈ PropsSI("V", "P", p, "Q", 0, f) rtol = 1e-12
        @test λ_bubble(r, p) ≈ PropsSI("CONDUCTIVITY", "P", p, "Q", 0, f) rtol = 1e-12
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
        @test dderivative(f, x) ≈ fd(f, x, h) rtol = 1e-4
    end
    # dh/ds at constant p is the temperature
    @test dderivative(x -> h_ps(r, 15e5, x), s) ≈ T_ps(r, 15e5, s) rtol = 1e-8
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
    # saturated liquid of the blend (EntropyScaling): between the pure components at the same T,
    # within the EntropyScaling accuracy (liquid conductivity of a blend is below the linear mix)
    Tb = T_bubble(mx, p)
    for (f, key) in ((μ_bubble, "V"), (λ_bubble, "CONDUCTIVITY"))
        lo, hi = extrema(PropsSI(key, "T", Tb, "Q", 0, c) for c in ("R32", "R1234yf"))
        @test 0.8lo < f(mx, p) < 1.2hi
    end
    @test cp_bubble(mx, p) ≈ cp_pT(mx, p, Tb - 1e-3; phase = :liquid) rtol = 1e-3
    @test dderivative(q -> T_dew(mx, q), p) ≈ (T_dew(mx, p + 1) - T_dew(mx, p - 1)) / 2 rtol = 1e-4
    @test cp_twophase(ClapeyronRefrigerant("R290"), p) == Inf
    @test_throws ArgumentError ClapeyronRefrigerant(["R32", "R1234yf"])
end

# T_ps with phase = :vapor: Newton in T at constant p (no Clapeyron flash: ≈ 16 s less to compile),
# same result as CoolProp and as the flash, exact AD, concrete types; other states fall back.
@testset "T_ps: vapour Newton" begin
    for f in ("R290", "R134a", "R32", "R744")
        r = ClapeyronRefrigerant(f)
        for pr in (0.05, 0.2, 0.5, 0.8), ΔT in (0.5, 10., 40.)
            p = max(pr * r.pcrit, 3r.ptriple)          # R744: triple point at 5.2 bar
            # isentropic compression from p/3 (dew + ΔT) to p: superheated outlet
            pin = max(p / 3, 1.1r.ptriple)
            s = s_pT(r, pin, T_dew(r, pin) + ΔT; phase = :vapor)
            @test T_ps(r, p, s; phase = :vapor) ≈ PropsSI("T", "P", p, "S", s, f) rtol = 1e-9
        end
    end
    r = ClapeyronRefrigerant("R290")
    p, s = 15e5, s_pT(r, 3.5e5, T_dew(r, 3.5e5) + 5; phase = :vapor)
    T = T_ps(r, p, s; phase = :vapor)
    @test T ≈ T_ps(r, p, s) rtol = 1e-12                                  # Clapeyron flash
    @test T_ps(r, p, s; phase = :vapor, T0 = 330.) ≈ T rtol = 1e-12
    @test (@inferred T_ps(r, p, s; phase = :vapor)) isa Float64
    @test (@inferred h_ps(r, p, s; phase = :vapor)) isa Float64
    # exact derivatives: ∂T/∂s = T/cp; ∂T/∂p against central differences
    @test dderivative(x -> T_ps(r, p, x; phase = :vapor), s) ≈ T / cp_pT(r, p, T; phase = :vapor) rtol = 1e-10
    @test dderivative(q -> T_ps(r, q, s; phase = :vapor), p) ≈
          (T_ps(r, p + 10, s; phase = :vapor) - T_ps(r, p - 10, s; phase = :vapor)) / 20 rtol = 1e-6
    # fallbacks to the flash: wet state (s below the dew entropy) and supercritical pressure
    sw = (s_dew(r, p) + PropsSI("S", "P", p, "Q", 0, "R290")) / 2
    @test T_ps(r, p, sw; phase = :vapor) ≈ T_dew(r, p) rtol = 1e-8
    r744 = ClapeyronRefrigerant("R744")
    s744 = PropsSI("S", "P", 9e6, "T", 380., "R744")
    @test T_ps(r744, 9e6, s744; phase = :vapor) ≈ 380. rtol = 1e-6
end
