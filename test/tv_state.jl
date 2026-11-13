# (T, v) state layer for ClapeyronRefrigerant (feature/property-interface, step 2)
using ForwardDiff

@testset "(T, v) state layer" begin
    for r in (ClapeyronRefrigerant("R290"), ClapeyronRefrigerant("R290"; transport = EntropyScalingTransport()),
              ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C"))
        p = 4.7e5
        for (T, phase) in ((T_dew(r, p) + 5, :vapor), (T_bubble(r, p) - 5, :liquid))
            v = volume_pT(r, p, T; phase = phase)
            # explicit (T, v) evaluations reproduce the (p, T) functions
            @test p_Tv(r, T, v) ≈ p rtol = 1e-12
            @test abs(pressure_residual(r, T, v, p)) < 1e-12
            @test h_Tv(r, T, v)  ≈ h_pT(r, p, T; phase = phase)  rtol = 1e-12
            @test s_Tv(r, T, v)  ≈ s_pT(r, p, T; phase = phase)  rtol = 1e-12
            @test cp_Tv(r, T, v) ≈ cp_pT(r, p, T; phase = phase) rtol = 1e-12
            @test ρ_Tv(r, T, v)  ≈ ρ_pT(r, p, T; phase = phase)  rtol = 1e-12
            if !(phase === :liquid && r.transport isa CoolPropTransport)   # (p, T) liquid via CoolProp: PT flash
                @test μ_Tv(r, T, v) ≈ μ_pT(r, p, T; phase = phase) rtol = 1e-10
                @test λ_Tv(r, T, v) ≈ λ_pT(r, p, T; phase = phase) rtol = 1e-10
            end
            # log-scale volume map: co-volume to ∞
            @test volume_from_R(r, volume_to_R(r, v)) ≈ v rtol = 1e-14
            @test volume_from_R(r, -20.) > r.vmin
            # AD in (T, v) agrees with central differences
            for f in (p_Tv, h_Tv, s_Tv, cp_Tv, μ_Tv, λ_Tv)
                g = ForwardDiff.gradient(x -> f(r, x[1], x[2]), [T, v])
                @test g[1] ≈ (f(r, T + 1e-3, v) - f(r, T - 1e-3, v)) / 2e-3 rtol = 1e-5
                @test g[2] ≈ (f(r, T, v * (1 + 1e-6)) - f(r, T, v * (1 - 1e-6))) / (2e-6 * v) rtol = 1e-5
            end
        end
    end

    # a (T, v) state closed by Newton on the pressure residual, from an ideal-gas guess
    r = ClapeyronRefrigerant("R290")
    p, T = 4.7e5, T_dew(r, 4.7e5) + 5
    u = volume_to_R(r, 8.314462618 * T / p)
    for _ in 1:20
        f = pressure_residual(r, T, volume_from_R(r, u), p)
        abs(f) < 1e-13 && break
        u -= f / ForwardDiff.derivative(x -> pressure_residual(r, T, volume_from_R(r, x), p), u)
    end
    @test volume_from_R(r, u) ≈ volume_pT(r, p, T; phase = :vapor) rtol = 1e-10

    # CoolProp transport through the cached AbstractState equals PropsSI
    @test μ_pT(r, p, T) == PropsSI("V", "P", p, "T", T, "R290")
    @test λ_dew(r, p) == PropsSI("CONDUCTIVITY", "P", p, "Q", 1, "R290")
end

@testset "explicit saturation at T (superancillaries)" begin
    for f in ("R290", "R32", "R1234yf", "R134a")
        r = ClapeyronRefrigerant(f)
        @test r.superancillary
        pure = r.eos.pures[1]
        for T in (0.6, 0.7, 0.8, 0.9, 0.97) .* r.Tcrit
            p, vl, vv = saturation_T(r, T)
            # phase equilibrium: equal pressure and Gibbs energy of both phases
            @test p_Tv(r, T, vl) ≈ p rtol = 1e-9
            @test p_Tv(r, T, vv) ≈ p rtol = 1e-9
            @test h_Tv(r, T, vv) - T * s_Tv(r, T, vv) ≈ h_Tv(r, T, vl) - T * s_Tv(r, T, vl) atol = 1e-9 * r.Tcrit * 8.314 / r.M
            # agrees with CoolProp and with the iterative solver (fallback path)
            @test p ≈ PropsSI("P", "T", T, "Q", 1, f) rtol = 1e-9
            ps, vls, vvs = OptimalGshpFlow._saturation_solve(pure, T)
            @test p ≈ ps rtol = 1e-9
            @test vl ≈ vls rtol = 1e-9
            @test vv ≈ vvs rtol = 1e-9
            # exact AD: Clausius-Clapeyron
            @test ForwardDiff.derivative(t -> p_sat(r, t), T) ≈ r.M * (h_Tv(r, T, vv) - h_Tv(r, T, vl)) / (T * (vv - vl)) rtol = 1e-9   # h per kg, v molar
        end
        @test isnan(p_sat(r, r.Tcrit + 1))
        @test r.Ttriple < 0.6 * r.Tcrit
    end
    # consistent with the pressure-based dew point
    r = ClapeyronRefrigerant("R290")
    p, _, vv = saturation_T(r, 268.15)
    @test T_dew(r, p) ≈ 268.15 rtol = 1e-10
    @test h_Tv(r, 268.15, vv) ≈ h_dew(r, p) rtol = 1e-9
    @test_throws ArgumentError saturation_T(ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785]), 270.)
end

# R744: triple point (216.6 K) above 0.6 Tc (182 K). Superancillary detection must only sample
# inside (Ttriple, Tcrit), and the dew point at a given pressure must not be NaN
# (Clapeyron 0.6.30 saturation_temperature returns NaN for R744 at every pressure).
@testset "R744 saturation" begin
    r = ClapeyronRefrigerant("R744")
    @test r.superancillary
    for p in (6e5, 2e6, 3.5e6, 6e6)
        @test T_dew(r, p) ≈ PropsSI("T", "P", p, "Q", 1, "R744") rtol = 1e-9
        @test T_bubble(r, p) ≈ PropsSI("T", "P", p, "Q", 0, "R744") rtol = 1e-9
        @test h_dew(r, p) - h_bubble(r, p) ≈ PropsSI("H", "P", p, "Q", 1, "R744") - PropsSI("H", "P", p, "Q", 0, "R744") rtol = 1e-7
        @test ForwardDiff.derivative(q -> T_dew(r, q), p) ≈ (T_dew(r, p + 1) - T_dew(r, p - 1)) / 2 rtol = 1e-6
    end
end
