# Refrigerant-side condensation in a brazed plate heat exchanger (Longo 2010; step 9b)
@testset "condensation heat transfer (Longo 2010)" begin
    # Nusselt film condensation written with the heat flux: same as the classical ΔT form
    λl, ρl, μl, Δh, L = 0.09, 480.0, 9e-5, 3.3e5, 0.278
    q = 1.2e4
    h = nusselt_film_htc(λl, ρl, μl, Δh, L, q)
    @test h ≈ 0.943 * (λl^3 * ρl^2 * 9.81 * Δh / (μl * L * (q / h)))^0.25 rtol = 1e-12
    @test nusselt_film_htc(λl, ρl, μl, Δh, L, 2q) < h                          # thicker film

    # Akers: the closed-form mean over the quality equals the numerical mean of the local value
    G, ρv, cpl, dh = 15.0, 20.0, 2700.0, 0.004
    local_h(x) = 5.03 * λl / dh * cbrt(reynolds_equivalent(G, x, ρl, ρv, μl, dh)) * cbrt(cpl * μl / λl)
    n = 2000
    @test akers_htc(G, 1.0, 0.0, ρl, ρv, μl, λl, cpl, dh) ≈ sum(local_h((i - 0.5)/n) for i in 1:n) / n rtol = 1e-6
    @test akers_htc(G, 0.4, 0.4, ρl, ρv, μl, λl, cpl, dh) ≈ local_h(0.4) rtol = 1e-12

    # Longo's BPHE (Table 1: L 278 mm, W 72 mm, b 2 mm, Φ 1.24, 65°, 4 refrigerant channels, 8
    # effective plates), R290 at 30 °C, x 1 → 0; Longo refers h and q to the projected area (× Φ).
    # Fig. 4: ≈ 1750 W/(m² K) up to G ≈ 15 (gravity), ≈ 2700 at G = 28 (Akers over-predicts, Fig. 8)
    r = SPL_R290
    Φ, Wp, b = 1.24, 0.072, 0.002
    S = 8 * 0.278 * Wp
    phx = OptimalGshpFlow.PlateHeatExchanger(Φ * S, 0.278, Wp, 4b, 9, Φ, b, 65.)
    p = PropsSI("P", "T", 303.15, "Q", 0, "R290")
    Δh_lv = h_dew(r, p) - h_bubble(r, p)
    longo(Gr) = (m = Gr * 4 * Wp * b; condensation_htc(Longo2010(), r, phx, m, p, 1.0, 0.0, m * Δh_lv / S / Φ))
    c6, c28 = longo(6.0), longo(28.0)
    @test c6.Re_eq < 1500 && c6.h == c6.gravity                                 # gravity controlled
    @test c28.Re_eq > 1800 && c28.h == c28.forced                               # forced convection
    @test Φ * c6.h ≈ 1750 rtol = 0.15
    @test Φ * c28.h ≈ 2700 rtol = 0.20
    c15 = longo(15.0)                                                           # in the C¹ blend
    @test min(c15.gravity, c15.forced) ≤ c15.h ≤ max(c15.gravity, c15.forced)
    @test dderivative(q -> condensation_htc(Longo2010(), r, phx, 6 * 4 * Wp * b, p, 1.0, 0.0, q).h, 1e4) < 0

    # blends: Silver–Bell–Ghaly lowers the film coefficient (vapour-side sensible resistance)
    mx = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
    pm = OptimalGshpFlow._p_bubble(mx, 303.15)
    cm = condensation_htc(Longo2010(), mx, phx, 0.01, pm, 1.0, 0.0, 1e4)
    @test cm.h < cm.gravity && cm.Re_eq < 1500
end
