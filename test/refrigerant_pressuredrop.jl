# Refrigerant-side pressure drop in a brazed plate heat exchanger (glide, step 5)
@testset "refrigerant pressure drop" begin
    r = SPL_R290
    T = 268.15
    p = PropsSI("P", "T", T, "Q", 1, "R290")
    # saturated liquid density and surface tension through the cached CoolProp state
    @test ρ_bubble(r, p) ≈ PropsSI("D", "P", p, "Q", 0, "R290") rtol = 1e-12
    @test surface_tension(r, p) ≈ PropsSI("I", "P", p, "Q", 0, "R290") rtol = 1e-12
    @test dderivative(q -> surface_tension(r, q), p) < 0          # σ falls with T (and p)

    phx = PlateHeatExchanger(1.1, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)
    ρl, ρv, σ = ρ_bubble(r, p), ρ_dew(r, p), surface_tension(r, p)
    m = 0.02
    G = m / (phx.W * phx.b * (phx.np - 1)/2)
    @test channel_mass_flux(phx, m) ≈ G
    xm = 0.65
    ρm = 1/(xm/ρv + (1 - xm)/ρl)
    @test homogeneous_density(ρl, ρv, xm) ≈ ρm

    # Amalfi et al. 2016 (part 2 eqs. 7-8, part 1 eq. 18)
    dh = 2phx.b/phx.ϕ
    We = G^2*dh/(ρm*σ); Bd = (ρl - ρv)*9.81*dh^2/σ
    f = (2.125*(30/70)^9.993 + 0.955) * 15.698 * We^-0.475 * Bd^0.255 * (ρl/ρv)^-0.571
    @test friction_pressure_drop(Amalfi2016(), phx, G, xm, ρl, ρv, σ, 0.4) ≈ 2f*0.4*G^2/(dh*ρm)

    # Longo-type fit: Δp_f = ζ G²/(2ρm) on its own plate, scaled with L/(2b) on another one
    k = LONGO_FITS.r404a_vaporisation
    longo = PlateHeatExchanger(0.15, 0.278, 0.072; b = 0.002, β = 65.)
    @test friction_pressure_drop(k, longo, G, xm, ρl, ρv, σ, 0.278) ≈ 1667 * G^2/(2ρm)
    @test friction_pressure_drop(k, phx, G, xm, ρl, ρv, σ, 0.4) ≈ 1667 * (0.4/0.004)/(0.278/0.004) * G^2/(2ρm)

    # other terms
    @test acceleration_pressure_drop(G, 0.3, 0.3, ρl, ρv) == 0
    @test acceleration_pressure_drop(G, 0.3, 1.0, ρl, ρv) ≈ G^2 * 0.7 * (1/ρv - 1/ρl)
    @test port_pressure_drop(m, 500., 10.; D_port = 0.02) ≈ 0.75 * (m/(π*0.01^2))^2 * (1/1000 + 1/20)
    d = twophase_pressure_drop(Amalfi2016(), r, phx, m, p, 0.3, 1.0; L = 0.4)
    @test d.total ≈ d.friction + d.acceleration + d.gravity
    @test d.gravity ≈ homogeneous_density(ρl, ρv, 0.65) * 9.81 * 0.4
    @test twophase_pressure_drop(Amalfi2016(), r, phx, m, p, 0.3, 1.0; L = 0.4, gravity = false).gravity == 0
    @test vapour_pressure_drop(phx, m, 5., 8e-6, 0.1) ≈
          OptimalGshpFlow.f_Martin(G*dh/8e-6, 30.) * 0.1/dh * G^2/(2*5.)

    # number-type generic (Newton on the evaporator unknowns), inferred
    @test isfinite(dderivative(q -> twophase_pressure_drop(Amalfi2016(), r, phx, m, q, 0.3, 1.0; L = 0.4).total, p))
    @test isfinite(dderivative(mm -> twophase_pressure_drop(k, r, phx, mm, p, 0.3, 1.0; L = 0.4).total, m))
    @test (@inferred twophase_pressure_drop(Amalfi2016(), r, phx, m, p, 0.3, 1.0; L = 0.4)).total isa Float64

    # Clapeyron backend: same σ (CoolProp transport)
    cla = ClapeyronRefrigerant("R290")
    @test surface_tension(cla, p) ≈ surface_tension(r, p) rtol = 1e-12
end
