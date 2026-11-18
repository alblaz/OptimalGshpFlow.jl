# compressor_power / condenser_power: any refrigerant backend (feature/optimal-flow-components)
@testset "compressor and condenser" begin
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    spl = SPL_R290
    cla = ClapeyronRefrigerant("R290")
    p = EvaporatorParameters(refrigerant = spl, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei,
                             Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = 278.)
    pe = PropsSI("P", "T", 270., "Q", 1., "R290")
    u = [to_R(pressure_transform(spl), pe), 5000., 300., 0., 0., 0., 274., 276.]
    es = evaporator_model(u, p)
    Tc, ηis = 273.15 + 40, 0.7

    # spline backend: the original CoolProp computation
    pc = PropsSI("P", "T", Tc, "Q", 0, "R290")
    pe_ = PropsSI("P", "T", es.twophase_region.Tcl, "Q", 1, "R290")
    Tsuc = es.superheat_region.Tch
    h1 = PropsSI("H", "T", Tsuc, "P", pe_, "R290"); s1 = PropsSI("S", "T", Tsuc, "P", pe_, "R290")
    Δh = (PropsSI("H", "S", s1, "P", pc, "R290") - h1) / ηis
    ρ = PropsSI("D", "T", Tsuc, "P", pe_, "R290")
    W_ref = ρ * es.twophase_region.Vc * Δh / (1 - 0.07)
    Q1_ref = ρ * es.twophase_region.Vc * (h1 + Δh - PropsSI("H", "T", Tc, "Q", 0, "R290"))
    @test compressor_power(Tc, ηis, es, spl) ≈ W_ref rtol = 1e-13
    @test condenser_power(Tc, ηis, es, spl) ≈ Q1_ref rtol = 1e-13
    for r in (spl, cla)
        @test compressor_condenser_power(Tc, ηis, es, r) == (W = compressor_power(Tc, ηis, es, r), Q1 = condenser_power(Tc, ηis, es, r))
    end
    @test (@inferred compressor_condenser_power(Tc, ηis, es, spl)) isa NamedTuple{(:W, :Q1), Tuple{Float64, Float64}}

    # Clapeyron backend: same thermodynamics (Helmholtz EoS of the CoolProp fluid files)
    @test compressor_power(Tc, ηis, es, cla) ≈ W_ref rtol = 1e-6
    @test condenser_power(Tc, ηis, es, cla) ≈ Q1_ref rtol = 1e-6
    @test compressor_power(Tc, ηis, es, cla; ξ = 0.) ≈ W_ref * (1 - 0.07) rtol = 1e-6
    # energy balance: condenser = evaporator side + compressor work (no losses, ξ = 0)
    Q2 = es.twophase_region.Q + es.superheat_region.Q
    @test isfinite(Q2) && condenser_power(Tc, ηis, es, cla) > compressor_power(Tc, ηis, es, cla; ξ = 0.)

    # a zeotropic mixture runs through the interface (bubble pressure at Tc)
    mix = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785])
    @test isfinite(compressor_power(Tc, ηis, es, mix)) && compressor_power(Tc, ηis, es, mix) > 0
end
