# Closed heat-pump cycle (step 9d): evaporator + compressor + condenser + borehole
@testset "closed cycle" begin
    r = SPL_R290
    # generic compression (inside residuals) against the CoolProp one of the post-processing
    pe, pc = PropsSI("P", "T", 273.15, "Q", 1, "R290"), PropsSI("P", "T", 313.15, "Q", 0, "R290")
    c0, c1 = OptimalGshpFlow.compression(r, pe, 276.65, 0.0, pc, 0.7), compression_state(r, pe, 276.65, 0.0, pc, 0.7)
    @test c1.h_in ≈ c0.h_in rtol = 1e-6
    @test c1.h_out ≈ c0.h_out rtol = 1e-5
    @test vapour_h_ps(r, pc, s_pT(r, pc, 330.0)) ≈ h_pT(r, pc, 330.0) rtol = 1e-10

    brine = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), "INCOMP::MEA-28%")
    water = SecondaryFluid(range(273.15 + 5, 273.15 + 70, length = 40), "Water")
    ev = EvaporatorDesign(6000., 273., 273 + 40, "R290", "INCOMP::MEA-28%"; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    cd = condenser_design(7641.0, 313.15, 303.15, 308.15, water)
    evap = Evaporator(refrigerant = r, secondaryfluid = brine, plateheatexchanger = PlateHeatExchanger(ev.A, 0.5, 0.15),
                      hei = nothing, hei_nominal = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 3.5)
    @test keys(ports(evap)) == (:Ts_in, :Ts_out, :hei, :p_suc, :T_suc, :m)
    comp = Compressor(refrigerant = r, ηis = 0.7)
    cond = Condenser(refrigerant = r, secondaryfluid = water, plateheatexchanger = PlateHeatExchanger(cd.A, 0.5, 0.15),
                     Vs = cd.Vs, Qscale = 7641.0)
    bhx = BoreholeHX(borehole = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.), (-0.04, 0.))),
                     secondaryfluid = brine, Vs = ev.Vs, λground = 3.)
    # concrete system type (mixed node transforms; abstract before B16: everything solved with it
    # was inferred on an abstract type, Clapeyron included)
    S = @inferred heatpump_cycle(evaporator = evap, compressor = comp, condenser = cond, borehole = bhx, Tb = 281.15, Tw_in = 303.15)
    @test length(S) == 19
    u0 = heatpump_guess(S; Te = 273.15, Tc = 313.15, T_cold = 275., T_warm = 278., Tw_out = 308.15)
    u, its = newton_solve((du, u) -> OptimalGshpFlow.residual!(du, u, S), u0; tol = 1e-11)
    @test its > 0
    x = physical(S, u); n = x.nodes
    pv = port_values(S, u, :compressor)
    cs = OptimalGshpFlow._compression_state(comp, pv)
    Q2, Q1 = x.evaporator.Qe + x.evaporator.Qsup, x.condenser.Qd + x.condenser.Qc
    # closed loop: energy balance, valve isenthalpic (saturated liquid), outlet state, mass flow
    @test Q2 + n.m * (cs.h_out - cs.h_in) ≈ Q1 rtol = 1e-10
    @test n.h_liq ≈ h_bubble(r, n.pc) rtol = 1e-12
    @test n.T_suc ≈ T_dew(r, n.p_suc) + 3.5 rtol = 1e-12
    @test n.m ≈ ρ_superheated(r, n.p_suc, 3.5) * ev.Vref rtol = 1e-10
    # the open-cycle post-processing at the same state (CoolProp compression, Tc = T_bubble(pc))
    es = evaporator_solution(evap, x.evaporator, port_values(S, u, :evaporator))
    oc = compressor_condenser_power(T_bubble(r, n.pc), 0.7, es, r)
    @test compressor_power(comp, pv) ≈ oc.W rtol = 1e-5
    @test Q1 ≈ oc.Q1 rtol = 1e-5
    # condensing temperature 42.9 °C with LongoJH (lower h than the Akers branch of Longo2010 at Re_eq ≈ 1650)
    @test 270 < T_dew(r, n.p_suc) < n.T_cold < n.T_warm < 281.15 && 311 < T_dew(r, n.pc) < 318 && 307 < n.Tw_out < 309
    @test isempty(validity_issues(S, u))
    @test_throws ArgumentError heatpump_cycle(evaporator = Evaporator(refrigerant = r, secondaryfluid = brine,
            plateheatexchanger = evap.plateheatexchanger, hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 3.5),
            compressor = comp, condenser = cond, borehole = bhx, Tb = 281.15, Tw_in = 303.15)
    @test_throws ArgumentError Evaporator(refrigerant = r, secondaryfluid = brine, plateheatexchanger = evap.plateheatexchanger,
                                          hei = nothing, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 3.5)
end
