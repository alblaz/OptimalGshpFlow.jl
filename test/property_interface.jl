# Refrigerant property interface (feature/property-interface, step 1).
#
# The legacy residual functions below are verbatim copies of `evaporator_model!` and
# `evaporatorborehole_model!` on main at bfb6816 (before the port to the interface), renamed and
# with module-qualified internals, and with the pressure mapped by `pressure_transform`
# (feature/bounded-variables: originally a linear sigmoid on `(pr_min, pcrit)`, see B4). The ported
# models must reproduce them bit for bit with the spline backend.
using ForwardDiff
import OptimalGshpFlow: evaporator_uvalue_twophase, evaporator_uvalue_superheat, uniform_Tb_coeffs
import OptimalGshpFlow.Parameters: @unpack

function legacy_evaporator_model!(du, u, p)

    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, Ts_h = p
    A = plateheatexchanger.A
    ## u is a vector that contains the unknonws of the problem    
    ## u = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m]

    # first we compute model variables that are considered locally "constant": 
    # These are Cs_e, Ue, Cs_sup, Cref_sup, Cr,  U_sup,

    # two-phase
    pe = from_R(pressure_transform(refrigerant), u[1])
    xs = sigmoid_R_to_01(u[6])
    Qe = u[2]
    Tsm_e = (u[7] + u[8])/2

    Ue, hs_e, href_e  =  evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    Cs_e = Cmin_e = secondaryfluid.ρ(Tsm_e )* secondaryfluid.cp(Tsm_e) * Vs   #  capacity of secondary fluid 

    # superheat
    Tsm_sup = (Ts_h + u[8])/2   
    Usup, hs_sup, href_sup = evaporator_uvalue_superheat(Tsm_sup, pe,  Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    ρs_sup = secondaryfluid.ρ(Tsm_sup)      
    cps_sup = secondaryfluid.cp(Tsm_sup) 
    Cs_sup = ρs_sup * cps_sup * Vs
    hsup =  refrigerant.superheated_vapor.h(log10(pe),ΔTsup)     # fixed super heating of 5 degrees (find H for that super heating degree)
    ρref_sup  = refrigerant.superheated_vapor.ρ(log10(pe), ΔTsup)
    hes = refrigerant.saturated_vapor.h(log10(pe))
    Cref_sup = ρref_sup*Vref*(hsup - hes)/ΔTsup   #Cref_sup (mean capacity of the refrigerant in the superheating region)
    Cmin_sup = min(Cref_sup,Cs_sup)
    Cmax_sup = max(Cref_sup,Cs_sup)
    Cr = Cmin_sup/Cmax_sup    

    # Heat balance equations (two-phase)
    du[1] = Cs_e * (u[8] - u[7]) - u[2]  # heat balance secondary fluid
    du[2] = refrigerant.superheated_vapor.ρ(log10(from_R(pressure_transform(refrigerant), u[1])), ΔTsup) * Vref *  (refrigerant.saturated_vapor.h(log10(from_R(pressure_transform(refrigerant), u[1])))  - hei) - u[2]   # heat balance refrigerant
    du[3] = sigmoid_R_to_01(u[4]) *Cmin_e*(u[8] - refrigerant.saturated_vapor.T(log10(from_R(pressure_transform(refrigerant), u[1]))) ) - u[2]        # ϵ-NTU method 
    du[4] =  1 - exp(-Ue * A * sigmoid_R_to_01(u[6]) / Cmin_e)  -  sigmoid_R_to_01(u[4])

    # heat balance equations (superheater)
    du[5] = Cs_sup* (Ts_h - u[8]) - u[3]  # heat balance secondary fluid
    du[6] = Cref_sup * ΔTsup - u[3]
    du[7] = sigmoid_R_to_01(u[5]) * Cmin_sup * (Ts_h - refrigerant.saturated_vapor.T(log10(from_R(pressure_transform(refrigerant), u[1])))  ) - u[3]     # ϵ-NTU 
    du[8] = (1 - exp((-Usup*A*(1- sigmoid_R_to_01(u[6]))/Cmin_sup*(1-Cr))))/(1 - Cr *(exp(-Usup*A*(1- sigmoid_R_to_01(u[6]) )/Cmin_sup*(1-Cr)))) - sigmoid_R_to_01(u[5])  

    return nothing
end

function legacy_evaporatorborehole_model!(du, u, p)

    @unpack refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, borehole, λground, Tb = p
    A = plateheatexchanger.A
    ## u is a vector that contains the unknonws of the problem    
    ## u = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, Tsf_h]

    # first we compute model variables that are considered locally "constant": 
    # These are Cs_e, Ue, Cs_sup, Cref_sup, Cr,  U_sup,

    # two-phase
    pe = from_R(pressure_transform(refrigerant), u[1])
    xs = sigmoid_R_to_01(u[6])
    Qe = u[2]
    Tsm_e = (u[7] + u[8])/2

    Ue, hs_e, href_e  =  evaporator_uvalue_twophase(Tsm_e, Vs, Qe, pe, xs, secondaryfluid, refrigerant, plateheatexchanger)
    Cs_e = Cmin_e = secondaryfluid.ρ(Tsm_e )* secondaryfluid.cp(Tsm_e) * Vs   #  capacity of secondary fluid 

    # superheat
    Tsm_sup = (u[9] + u[8])/2   
    Usup, hs_sup, href_sup = evaporator_uvalue_superheat(Tsm_sup, pe,  Vs, Vref, ΔTsup, secondaryfluid, refrigerant, plateheatexchanger)
    ρs_sup = secondaryfluid.ρ(Tsm_sup)      
    cps_sup = secondaryfluid.cp(Tsm_sup) 
    Cs_sup = ρs_sup * cps_sup * Vs
    hsup =  refrigerant.superheated_vapor.h(log10(pe),ΔTsup)     # fixed super heating of 5 degrees (find H for that super heating degree)
    ρref_sup  = refrigerant.superheated_vapor.ρ(log10(pe), ΔTsup)
    hes = refrigerant.saturated_vapor.h(log10(pe))
    Cref_sup = ρref_sup*Vref*(hsup - hes)/ΔTsup   #Cref_sup (mean capacity of the refrigerant in the superheating region)
    Cmin_sup = min(Cref_sup,Cs_sup)
    Cmax_sup = max(Cref_sup,Cs_sup)
    Cr = Cmin_sup/Cmax_sup    

    # Heat balance equations (two-phase)
    du[1] = Cs_e * (u[8] - u[7]) - u[2]  # heat balance secondary fluid
    du[2] = refrigerant.superheated_vapor.ρ(log10(from_R(pressure_transform(refrigerant), u[1])), ΔTsup) * Vref *  (refrigerant.saturated_vapor.h(log10(from_R(pressure_transform(refrigerant), u[1])))  - hei) - u[2]   # heat balance refrigerant
    du[3] = sigmoid_R_to_01(u[4]) *Cmin_e*(u[8] - refrigerant.saturated_vapor.T(log10(from_R(pressure_transform(refrigerant), u[1]))) ) - u[2]        # ϵ-NTU method 
    du[4] =  1 - exp(-Ue * A * sigmoid_R_to_01(u[6]) / Cmin_e)  -  sigmoid_R_to_01(u[4])

    # heat balance equations (superheater)
    du[5] = Cs_sup* (u[9] - u[8]) - u[3]  # heat balance secondary fluid
    du[6] = Cref_sup * ΔTsup - u[3]
    du[7] = sigmoid_R_to_01(u[5]) * Cmin_sup * (u[9] - refrigerant.saturated_vapor.T(log10(from_R(pressure_transform(refrigerant), u[1])))  ) - u[3]     # ϵ-NTU 
    du[8] = (1 - exp((-Usup*A*(1- sigmoid_R_to_01(u[6]))/Cmin_sup*(1-Cr))))/(1 - Cr *(exp(-Usup*A*(1- sigmoid_R_to_01(u[6]) )/Cmin_sup*(1-Cr)))) - sigmoid_R_to_01(u[5])  #1/(1 + exp(-u[5])) # sigmoid_R_to_01(u[5])

    # borehole model 
    # kin = uniform_Tb_coeffs(bh, λground, secondaryfluid.ρ(Tsm_e )*Vs, (u[2] + u[9])/2, secondaryflui)
    # du[9] = kin*u[2] - 1*u[9] + (1-kin)* Tb

    kin = uniform_Tb_coeffs(borehole, λground, secondaryfluid.ρ(Tsm_e)*Vs, (u[7] + u[9])/2, secondaryfluid)
    du[9] = kin*u[7] - u[9] + (1 - kin)*Tb

    return nothing
end

# plain Newton with a ForwardDiff Jacobian (NonlinearSolve is not a test dependency)
function newton_solve(f!, u; maxit = 30, tol = 1e-8)
    du = similar(u)
    for it in 1:maxit
        J = ForwardDiff.jacobian(f!, du, u)
        f!(du, u)
        maximum(abs, du) < tol && return u, it
        u = u - J \ du
    end
    return u, -1
end

@testset "refrigerant property interface" begin
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.), (-0.04, 0.)))
    spl = Refrigerant("R290")
    cla = ClapeyronRefrigerant("R290")
    @test spl isa AbstractRefrigerant && cla isa AbstractRefrigerant

    evp(r) = EvaporatorParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
                 hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = 278.)
    ebp(r) = EvaporatorBoreholeParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
                 hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., borehole = bh, λground = 3.5, Tb = 283.15)
    pe0 = PropsSI("P", "T", ev.Te, "Q", 1., "R290")
    guess(r) = [to_R(pressure_transform(r), pe0), 0.9ev.Qnom, 0.1ev.Qnom, 0., 0., 0., 274., 276.]

    @testset "spline backend" begin
        p = 4.7e5; lp = log10(p)
        @test T_dew(spl, p) === spl.saturated_vapor.T(lp)
        @test h_superheated(spl, p, 5.) === spl.superheated_vapor.h(lp, 5.)
        @test h_pT(spl, p, T_dew(spl, p) + 5.) ≈ h_superheated(spl, p, 5.) rtol = 1e-12
        @test_throws ArgumentError h_pT(spl, p, 250.; phase = :liquid)
        @test critical_pressure(spl) == spl.pcrit && molar_mass(spl) == spl.M
        # the superheated functions of both backends agree (spline error only)
        for f in (h_superheated, ρ_superheated, cp_superheated, μ_superheated, λ_superheated)
            @test f(spl, p, 5.) ≈ f(cla, p, 5.) rtol = 1e-3
        end
    end

    @testset "models: spline residuals identical to main" begin
        u8 = guess(spl)
        for u in (u8, u8 .+ [0.1, -300., 50., 0.3, -0.4, 0.5, -0.5, 0.7])
            d1, d2 = zeros(8), zeros(8)
            evaporator_model!(d1, u, evp(spl)); legacy_evaporator_model!(d2, u, evp(spl))
            @test d1 == d2
            u9 = [u; 279.]
            d1, d2 = zeros(9), zeros(9)
            evaporatorborehole_model!(d1, u9, ebp(spl)); legacy_evaporatorborehole_model!(d2, u9, ebp(spl))
            @test d1 == d2
        end
    end

    @testset "models: Clapeyron backend solves to the spline solution" begin
        us, its = newton_solve((d, u) -> evaporator_model!(d, u, evp(spl)), guess(spl))
        @test its > 0
        uc, itc = newton_solve((d, u) -> evaporator_model!(d, u, evp(cla)), guess(cla))
        @test itc > 0
        pes = from_R(pressure_transform(spl), us[1])
        pec = from_R(pressure_transform(cla), uc[1])
        @test pec ≈ pes rtol = 1e-5
        @test uc[2:end] ≈ us[2:end] rtol = 1e-5
        ss, sc = evaporator_model(us, evp(spl)), evaporator_model(uc, evp(cla))
        @test sc.twophase_region.Tcl ≈ ss.twophase_region.Tcl atol = 1e-3
        @test sc.superheat_region.Q ≈ ss.superheat_region.Q rtol = 1e-4

        ubs, i1 = newton_solve((d, u) -> evaporatorborehole_model!(d, u, ebp(spl)), [us; 279.])
        ubc, i2 = newton_solve((d, u) -> evaporatorborehole_model!(d, u, ebp(cla)), [uc; 279.])
        @test i1 > 0 && i2 > 0
        @test ubc[2:end] ≈ ubs[2:end] rtol = 1e-5
    end
end

@testset "evaporator property formulations" begin
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.), (-0.04, 0.)))
    r   = ClapeyronRefrigerant("R290")
    pe  = EvaporatorParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
              hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = 278.)
    pb  = EvaporatorBoreholeParameters(refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
              hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., borehole = bh, λground = 3.5, Tb = 283.15)
    guess(F, Ts_h) = evaporator_unknowns(F, r, 5.; Te = ev.Te, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom,
                                         Ts_l = 274., Ts_m = 276., Ts_h = Ts_h)

    for (sys, p, nb, Ts_h, f!, post) in (("evaporator", pe, 8, nothing, evaporator_model!, evaporator_model),
                                         ("evaporator + borehole", pb, 9, 279., evaporatorborehole_model!, evaporatorborehole_model))
        sols = map((NestedPressure(), SaturationTemperature(), ExplicitTv())) do F
            u0 = guess(F, Ts_h)
            @test length(u0) == nb + n_extra_unknowns(F)
            u, it = newton_solve((d, x) -> f!(d, x, p, F), u0)
            @test it > 0
            un = nested_unknowns(F, r, u, nb)
            # the nested residuals vanish at the converted solution
            d = zeros(nb); f!(d, un, p)
            @test maximum(abs, d) < 1e-6
            un
        end
        @test sols[2] ≈ sols[1] rtol = 1e-9
        @test sols[3] ≈ sols[1] rtol = 1e-9
        @test (post(sols[3], p); true)
    end

    # evaporator_unknowns ↔ nested_unknowns round trip
    u = evaporator_unknowns(ExplicitTv(), r, 5.; Te = 270., Qe = 5e3, Qsup = 300., Ts_l = 274., Ts_m = 276.)
    @test from_R(pressure_transform(r), nested_unknowns(ExplicitTv(), r, u, 8)[1]) ≈ p_sat(r, 270.) rtol = 1e-12
    @test temperature_bounds(r) == (r.Ttriple, r.Tcrit)

    # transforms of the unknowns: evaporator_unknowns maps physical values with them
    @test pressure_bounds(r) == (r.ptriple, r.pcrit)
    @test r.ptriple ≈ PropsSI("ptriple", "R290") rtol = 1e-6
    @test pressure_transform(r) isa LogBox && temperature_transform(r) isa Box
    for F in (NestedPressure(), SaturationTemperature(), ExplicitTv()), Ts_h in (nothing, 279.)
        ts = evaporator_transforms(F, r; borehole = Ts_h !== nothing)
        u = evaporator_unknowns(F, r, 5.; Te = 270., Qe = 5e3, Qsup = 300., ϵe = 0.6, ϵsup = 0.7, xs = 0.8,
                                Ts_l = 274., Ts_m = 276., Ts_h = Ts_h)
        @test length(ts) == length(u) == 8 + (Ts_h !== nothing) + n_extra_unknowns(F)
        x = map(from_R, ts, u)
        @test x[1] ≈ (F isa NestedPressure ? p_sat(r, 270.) : 270.) rtol = 1e-12
        @test collect(x[2:8]) ≈ [5e3, 300., 0.6, 0.7, 0.8, 274., 276.] rtol = 1e-12
        F isa ExplicitTv && @test x[end] ≈ volume_pT(r, p_sat(r, 270.), 272.5; phase = :vapor) rtol = 1e-12
    end
    # the Te-based formulations need a ClapeyronRefrigerant
    spl = Refrigerant("R290")
    ps = EvaporatorParameters(refrigerant = spl, secondaryfluid = sf, plateheatexchanger = phx,
             hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = 278.)
    @test_throws MethodError evaporator_model!(zeros(8), zeros(8), ps, SaturationTemperature())
end
