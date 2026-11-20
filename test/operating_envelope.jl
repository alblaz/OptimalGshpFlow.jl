# feature/evaporator-operating-envelope: Te lower bound, validity reports, continuation.
using NonlinearSolve

@testset "operating envelope" begin
    sfname = "INCOMP::MEA-28%"
    @testset "secondary fluid: freezing point and table range" begin
        Tf = freezing_point(sfname)
        @test Tf ≈ 273.15 - 18.2355 atol = 0.01
        @test isnan(freezing_point("Water"))
        # NB: not `sf`: an inner @testset is a closure sharing the outer scope, so assigning the outer
        # `sf` here would box it (Core.Box) and every component built from it would be inferred as Any
        sf1 = SecondaryFluid(sfname)
        @test sf1.Tfreeze == Tf && sf1.Tmin ≈ Tf + 0.1 && sf1.Tmax ≈ 303.15
        @test sf1.ρ(sf1.Tmin) ≈ PropsSI("D", "T", sf1.Tmin, "P", 101325, sfname) rtol = 1e-10
        sf2 = SecondaryFluid(range(260., 300., length = 20), sfname)
        @test (sf2.Tmin, sf2.Tmax, sf2.Tfreeze) == (260., 300., Tf)
        @test_throws ArgumentError SecondaryFluid("Water")          # no freezing point: Tmin needed
    end

    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.), (-0.04, 0.)))
    r   = ClapeyronRefrigerant("R290")
    evap(; Vs = ev.Vs, kw...) = Evaporator(; refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
                                hei = ev.hei, Vs = Vs, Vref = ev.Vref, ΔTsup = 5., kw...)

    @testset "lower bound of the evaporation temperature" begin
        @test default_Te_min(r) ≈ T_dew(r, 0.01r.pcrit) rtol = 1e-12          # −60 °C
        r744 = ClapeyronRefrigerant("R744")
        @test isfinite(default_Te_min(r744)) && default_Te_min(r744) > r744.Ttriple   # 0.01 pcrit < ptriple
        @test default_Te_min(SPL_R290) ≈ default_Te_min(r) atol = 0.05
        e = evap()
        @test e.Te_min == default_Te_min(r)
        @test bounds(unknowns(e).Te) == (e.Te_min, r.Tcrit)
        @test from_R(unknowns(e).Te, -50.) ≥ e.Te_min
        en = evap(formulation = NestedPressure(), Te_min = 250.)
        @test bounds(unknowns(en).pe)[1] ≈ p_sat(r, 250.) rtol = 1e-10
        @test_throws ArgumentError evap(Te_min = 400.)                         # above Tcrit
    end

    @testset "validity reports" begin
        e = evap()
        sys = ComponentSystem((evaporator = e,);
            nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)), fixed = (T_in = 278.,))
        u = initial_guess(sys; evaporator = (Te = 268., Qe = 5e3, Qsup = 300., Ts_m = 276.), nodes = (T_out = 274.,))
        @test isempty(validity_issues(sys, u))
        # brine frozen and outside its table, two-phase pinch
        u = initial_guess(sys; evaporator = (Te = 250., Qe = 5e3, Qsup = 300., Ts_m = 276.), nodes = (T_out = 250.,))
        msgs = [i.message for i in validity_issues(sys, u)]
        @test "secondary fluid below its freezing point" in msgs
        @test "secondary fluid outside its property table (properties extrapolated flat)" in msgs
        @test "two-phase pinch: Ts_out - T_r,in below ΔTpinch" in msgs
        @test all(i.component == :evaporator for i in validity_issues(sys, u))
        # superheat pinch, Te at its bound, zone vanishing
        u = initial_guess(sys; evaporator = (Te = 273.1, Qe = 5e3, Qsup = 300., Ts_m = 276., xs = 1 - 1e-5), nodes = (T_out = 274.,))
        msgs = [i.message for i in validity_issues(sys, u)]
        @test "superheat pinch: Ts_in - (Te + ΔTsup) below ΔTpinch" in msgs
        @test "superheat zone almost vanished" in msgs
        u = initial_guess(sys; evaporator = (Te = e.Te_min + 1e-3, Qe = 5e3, Qsup = 300., Ts_m = 276.), nodes = (T_out = 274.,))
        @test "evaporation temperature at its lower bound Te_min" in [i.message for i in validity_issues(sys, u)]
        # borehole
        b = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = ev.Vs, λground = 3.5)
        @test length(validity(b, (;), (Tin = 250., Tout = 280., Tb = 283.15))) == 2
    end

    @testset "continuation in flow" begin
        # NB: `build` must not call the keyword closure `evap`: capturing a local keyword function in
        # another closure boxes it (Core.Box), every `evap(...)` in this block is then inferred as Any,
        # and compiling the block with an abstract refrigerant takes practically forever.
        build(λ) = ComponentSystem((evaporator = Evaporator(; refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
                                                            hei = ev.hei, Vs = λ * ev.Vs, Vref = ev.Vref, ΔTsup = 5.),
                                    borehole = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = λ * ev.Vs, λground = 3.5));
            nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin), T_warm = (:evaporator => :Ts_in, :borehole => :Tout),
                     Tb = (:borehole => :Tb,)), fixed = (Tb = 283.15,))
        u0 = initial_guess(build(1.0); evaporator = (Te = ev.Te, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom, Ts_m = 276.),
                           nodes = (T_cold = 274., T_warm = 279.))
        λs = [1.0, 0.7, 0.4, 0.2]
        res = continuation(build, λs, u0, TrustRegion(); abstol = 1e-10)
        @test length(res) == 4 && all(x -> x.ok, res)
        @test [x.λ for x in res] == λs
        # concrete element type: code that post-processes the results is inferred (an abstract
        # `NamedTuple[]` made a study with ClapeyronRefrigerant compile practically forever)
        @test isconcretetype(eltype(res)) && eltype(res).parameters[2] == Tuple{Float64, Bool, typeof(u0), Int}
        # same roots as direct solves started from the continuation results
        for x in res
            s = solve(NonlinearProblem(build(x.λ), x.u), NewtonRaphson(); abstol = 1e-10)
            @test physical(build(x.λ), s.u).evaporator.Te ≈ physical(build(x.λ), x.u).evaporator.Te atol = 1e-7
        end
        # a failing first point stops the continuation
        bad = continuation(build, λs, fill(NaN, length(u0)), NewtonRaphson(); maxiters = 5)
        @test length(bad) == 1 && !bad[1].ok
    end
end
