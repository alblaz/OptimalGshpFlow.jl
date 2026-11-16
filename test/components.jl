# Component layer (feature/component-core): same equations as the original models, scaled.
using NonlinearSolve, ForwardDiff

@testset "components" begin
    sfname = "INCOMP::MEA-28%"
    sf  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), sfname)
    ev  = EvaporatorDesign(6000., 273., 273 + 40, "R290", sfname; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15)
    bh  = SingleUPipeBorehole(H = 200., D = 4., λg = 3., pipe_position = ((0.04, 0.), (-0.04, 0.)))
    cla = ClapeyronRefrigerant("R290")
    spl = SPL_R290
    evap(r, F = default_formulation(r); kw...) = Evaporator(; refrigerant = r, secondaryfluid = sf, plateheatexchanger = phx,
        hei = ev.hei, Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., formulation = F, kw...)
    bhx = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = ev.Vs, λground = 3.5, Tb = 283.15)
    standalone(e) = ComponentSystem((evaporator = e,);
        nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)), fixed = (T_in = 278.,))
    coupled(e) = ComponentSystem((evaporator = e, borehole = bhx);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin), T_warm = (:evaporator => :Ts_in, :borehole => :Tout)))
    g = (Te = ev.Te, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom, Ts_m = 276.)

    @testset "defaults and declarations" begin
        @test default_formulation(cla) == SaturationTemperature()
        @test default_formulation(spl) == NestedPressure()
        @test default_formulation(ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785])) == NestedPressure()
        @test nominal_duty(cla, ev.Vref, ev.hei) ≈ ev.Qnom rtol = 0.05      # design point Te = 273 K
        for F in (NestedPressure(), SaturationTemperature(), ExplicitTv())
            e = evap(cla, F)
            @test keys(unknowns(e))[1] == (F isa NestedPressure ? :pe : :Te)
            @test length(unknowns(e)) == length(residual_scales(e)) - 1 == 7 + n_extra_unknowns(F)
            @test keys(ports(e)) == (:Ts_in, :Ts_out)
        end
    end

    @testset "system assembly checks" begin
        e = evap(cla)
        @test length(standalone(e)) == 8
        @test length(coupled(e)) == 9
        @test length(coupled(evap(cla, ExplicitTv()))) == 11
        # unconnected port, port in two nodes, unknown fixed node, non-square system
        @test_throws ArgumentError ComponentSystem((evaporator = e,); nodes = (T_out = (:evaporator => :Ts_out,),))
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (a = (:evaporator => :Ts_out,), b = (:evaporator => :Ts_out, :evaporator => :Ts_in)))
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)), fixed = (T_x = 1.,))
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)))   # 9 unknowns, 8 residuals
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (T_out = (:evaporator => :Ts_out, :evaporator => :Tx), T_in = (:evaporator => :Ts_in,)), fixed = (T_in = 1.,))
    end

    # The system and its unknowns must be inferred from the component types: with an abstract
    # system, every caller (e.g. a continuation loop in a function) infers NonlinearSolve and
    # Clapeyron on abstract types, ≈ 5 min of compilation in the optimal-flow study.
    @testset "system and unknowns are inferred" begin
        e = evap(cla)
        S = @inferred coupled(e)
        @test isconcretetype(typeof(S))
        @test isconcretetype(typeof(@inferred standalone(e)))
        u = @inferred initial_guess(S; evaporator = g, nodes = (T_cold = 276., T_warm = 279.))
        @test length(u) == length(S)
    end

    @testset "residuals = original model residuals / scales ($(F), $(nameof(typeof(r))))" for (r, F) in
            ((cla, NestedPressure()), (cla, SaturationTemperature()), (cla, ExplicitTv()), (spl, NestedPressure()))
        e = evap(r, F)
        nx = n_extra_unknowns(F)
        # standalone, at a point away from the solution
        sys = standalone(e)
        x = (Te = 270., Qe = 5e3, Qsup = 300., Ts_m = 276., ϵe = 0.6, ϵsup = 0.7, xs = 0.8)
        u = initial_guess(sys; evaporator = x, nodes = (T_out = 274.,))
        uo = evaporator_unknowns(F, r, 5.; Te = 270., Qe = 5e3, Qsup = 300., ϵe = 0.6, ϵsup = 0.7, xs = 0.8, Ts_l = 274., Ts_m = 276.)
        du, duo = zeros(8 + nx), zeros(8 + nx)
        residual!(du, u, sys)
        evaporator_model!(duo, uo, EvaporatorParameters(e, 278.), F)
        @test du .* collect(residual_scales(e)) ≈ duo rtol = 1e-12
        # physical values back
        xp = physical(sys, u)
        @test xp.evaporator.Qe ≈ 5e3 rtol = 1e-14
        @test xp.nodes == (T_out = 274., T_in = 278.)
        @test port_values(sys, u, :evaporator) == (Ts_in = 278., Ts_out = 274.)
        # Duals pass through
        J = djacobian((d, v) -> residual!(d, v, sys), similar(du), u)
        @test all(isfinite, J) && size(J) == (8 + nx, 8 + nx)
    end

    @testset "coupled system: same root as evaporatorborehole_model!" begin
        e = evap(cla)
        sys = coupled(e)
        u0 = initial_guess(sys; evaporator = g, nodes = (T_cold = 274., T_warm = 279.))
        sol = solve(NonlinearProblem(sys, u0), NewtonRaphson(); abstol = 1e-10)
        @test NonlinearSolve.SciMLBase.successful_retcode(sol)
        x = physical(sys, sol.u)

        F = SaturationTemperature()
        pb = EvaporatorBoreholeParameters(refrigerant = cla, secondaryfluid = sf, plateheatexchanger = phx, hei = ev.hei,
                 Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., borehole = bh, λground = 3.5, Tb = 283.15)
        uo = evaporator_unknowns(F, cla, 5.; Te = ev.Te, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom, Ts_l = 274., Ts_m = 276., Ts_h = 279.)
        uso, its = newton_solve((d, u) -> evaporatorborehole_model!(d, u, pb, F), uo)   # as in property_interface.jl
        @test its > 0
        xo = map(from_R, evaporator_transforms(F, cla; borehole = true), uso)
        # the borehole component takes the fluid density at its own mean temperature (see ROADMAP)
        @test x.evaporator.Te ≈ xo[1] atol = 1e-3
        @test x.evaporator.Qe ≈ xo[2] rtol = 1e-3
        @test x.nodes.T_cold ≈ xo[7] atol = 1e-3
        @test x.nodes.T_warm ≈ xo[9] atol = 1e-3

        # post-processing through the original solution types
        es = evaporator_solution(e, x.evaporator, port_values(sys, sol.u, :evaporator))
        @test es.twophase_region.Tcl ≈ x.evaporator.Te rtol = 1e-12
        @test es.twophase_region.Q ≈ x.evaporator.Qe rtol = 1e-12
        bs = borehole_solution(bhx, port_values(sys, sol.u, :borehole))
        @test bs isa BoreholeSolution
    end

    @testset "unscaled system (Qscale = 1) has the same root" begin
        s1 = coupled(evap(cla; Qscale = 1.0)); s2 = coupled(evap(cla))
        sols = map((s1, s2)) do sys
            u0 = initial_guess(sys; evaporator = g, nodes = (T_cold = 274., T_warm = 279.))
            physical(sys, solve(NonlinearProblem(sys, u0), NewtonRaphson(); abstol = 1e-10).u)
        end
        @test sols[1].evaporator.Te ≈ sols[2].evaporator.Te atol = 1e-8
        @test sols[1].evaporator.Qe ≈ sols[2].evaporator.Qe rtol = 1e-9
    end
end
