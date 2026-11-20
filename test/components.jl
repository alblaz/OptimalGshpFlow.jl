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
    bhx = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = ev.Vs, λground = 3.5)
    standalone(e) = ComponentSystem((evaporator = e,);
        nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)), fixed = (T_in = 278.,))
    coupled(e; Tb = 283.15) = ComponentSystem((evaporator = e, borehole = bhx);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin), T_warm = (:evaporator => :Ts_in, :borehole => :Tout),
                 Tb = (:borehole => :Tb,)), fixed = (Tb = Tb,))
    g = (Te = ev.Te, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom, Ts_m = 276.)

    @testset "defaults and declarations" begin
        @test default_formulation(cla) == SaturationTemperature()
        @test default_formulation(spl) == NestedPressure()
        mix = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785])
        @test default_formulation(mix) == NestedPressure()
        @test default_refrigerant_pressure_drop(cla) == default_refrigerant_pressure_drop(spl) == Amalfi2016()
        @test default_refrigerant_pressure_drop(mix) == Amalfi2016()
        @test nominal_duty(cla, ev.Vref, ev.hei) ≈ ev.Qnom rtol = 0.05      # design point Te = 273 K
        for F in (NestedPressure(), SaturationTemperature(), ExplicitTv())
            e = evap(cla, F)
            @test keys(unknowns(e))[1] == (F isa NestedPressure ? :pe : :Te)
            @test length(unknowns(e)) == length(residual_scales(e)) - 1 == 6 + n_extra_unknowns(F)
            @test keys(unknowns(e)) == keys(unknowns(evap(cla, F; refrigerant_pressure_drop = nothing)))
            @test keys(ports(e)) == (:Ts_in, :Ts_out)
        end
    end

    @testset "system assembly checks" begin
        e = evap(cla)
        @test length(standalone(e)) == 7
        @test length(coupled(e)) == 8
        @test length(coupled(evap(cla, ExplicitTv()))) == 10
        # unconnected port, port in two nodes, unknown fixed node, non-square system
        @test_throws ArgumentError ComponentSystem((evaporator = e,); nodes = (T_out = (:evaporator => :Ts_out,),))
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (a = (:evaporator => :Ts_out,), b = (:evaporator => :Ts_out, :evaporator => :Ts_in)))
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)), fixed = (T_x = 1.,))
        @test_throws ArgumentError ComponentSystem((evaporator = e,);
            nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)))   # 8 unknowns, 7 residuals
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
        # without refrigerant pressure drop the component has the original equations, with the
        # two-phase effectiveness explicit (ϵe = 1 - exp(-NTU)) instead of an unknown
        e = evap(r, F; refrigerant_pressure_drop = nothing)
        nx = n_extra_unknowns(F)
        # standalone, at a point away from the solution
        sys = standalone(e)
        x = (Te = 270., Qe = 5e3, Qsup = 300., Ts_m = 276., ϵsup = 0.7, xs = 0.8)
        u = initial_guess(sys; evaporator = x, nodes = (T_out = 274.,))
        pe = OptimalGshpFlow._p_dew(r, 270.)
        U = OptimalGshpFlow.evaporator_uvalue_twophase(275., e.Vs, 5e3, pe, 0.8, sf, r, phx)[1]
        ϵe = 1 - exp(-U * phx.A * 0.8 / (sf.ρ(275.) * sf.cp(275.) * e.Vs))
        uo = evaporator_unknowns(F, r, 5.; Te = 270., Qe = 5e3, Qsup = 300., ϵe = ϵe, ϵsup = 0.7, xs = 0.8, Ts_l = 274., Ts_m = 276.)
        du, duo = zeros(7 + nx), zeros(8 + nx)
        residual!(du, u, sys)
        evaporator_model!(duo, uo, EvaporatorParameters(e, 278.), F)
        @test abs(duo[4]) < 1e-12
        @test du .* collect(residual_scales(e)) ≈ duo[[1:3; 5:8 + nx]] rtol = 1e-10
        # physical values back
        xp = physical(sys, u)
        @test xp.evaporator.Qe ≈ 5e3 rtol = 1e-14
        @test xp.nodes == (T_out = 274., T_in = 278.)
        @test port_values(sys, u, :evaporator) == (Ts_in = 278., Ts_out = 274.)
        # Duals pass through
        J = djacobian((d, v) -> residual!(d, v, sys), similar(du), u)
        @test all(isfinite, J) && size(J) == (7 + nx, 7 + nx)
    end

    @testset "coupled system: same root as evaporatorborehole_model!" begin
        e = evap(cla; refrigerant_pressure_drop = nothing)
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

    # Two-phase zone with the refrigerant pressure drop (Qiao 2025, docs/dev/glide.md)
    @testset "effectiveness with pressure drop (Qiao 2025)" begin
        @test effectiveness_counterflow(2.0, 0.0) ≈ 1 - exp(-2)
        @test effectiveness_counterflow(2.0, 1.0) ≈ (1 - exp(-2)) * (1 + 1/2) - exp(-2)     # eq. 58
        @test effectiveness_counterflow(2.0, 0.2, 1e-9) ≈ effectiveness_counterflow(2.0, 0.2) rtol = 1e-7
        @test effectiveness_counterflow(4.0, 0.17) > 1          # the refrigerant enters warmer than it leaves
        # reference: discretised zone, refrigerant temperature linear in the area, RK4 on the fluid
        # (which flows from a = 1 to a = 0, a = A/A_zone along the refrigerant)
        Trin, ΔTsat, Tfin, NTU = 268.0, -0.8, 276.0, 2.5
        Tr(a) = Trin + ΔTsat * a
        f(a, T) = NTU * (T - Tr(a))
        n = 2000; h = -1/n; Tf = Tfin; a = 1.0
        for _ in 1:n
            k1 = f(a, Tf); k2 = f(a + h/2, Tf + h/2 * k1); k3 = f(a + h/2, Tf + h/2 * k2); k4 = f(a + h, Tf + h * k3)
            Tf += h/6 * (k1 + 2k2 + 2k3 + k4); a += h
        end
        @test (Tf - Tfin) / (Trin - Tfin) ≈ effectiveness_counterflow(NTU, ΔTsat / (Trin - Tfin)) rtol = 1e-8
        # glide (φ = C_f/C_r): balanced counter-flow φ = 1 is finite and continuous (NTU/(1 + NTU)
        # for γ = 0); a blend in a brine-flow sweep crosses it
        @test effectiveness_counterflow(2.0, 0.0, 1.0) ≈ 2/3 rtol = 1e-14
        for γ in (0.0, 0.2), φ in (1 - 1e-7, 1 + 1e-7)
            @test effectiveness_counterflow(2.0, γ, φ) ≈ effectiveness_counterflow(2.0, γ, 1.0) rtol = 1e-6
        end
        @test effectiveness_counterflow(2.0, 0.0, 0.5) ≈ (1 - exp(-1)) / (1 - 0.5exp(-1)) rtol = 1e-14
        @test effectiveness_counterflow(2.0, 0.0, 2.0) ≈ (1 - exp(2)) / (1 - 2exp(2)) rtol = 1e-14
        # glide and pressure drop together: the zone integrated from the refrigerant inlet (a = 0),
        # with the fluid outlet from ε, returns the fluid inlet temperature at a = 1
        # (dTf/da = NTU (Tf − Tr), dTr/da = ΔTsat + φ NTU (Tf − Tr))
        # (names unique in the file: a local `g(y) = …` here would overwrite the guess `g` above)
        for φg in (0.6, 1.0, 1.9)
            εg = effectiveness_counterflow(NTU, ΔTsat / (Trin - Tfin), φg)
            glide_rhs(yg) = (d = NTU * (yg[1] - yg[2]); (d, ΔTsat + φg * d))
            yg = (Tfin + εg * (Trin - Tfin), Trin); hg = 1/n
            for _ in 1:n
                k1 = glide_rhs(yg); k2 = glide_rhs(yg .+ hg/2 .* k1); k3 = glide_rhs(yg .+ hg/2 .* k2)
                k4 = glide_rhs(yg .+ hg .* k3)
                yg = yg .+ hg/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
            end
            @test yg[1] ≈ Tfin atol = 1e-9
        end
    end

    # Zeotropic blend: glide in the two-phase zone. Float64 evaluations only (a solve compiles the
    # mixture code for Duals, ≈ 70 s; see examples/blend_glide.jl)
    @testset "blend glide (R454C), zone quantities" begin
        mx = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
        hei = h_bubble(mx, OptimalGshpFlow._p_bubble(mx, 313.15))
        e = Evaporator(; refrigerant = mx, secondaryfluid = sf, plateheatexchanger = phx, hei = hei,
                       Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5., refrigerant_pressure_drop = nothing)
        pe = 4.5e5
        x = (pe = pe, Qe = 5500., Qsup = 200., ϵsup = 0.5, xs = 0.5, Ts_m = 277.9)
        p = (Ts_in = 278.15, Ts_out = 275.1)
        z = evaporator_zones(e, x, p).twophase
        Te = T_dew(mx, pe)
        @test z.Te ≈ Te && z.Tr_in ≈ T_ph(mx, pe, hei)
        @test T_bubble(mx, pe) < z.Tr_in < Te - 4                     # ≈ 5 K glide over the zone
        m = ρ_pT(mx, pe, Te + 5.; phase = :vapor) * ev.Vref
        @test z.φ ≈ z.Cs * (Te - z.Tr_in) / (m * (h_dew(mx, pe) - hei)) rtol = 1e-12
        @test z.γ == 0 && z.Δp == 0
        # two-phase pinch against the refrigerant inlet: brine leaving below Te is fine above T_r,in
        pinch(Ts_out) = any(i -> occursin("two-phase pinch", i.message), validity(e, x, (Ts_in = p.Ts_in, Ts_out = Ts_out)))
        @test !pinch(Te - 1) && pinch(z.Tr_in)
        # with the refrigerant pressure drop (Amalfi): the inlet at p_in = pe + Δp; the saturation
        # shift Te − T_dew(p_in) in γ, the glide T_dew(p_in) − T_r,in in φ
        e1 = Evaporator(; refrigerant = mx, secondaryfluid = sf, plateheatexchanger = phx, hei = hei,
                        Vs = ev.Vs, Vref = ev.Vref, ΔTsup = 5.)
        z1 = evaporator_zones(e1, x, p).twophase
        p_in = pe + z1.Δp
        @test 1e3 < z1.Δp < 1e4
        @test z1.Tr_in ≈ T_ph(mx, p_in, hei) && z1.Tr_in > z.Tr_in
        @test z1.γ ≈ (Te - T_dew(mx, p_in)) / (z1.Tr_in - x.Ts_m) rtol = 1e-12
        @test z1.φ ≈ z1.Cs * (T_dew(mx, p_in) - z1.Tr_in) / (m * (h_dew(mx, pe) - hei)) rtol = 1e-12
        # blend surface tension: mole-fraction average of the components at the bubble temperature
        Tb = T_bubble(mx, pe)
        σi = [PropsSI("I", "T", Tb, "Q", 0, c) for c in ("R32", "R1234yf")]
        @test surface_tension(mx, pe) ≈ sum(mx.z .* σi) rtol = 1e-10
        # nucleate boiling of the blend (Thome & Shakir 1987): F = 1 for pure fluids; for R454C
        # F = [1 + (h/q) ΔT_bp (1 − e^{−q/(ρ_l Δh_lv β_L)})]⁻¹, → [1 + h ΔT_bp/(ρ_l Δh_lv β_L)]⁻¹ at q → 0
        @test mixture_boiling_factor(cla, pe, 2e4, 3e3) == mixture_boiling_factor(spl, pe, 2e4, 3e3) == 1
        ΔTbp, ρl, Δhlv = Te - Tb, ρ_bubble(mx, pe), h_dew(mx, pe) - h_bubble(mx, pe)
        qb, hb = 2e4, 3e3
        Fb = mixture_boiling_factor(mx, pe, qb, hb)
        @test Fb ≈ inv(1 + hb / qb * ΔTbp * (1 - exp(-qb / (ρl * Δhlv * 3e-4)))) rtol = 1e-12
        @test 0.6 < Fb < 0.9
        @test mixture_boiling_factor(mx, pe, 0.0, hb) ≈ inv(1 + hb * ΔTbp / (ρl * Δhlv * 3e-4)) rtol = 1e-12
        @test mixture_boiling_factor(mx, pe, 1e-9, hb) ≈ mixture_boiling_factor(mx, pe, 0.0, hb) rtol = 1e-9
        # the zone's U value includes it: refrigerant side h_ideal F in series with the brine side
        U, hs, href = OptimalGshpFlow.evaporator_uvalue_twophase((p.Ts_out + x.Ts_m)/2, ev.Vs, x.Qe, pe, x.xs, sf, mx, phx)
        qz = x.Qe / (phx.A * x.xs)
        prz = pe / OptimalGshpFlow.critical_pressure(mx)
        h_ideal = 1.5 * 55 * prz^0.12 / ((-log10(prz))^0.55 * (OptimalGshpFlow.molar_mass(mx) * 1000)^0.5) * qz^0.67
        @test href ≈ h_ideal * mixture_boiling_factor(mx, pe, qz, h_ideal) rtol = 1e-12
        @test U ≈ inv(1/hs + 1/href) rtol = 1e-12
    end

    @testset "evaporator with refrigerant pressure drop" begin
        e0, e1 = evap(cla; refrigerant_pressure_drop = nothing), evap(cla)
        @test e1.pressure_drop isa Amalfi2016
        S1 = @inferred coupled(e1)
        # straight-line code: a closure (map ... do) calling the local keyword function `coupled`
        # boxed it, and NonlinearSolve and Clapeyron were inferred on abstract types (25 min)
        S0 = coupled(e0)
        sol0 = solve(NonlinearProblem(S0, initial_guess(S0; evaporator = g, nodes = (T_cold = 274., T_warm = 279.))),
                     NewtonRaphson(); abstol = 1e-10)
        sol1 = solve(NonlinearProblem(S1, initial_guess(S1; evaporator = g, nodes = (T_cold = 274., T_warm = 279.))),
                     NewtonRaphson(); abstol = 1e-10)
        @test NonlinearSolve.SciMLBase.successful_retcode(sol0) && NonlinearSolve.SciMLBase.successful_retcode(sol1)
        x0, x1 = physical(S0, sol0.u), physical(S1, sol1.u)
        dp0 = evaporator_pressure_drop(e0, x0.evaporator, port_values(S0, sol0.u, :evaporator))
        dp = evaporator_pressure_drop(e1, x1.evaporator, port_values(S1, sol1.u, :evaporator))
        @test dp0.Δp_twophase == 0 && dp0.Tr_in == dp0.Te
        @test 0 < dp.Δp_twophase < 2e4 && dp.Tr_in > dp.Te && dp.γ > 0 && dp.Δp_superheat > 0
        # the refrigerant enters the two-phase zone warmer: the evaporation temperature falls a little
        @test 0 < x0.evaporator.Te - x1.evaporator.Te < 1
        # zone quantities: the closed-form profiles of both zones reach the solved far-end temperatures
        pv1 = port_values(S1, sol1.u, :evaporator)
        zn = evaporator_zones(e1, x1.evaporator, pv1)
        @test zn.twophase.Tr_in ≈ dp.Tr_in && zn.twophase.Δp ≈ dp.Δp_twophase
        Nt, ΔTs = zn.twophase.NTU, zn.twophase.Te - zn.twophase.Tr_in
        D1 = (pv1.Ts_out - zn.twophase.Tr_in - ΔTs/Nt) * exp(Nt) + ΔTs/Nt          # brine − refrigerant at xs
        @test D1 ≈ x1.evaporator.Ts_m - zn.twophase.Te atol = 1e-8
        k = zn.superheat.NTU_s - zn.superheat.NTU_r
        Ik = (x1.evaporator.Ts_m - zn.twophase.Te) * (exp(k) - 1) / k
        @test zn.twophase.Te + zn.superheat.NTU_r * Ik ≈ zn.twophase.Te + e1.ΔTsup atol = 1e-6
        @test x1.evaporator.Ts_m + zn.superheat.NTU_s * Ik ≈ pv1.Ts_in atol = 1e-6
        # the suction pressure reaches the compressor
        es = evaporator_solution(e1, x1.evaporator, port_values(S1, sol1.u, :evaporator))
        @test compressor_power(313.15, 0.7, es, cla; Δp_suction = 1e3) > compressor_power(313.15, 0.7, es, cla)
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

    # The ground temperature is a port of the borehole (a fixed node here): a new value is a new node
    # value of the same system (set_fixed), not a rebuilt component, e.g. at each time step of a
    # dynamic borehole model.
    @testset "ground temperature as a port" begin
        @test keys(ports(bhx)) == (:Tin, :Tout, :Tb)
        @test_throws ArgumentError BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = ev.Vs, λground = 3.5, Tb = 283.15)
        e = evap(cla)
        S = coupled(e)
        @test length(S) == 8                                     # a fixed node adds no unknown
        S2 = @inferred set_fixed(S; Tb = 280.15)
        @test typeof(S2) === typeof(S)                           # same type: nothing to recompile
        @test physical(S2, zeros(length(S2))).nodes.Tb == 280.15
        @test physical(S, zeros(length(S))).nodes.Tb == 283.15       # S unchanged
        u0 = initial_guess(S; evaporator = g, nodes = (T_cold = 274., T_warm = 279.))
        sol(sys) = physical(sys, solve(NonlinearProblem(sys, u0), NewtonRaphson(); abstol = 1e-10).u)
        x2 = sol(S2)
        @test x2.evaporator.Te ≈ sol(coupled(e; Tb = 280.15)).evaporator.Te rtol = 1e-12   # = rebuilt system
        @test x2.evaporator.Te < sol(S).evaporator.Te                                      # colder ground
        bs = borehole_solution(bhx, port_values(S2, u0, :borehole))
        @test bs.Tb == 280.15
        # only fixed nodes can be set
        @test_throws ArgumentError set_fixed(S; T_cold = 270.)
        @test_throws ArgumentError set_fixed(S; T_x = 270.)
    end
end

# initial_guess with more component types than inference keeps in a union: each `guess` must be
# called with concrete types (a loop over the blocks widened them to an abstract type and inferred
# the guesses, i.e. property code and Clapeyron, on `Any`: 560 s for the closed cycle of step 9d).
# Four cheap dummy component types; no runtime dispatch shows as no allocation beyond `u`.
struct GuessDummy{N} <: OptimalGshpFlow.AbstractComponent end
OptimalGshpFlow.unknowns(::GuessDummy) = (x = OptimalGshpFlow.Linear(1.0, 2.0), y = OptimalGshpFlow.Unit())
OptimalGshpFlow.ports(::GuessDummy) = (a = OptimalGshpFlow.Free(),)
OptimalGshpFlow.residual_scales(::GuessDummy) = (1.0, 1.0)
OptimalGshpFlow.residuals!(r, ::GuessDummy, x, p) = (r .= x.x - p.a; nothing)
OptimalGshpFlow.guess(::GuessDummy; x, y = 0.5) = (x = x, y = y)
guess_dummies(S) = initial_guess(S; d1 = (x = 1.0,), d2 = (x = 2.0,), d3 = (x = 3.0,), d4 = (x = 4.0, y = 0.5))

@testset "initial_guess: concrete guesses with four component types" begin
    S = ComponentSystem((d1 = GuessDummy{1}(), d2 = GuessDummy{2}(), d3 = GuessDummy{3}(), d4 = GuessDummy{4}());
                        nodes = (a = (:d1 => :a, :d2 => :a, :d3 => :a, :d4 => :a),), fixed = (a = 1.0,))
    u = guess_dummies(S)
    @test u ≈ [0.0, 0.0, 0.5, 0.0, 1.0, 0.0, 1.5, 0.0] atol = 1e-14
    @test (@allocated guess_dummies(S)) <= (@allocated zeros(length(S)))
end

# validity_issues and port_values with four component types: same widening as initial_guess (B12).
# For the closed cycle of a blend (step 9d) the post-processing never finished: `validity` and
# the callers of `port_values` (evaporator_solution, …) were inferred on `Any`, Clapeyron included.
for N in 1:4                 # one method per type, as real components: no static dispatch on a widened type
    @eval OptimalGshpFlow.validity(::GuessDummy{$N}, x, p) = ()
end
pv_dummy2(S, u) = port_values(S, u, :d2)
validity_allocs(S, u) = (validity_issues(S, u); @allocated validity_issues(S, u))   # compiled first
@testset "validity_issues and port_values: concrete with four component types" begin
    S = ComponentSystem((d1 = GuessDummy{1}(), d2 = GuessDummy{2}(), d3 = GuessDummy{3}(), d4 = GuessDummy{4}());
                        nodes = (a = (:d1 => :a, :d2 => :a, :d3 => :a, :d4 => :a),), fixed = (a = 1.0,))
    u = guess_dummies(S)
    @test isempty(validity_issues(S, u))
    @test validity_allocs(S, u) < 200                  # the issues vector only (runtime dispatch: ≈ 1000 bytes)
    @test (@inferred pv_dummy2(S, u)) == (a = 1.0,)
end
