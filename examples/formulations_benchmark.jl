#= Sub-step 4 of feature/property-interface: nested vs explicit (T, v) evaporator formulations on the
42-case sweep (see docs/dev/ROADMAP.md). Needs NonlinearSolve in the environment, like the other examples.
Usage: include this file, then e.g. `rows, ctx = sweep2("R290"); report(rows); costs(rows, ctx)`. =#
# Sub-step 4: nested vs explicit (T, v) formulation of the evaporator, on the 42-case sweep.
using OptimalGshpFlow, NonlinearSolve, CoolProp, ForwardDiff, LinearAlgebra, Statistics, Printf

const SFN = "INCOMP::MEA-28%"
const SF  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), SFN)
const VFACS = 0.2:0.1:1.5
const TSRC  = (273.15, 278.15, 283.15)
const FORMS = (NestedPressure(), SaturationTemperature(), ExplicitTv())
const BH = SingleUPipeBorehole(H = 200., D = 4., rpi = 0.04 * (1 - 1/17) / 2, rb = 0.115, λg = 2.,
                               pipe_position = ((0.035, 0.), (-0.035, 0.)))

function setup(fluid)
    Tcond = fluid == "R744" ? 273.15 + 25 : 273.15 + 40     # R744: subcritical stand-in for the gas cooler
    ev  = EvaporatorDesign(6000., 273., Tcond, fluid, SFN; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)
    return ev, phx, ClapeyronRefrigerant(fluid)
end

params(sys, r, ev, phx, Vf, T) = sys === :evap ?
    EvaporatorParameters(refrigerant = r, secondaryfluid = SF, plateheatexchanger = phx, hei = ev.hei,
                         Vs = Vf * ev.Vs, Vref = ev.Vref, ΔTsup = 5., Ts_h = T) :
    EvaporatorBoreholeParameters(refrigerant = r, secondaryfluid = SF, plateheatexchanger = phx, hei = ev.hei,
                         Vs = Vf * ev.Vs, Vref = ev.Vref, ΔTsup = 5., borehole = BH, λground = 3., Tb = T)

nbase(sys) = sys === :evap ? 8 : 9
model!(sys) = sys === :evap ? evaporator_model! : evaporatorborehole_model!

fixed_guess(F, sys, r, ev, T) = evaporator_unknowns(F, r, 5.; Te = ev.Te, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom,
                                                    Ts_l = T - 4, Ts_m = T - 2, Ts_h = sys === :evap ? nothing : T)

function solve_case(F, sys, p, u0, alg)
    f! = (du, u, p) -> model!(sys)(du, u, p, F)
    prob = NonlinearProblem(f!, u0, p)
    t = @elapsed sol = try
        solve(prob, alg; abstol = 1e-8, reltol = 1e-12, maxiters = 200)
    catch
        nothing
    end
    ok = sol !== nothing && SciMLBase.successful_retcode(sol) && all(isfinite, sol.u)
    return (ok = ok, it = sol === nothing ? -1 : sol.stats.nsteps, t = t, u = ok ? sol.u : u0)
end

# physical summary of a solution, independent of the formulation
function summary(F, sys, r, p, u)
    un = nested_unknowns(F, r, u, nbase(sys))
    s = sys === :evap ? evaporator_model(un, p) : evaporatorborehole_model(un, p).evaporator
    return (Te = s.twophase_region.Tcl, Q = s.twophase_region.Q + s.superheat_region.Q, xs = s.xs)
end

# condition number of the Jacobian, raw and after row equilibration (rows scaled to unit max norm)
function conditioning(F, sys, p, u)
    n = length(u)
    J = ForwardDiff.jacobian((du, x) -> model!(sys)(du, x, p, F), zeros(n), u)
    Js = J ./ maximum(abs, J; dims = 2)
    return cond(J), cond(Js)
end

function sweep(fluid; algs = (("Newton", NewtonRaphson(; autodiff = AutoForwardDiff())),
                              ("TrustRegion", TrustRegion(; autodiff = AutoForwardDiff()))),
               systems = (:evap, :evapbh))
    ev, phx, r = setup(fluid)
    rows = []
    for sys in systems, F in FORMS, (aname, alg) in algs, T in TSRC
        # fixed guess
        for Vf in VFACS
            p = params(sys, r, ev, phx, Vf, T)
            res = solve_case(F, sys, p, fixed_guess(F, sys, r, ev, T), alg)
            push!(rows, (fluid = fluid, sys = sys, F = F, alg = aname, mode = :fixed, T = T, Vf = Vf, res...,
                         phys = res.ok ? summary(F, sys, r, p, res.u) : nothing))
        end
        # continuation: from the design flow up to 1.5, then from the design flow down to 0.2
        for branch in (filter(≥(1.0), VFACS), reverse(filter(≤(1.0), VFACS)))
            u = fixed_guess(F, sys, r, ev, T)
            for Vf in branch
                Vf == 1.0 && branch[1] != 1.0 && continue
                p = params(sys, r, ev, phx, Vf, T)
                res = solve_case(F, sys, p, u, alg)
                res.ok && (u = res.u)
                push!(rows, (fluid = fluid, sys = sys, F = F, alg = aname, mode = :cont, T = T, Vf = Vf, res...,
                             phys = res.ok ? summary(F, sys, r, p, res.u) : nothing))
            end
        end
    end
    return rows, (ev, phx, r)
end

function report(rows)
    for key in unique((x.fluid, x.sys, x.mode, x.alg) for x in rows)
        println("\n", join(key, " / "))
        for F in FORMS
            R = filter(x -> (x.fluid, x.sys, x.mode, x.alg) == key && x.F == F, rows)
            S = filter(x -> x.ok, R)
            @printf("  %-24s converged %2d/%2d  median iters %4s  max iters %4s  median time %7.2f ms\n",
                    string(F), length(S), length(R), isempty(S) ? "-" : string(median(getfield.(S, :it))),
                    isempty(S) ? "-" : string(maximum(getfield.(S, :it))), isempty(S) ? NaN : 1e3median(getfield.(S, :t)))
        end
    end
end

# same root across formulations?
function agreement(rows)
    worst = (dTe = 0.0, dQ = 0.0)
    for a in rows, b in rows
        (a.ok && b.ok && a.F isa NestedPressure && !(b.F isa NestedPressure)) || continue
        (a.fluid, a.sys, a.T, a.Vf, a.mode, a.alg) == (b.fluid, b.sys, b.T, b.Vf, b.mode, b.alg) || continue
        worst = (dTe = max(worst.dTe, abs(a.phys.Te - b.phys.Te)), dQ = max(worst.dQ, abs(a.phys.Q - b.phys.Q) / a.phys.Q))
    end
    return worst
end

# ---- second pass: guess = design state shifted to the source temperature ----
shifted_guess(F, sys, r, ev, T) = (ΔT = ev.Thh - ev.Te; evaporator_unknowns(F, r, 5.; Te = T - ΔT, Qe = 0.9ev.Qnom,
    Qsup = 0.1ev.Qnom, Ts_l = T - ΔT + ev.Thl - ev.Te, Ts_m = T - ΔT + ev.Thl - ev.Te + 2, Ts_h = sys === :evap ? nothing : T))

function sweep2(fluid; algs = (("Newton", NewtonRaphson(; autodiff = AutoForwardDiff())),
                               ("TrustRegion", TrustRegion(; autodiff = AutoForwardDiff()))), systems = (:evap, :evapbh))
    ev, phx, r = setup(fluid)
    rows = []
    for sys in systems, F in FORMS, (aname, alg) in algs, T in TSRC
        for Vf in VFACS
            p = params(sys, r, ev, phx, Vf, T)
            res = solve_case(F, sys, p, shifted_guess(F, sys, r, ev, T), alg)
            push!(rows, (fluid = fluid, sys = sys, F = F, alg = aname, mode = :shifted, T = T, Vf = Vf, res...,
                         phys = res.ok ? summary(F, sys, r, p, res.u) : nothing))
        end
        u1 = nothing
        for branch in (collect(filter(≥(1.0), VFACS)), collect(reverse(filter(<(1.0), VFACS))))
            u = u1 === nothing ? shifted_guess(F, sys, r, ev, T) : u1
            for Vf in branch
                p = params(sys, r, ev, phx, Vf, T)
                res = solve_case(F, sys, p, u, alg)
                res.ok && (u = res.u)
                Vf == 1.0 && res.ok && (u1 = res.u)
                push!(rows, (fluid = fluid, sys = sys, F = F, alg = aname, mode = :cont, T = T, Vf = Vf, res...,
                             phys = res.ok ? summary(F, sys, r, p, res.u) : nothing))
            end
        end
    end
    return rows, (ev, phx, r)
end

# residual / Jacobian cost and conditioning at the design-flow solutions
function costs(rows, ctx)
    ev, phx, r = ctx
    tm(f, n) = (f(); t = time_ns(); for _ in 1:n; f(); end; (time_ns() - t) / n / 1e3)
    for sys in (:evap, :evapbh), F in FORMS
        R = filter(x -> x.sys == sys && x.F == F && x.mode == :cont && x.alg == "Newton" && x.ok && x.Vf == 1.0, rows)
        isempty(R) && continue
        x = R[1]; p = params(sys, r, ev, phx, x.Vf, x.T); u = x.u; du = zeros(length(u))
        f! = (d, v) -> model!(sys)(d, v, p, F)
        tr = tm(() -> f!(du, u), 200); tj = tm(() -> ForwardDiff.jacobian(f!, du, u), 50)
        c = [conditioning(F, sys, params(sys, r, ev, phx, y.Vf, y.T), y.u) for y in
             filter(y -> y.sys == sys && y.F == F && y.mode == :cont && y.alg == "Newton" && y.ok, rows)]
        @printf("  %-7s %-24s n=%2d  residual %6.1f µs  Jacobian %7.1f µs  cond(J) median %.1e  cond(row-scaled J) median %.1e max %.1e\n",
                string(sys), string(F), length(u), tr, tj, median(first.(c)), median(last.(c)), maximum(last.(c)))
    end
end
