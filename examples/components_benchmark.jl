#= feature/component-core: effect of unknown/residual scaling on the 42-case sweep (see docs/dev/ROADMAP.md).
Needs NonlinearSolve in the environment. Usage: include this file, then
`rows = sweep("R290"); report(rows)`. =#
using OptimalGshpFlow, NonlinearSolve, CoolProp, Statistics, Printf

const SFN = "INCOMP::MEA-28%"
const SF  = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), SFN)
const VFACS = 0.2:0.1:1.5
const TSRC  = (273.15, 278.15, 283.15)
const BH = SingleUPipeBorehole(H = 200., D = 4., rpi = 0.04 * (1 - 1/17) / 2, rb = 0.115, λg = 2.,
                               pipe_position = ((0.035, 0.), (-0.035, 0.)))
const ALGS = (("Newton", NewtonRaphson()), ("TrustRegion", TrustRegion()), ("LevenbergMarquardt", LevenbergMarquardt()))

function design(fluid)
    Tcond = fluid == "R744" ? 273.15 + 25 : 273.15 + 40     # R744: subcritical stand-in for the gas cooler
    ev  = EvaporatorDesign(6000., 273., Tcond, fluid, SFN; U = 2000., ΔTs = 3., ΔT_l = 1.5)
    phx = PlateHeatExchanger(ev.A, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)
    return ev, phx, ClapeyronRefrigerant(fluid)
end

# :evap = evaporator with given inlet temperature T; :evapbh = evaporator + borehole with ground temperature T
function build(sys, scaled, r, ev, phx, Vf, T)
    e = Evaporator(; refrigerant = r, secondaryfluid = SF, plateheatexchanger = phx, hei = ev.hei, Vs = Vf * ev.Vs,
                   Vref = ev.Vref, ΔTsup = 5., (scaled ? (;) : (; Qscale = 1.0))...)
    sys === :evap && return ComponentSystem((evaporator = e,);
        nodes = (T_cold = (:evaporator => :Ts_out,), T_warm = (:evaporator => :Ts_in,)), fixed = (T_warm = T,))
    b = BoreholeHX(borehole = BH, secondaryfluid = SF, Vs = Vf * ev.Vs, λground = 3.)
    return ComponentSystem((evaporator = e, borehole = b);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin), T_warm = (:evaporator => :Ts_in, :borehole => :Tout),
                 Tb = (:borehole => :Tb,)), fixed = (Tb = T,))
end

# design state shifted to the source temperature (as `shifted_guess` in formulations_benchmark.jl)
function shifted_guess(S, sys, ev, T)
    ΔT = ev.Thh - ev.Te
    nodes = sys === :evap ? (T_cold = T - ΔT + ev.Thl - ev.Te,) : (T_cold = T - ΔT + ev.Thl - ev.Te, T_warm = T)
    return initial_guess(S; evaporator = (Te = T - ΔT, Qe = 0.9ev.Qnom, Qsup = 0.1ev.Qnom, Ts_m = T - ΔT + ev.Thl - ev.Te + 2),
                         nodes = nodes)
end

function solve_case(S, u0, alg)
    t = @elapsed sol = try
        solve(NonlinearProblem(S, u0), alg; abstol = 1e-9, maxiters = 200)
    catch
        nothing
    end
    ok = sol !== nothing && SciMLBase.successful_retcode(sol) && all(isfinite, sol.u)
    Te = ok ? physical(S, sol.u).evaporator.Te : NaN
    return (ok = ok, it = sol === nothing ? -1 : sol.stats.nsteps, t = t, u = ok ? sol.u : u0, Te = Te)
end

function sweep(fluid; systems = (:evap, :evapbh))
    ev, phx, r = design(fluid)
    rows = []
    for sys in systems, scaled in (false, true), (aname, alg) in ALGS, T in TSRC
        for Vf in VFACS
            S = build(sys, scaled, r, ev, phx, Vf, T)
            res = solve_case(S, shifted_guess(S, sys, ev, T), alg)
            push!(rows, (fluid = fluid, sys = sys, scaled = scaled, alg = aname, mode = :shifted, T = T, Vf = Vf, res...))
        end
        # continuation in flow: up from the design flow, then down from it
        u1 = nothing
        for branch in (collect(filter(≥(1.0), VFACS)), collect(reverse(filter(<(1.0), VFACS))))
            S = build(sys, scaled, r, ev, phx, 1.0, T)
            u = u1 === nothing ? shifted_guess(S, sys, ev, T) : u1
            for Vf in branch
                S = build(sys, scaled, r, ev, phx, Vf, T)
                res = solve_case(S, u, alg)
                res.ok && (u = res.u)
                Vf == 1.0 && res.ok && (u1 = res.u)
                push!(rows, (fluid = fluid, sys = sys, scaled = scaled, alg = aname, mode = :cont, T = T, Vf = Vf, res...))
            end
        end
    end
    return rows
end

function report(rows)
    for key in unique((x.fluid, x.sys, x.mode) for x in rows)
        println("\n", join(key, " / "))
        for (aname, _) in ALGS, scaled in (false, true)
            R = filter(x -> (x.fluid, x.sys, x.mode) == key && x.alg == aname && x.scaled == scaled, rows)
            S = filter(x -> x.ok, R)
            @printf("  %-19s %-9s converged %2d/%2d  median iters %4s  max iters %4s  median time %6.2f ms\n",
                    aname, scaled ? "scaled" : "unscaled", length(S), length(R),
                    isempty(S) ? "-" : string(median(getfield.(S, :it))),
                    isempty(S) ? "-" : string(maximum(getfield.(S, :it))), isempty(S) ? NaN : 1e3median(getfield.(S, :t)))
        end
    end
end

# same root with and without scaling, across algorithms
function agreement(rows)
    worst = 0.0
    for a in rows, b in rows
        (a.ok && b.ok && (a.fluid, a.sys, a.T, a.Vf) == (b.fluid, b.sys, b.T, b.Vf)) || continue
        worst = max(worst, abs(a.Te - b.Te))
    end
    return worst
end
