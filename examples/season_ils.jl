#= Heating season with a dynamic ground: infinite line source (ROADMAP step 8b).

The heat pump of `optimalflow_components.jl` runs continuously; the borehole wall temperature `Tb`
follows the history of the heat extracted (ILS, temporal superposition, explicit coupling: `Tb` at
the start of a step from the loads before it, then the heat pump is solved at that `Tb` with a warm
start, and its extraction is the load of the step). Two brine-flow strategies:
- a fixed flow (× design);
- `:optimal`: re-optimised once a day (COP1 maximum over the flow grid `FLOWS`, valid points).

Usage:
    include("examples/season_ils.jl")
    r = Refrigerant("R290")
    fixed = season(r; flow = 1.0); opt = season(r; flow = :optimal)
    report_season(fixed); report_season(opt)
    nodp = season(r; flow = 1.0, refrigerant_pressure_drop = nothing)    # without the refrigerant Δp
=#
include(joinpath(@__DIR__, "optimalflow_components.jl"))

const GROUND = InfiniteLineSource(λ = λGROUND, α = 1e-6, rb = BOREHOLE.rb, Tg = 273.15 + 8)
const DAY = 86400.

# solve at the flow factor Vf and ground temperature Tb from u (Newton, TrustRegion as fallback)
function solve_point(S, u)
    for alg in (NewtonRaphson(), TrustRegion())
        sol = solve(NonlinearProblem(S, u), alg; abstol = 1e-9)
        SciMLBase.successful_retcode(sol) && return sol.u
    end
    return nothing
end

# flow sweep at Tb: COP1 optimum on the grid (valid points) and the solution there
function optimal_point(r, Tb; sf = BRINE_PROPS, bh = BOREHOLE, flows = FLOWS,
                       refrigerant_pressure_drop = default_refrigerant_pressure_drop(r))
    hp(Vf) = heatpump(r, Vf; sf = sf, Tb = Tb, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop)
    u0 = first_guess(hp(first(flows)), Tb)
    res = continuation(hp, flows, u0, TrustRegion(); abstol = 1e-9)
    o = optimum(postprocess(r, res, Tb, sf, bh; refrigerant_pressure_drop = refrigerant_pressure_drop), :COP1)
    x = res[findfirst(x -> x.λ == o.row.Vf, res)]
    return x.λ, x.u
end

"""
    season(r; flow = 1.0, days = 180, dt = 3600, g = GROUND, sf = BRINE_PROPS, bh = BOREHOLE,
           refrigerant_pressure_drop = default_refrigerant_pressure_drop(r)) -> Vector

Continuous operation over `days`, time step `dt` [s]; `flow` a fixed flow factor or `:optimal`
(re-optimised daily). One row per step: time, `Tb`, flow, evaporation temperature, duties, powers,
COPs and validity issues. `refrigerant_pressure_drop = nothing`: the original equations.
"""
function season(r; flow = 1.0, days = 180, dt = 3600., g = GROUND, sf = BRINE_PROPS, bh = BOREHOLE,
                refrigerant_pressure_drop = default_refrigerant_pressure_drop(r))
    hp(Vf, Tb) = heatpump(r, Vf; sf = sf, Tb = Tb, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop)
    hist = UniformLoadHistory(g, dt)          # tabulated ILS responses: O(n) per step, no expint
    rows = []
    Vf, u = flow === :optimal ? optimal_point(r, g.Tg; sf = sf, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop) :
            (float(flow), last(continuation(λ -> hp(λ, g.Tg), 2.0:-0.05:float(flow),
                first_guess(hp(2.0, g.Tg), g.Tg), TrustRegion(); abstol = 1e-9)).u)
    S = hp(Vf, g.Tg)
    for n in 0:round(Int, days * DAY / dt) - 1
        t = n * dt
        Tb = borehole_wall_temperature(hist)
        if flow === :optimal && n > 0 && t % DAY == 0
            Vf, u = optimal_point(r, Tb; sf = sf, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop)
            S = hp(Vf, Tb)
        end
        S = set_fixed(S; T_ground = Tb)
        un = solve_point(S, u)
        un === nothing && (@warn "no solution at t = $(t/DAY) d, Tb = $(Tb - 273.15) °C"; break)
        u = un
        v = evaluate(S, u, r, sf)
        push_load!(hist, v.bs.Qb / bh.H)
        push!(rows, (t = t, Tb = Tb, Vf = Vf, Te = v.es.twophase_region.Tcl, qb = v.bs.Qb / bh.H, Δp = v.Δp,
                     v.perf..., issues = unique(i.message for i in validity_issues(S, u))))
    end
    return rows
end

"""
    report_season(rows; every = 30)

Every `every` days: `Tb`, flow, evaporation temperature, extraction per metre, COP1; then the
seasonal COP1 (heat delivered over compressor + pump energy) and the steps with validity issues.
"""
function report_season(rows; every = 30)
    @printf("%6s %8s %6s %8s %8s %8s %7s\n", "day", "Tb[°C]", "Vf", "Te[°C]", "q[W/m]", "Ppump[W]", "COP1")
    for x in rows
        x.t % (every * DAY) == 0 || continue
        @printf("%6.0f %8.2f %6.2f %8.2f %8.1f %8.1f %7.3f\n", x.t / DAY, x.Tb - 273.15, x.Vf, x.Te - 273.15,
                x.qb, x.Ppump, x.COP1)
    end
    Q1 = sum(x.Q1 for x in rows); E = sum(x.W + x.Ppump for x in rows)
    nbad = count(x -> !isempty(x.issues), rows)
    @printf("seasonal COP1 = %.3f, heat delivered %.0f kWh, final Tb = %.2f °C, %d steps with validity issues\n",
            Q1 / E, Q1 * (rows[2].t - rows[1].t) / 3.6e6, rows[end].Tb - 273.15, nbad)
end
