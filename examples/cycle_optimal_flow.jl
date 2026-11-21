#= Optimal brine flow with the closed cycle (after step 9d).

`study_tb` of `optimalflow_components.jl` (open cycle: evaporator + borehole, condensing temperature
and liquid enthalpy from a convention) redone on the closed cycle of `examples/cycle.jl`
(`heatpump_cycle`: the condensing pressure from the water side, floor heating water in at 30 °C,
fixed water flow). Same grids (`FLOWS`, `TBS`), same machine (`refrigerant_design(r)`), same
optimum (`optimum`: parabola through the best grid point). COP1 = Q1 / (W + Ppump), COP2 = Q2 / (W + Ppump),
Ppump: brine side only (borehole and evaporator), as the open study.

Usage:
    include("examples/cycle_optimal_flow.jl")
    r290 = ClapeyronRefrigerant("R290")
    closed = study_cycle_tb(r290); report_tb(closed)
    open = study_tb(r290); compare_optimal_flow(closed, open)
(one refrigerant per call; R454C ≈ 7 min: 39 flows × 15 Tb at ≈ 0.7 s per solve)
=#
include(joinpath(@__DIR__, "cycle.jl"))

# one row of the flow sweep: the fields of the open study's rows used by `optimum` and `report_tb`
function cycle_row(S, u, Vf)
    c = cycle_point(S, u)
    bs = borehole_solution(S.components.borehole, port_values(S, u, :borehole))
    return (Vf = Vf, Vs = Vf * DESIGN.Vs, Re = reynolds(bs, BRINE_PROPS), c..., COP2 = c.Q2 / (c.W + c.Ppump))
end

"""
    study_cycle_tb(r; Tbs = TBS, flows = FLOWS, Tw_in = TW_IN, design = refrigerant_design(r),
                   condensation = LongoJH(), log = devnull) -> Vector

The flow sweep of the closed cycle for each borehole wall temperature in `Tbs` (warm to cold): a
continuation in `Tb` at the first flow gives the start of each sweep, then a continuation in the
flow. One entry per `Tb`: `(Tb, rows, COP1, COP2)` with the optima of `optimum`, as `study_tb`.
"""
function study_cycle_tb(r; Tbs = TBS, flows = FLOWS, Tw_in = TW_IN, design = refrigerant_design(r),
                        condensation = LongoJH(), log = devnull)
    hp(Vf, Tb) = cycle_system(r; Vf = Vf, Tb = Tb, Tw_in = Tw_in, design = design, condensation = condensation)
    Vf1 = float(first(flows))
    u0 = cycle_guess(hp(Vf1, first(Tbs)), first(Tbs))
    _log(log, r.name, ": continuation in Tb")
    starts = continuation(Tb -> hp(Vf1, Tb), Tbs, u0, TrustRegion(); abstol = 1e-10)
    out = []
    for s in starts
        s.ok || (@warn "$(r.name): no closed cycle at Tb = $(s.λ - 273.15) °C, $Vf1 × design flow"; break)
        _log(log, r.name, " Tb ", s.λ - 273.15)
        res = continuation(Vf -> hp(Vf, s.λ), flows, s.u, TrustRegion(); abstol = 1e-10)
        for x in res
            x.ok || @warn "$(r.name): no closed cycle at $(x.λ) × design flow, Tb = $(s.λ - 273.15) °C"
        end
        # concrete rows (comprehension): the post-processing is inferred on the system type
        rows = [cycle_row(hp(x.λ, s.λ), x.u, x.λ) for x in res if x.ok]
        push!(out, (Tb = s.λ, rows = rows, COP1 = optimum(rows, :COP1), COP2 = optimum(rows, :COP2)))
    end
    return out
end

"""
    compare_optimal_flow(closed, open)

Optimal flow and COP1 of the closed cycle (`study_cycle_tb`) against the open study (`study_tb`)
at each `Tb`, and the COP1 lost by running the closed cycle at the open-study optimum.
"""
function compare_optimal_flow(closed, open)
    @printf("%7s | %6s %7s %6s | %6s %7s | %7s %8s %9s\n", "Tb[°C]", "Vf*", "COP1*", "Re", "Vf*", "COP1*",
            "ΔVf*", "ΔCOP1*", "loss[%]")
    @printf("%7s | %22s | %14s |\n", "", "closed cycle", "open cycle")
    for c in closed
        i = findfirst(o -> o.Tb ≈ c.Tb, open)
        i === nothing && continue
        o = open[i]
        # closed-cycle COP1 at the open optimum: linear interpolation on the closed grid
        pts = sort(filter(x -> isempty(x.issues), c.rows); by = x -> x.Vf)
        vfs, cops = [x.Vf for x in pts], [x.COP1 for x in pts]
        j = clamp(searchsortedlast(vfs, o.COP1.Vf), 1, length(vfs) - 1)
        w = (o.COP1.Vf - vfs[j]) / (vfs[j + 1] - vfs[j])
        cop_at_open = (1 - w) * cops[j] + w * cops[j + 1]
        @printf("%7.1f | %6.3f %7.3f %6.0f | %6.3f %7.3f | %+7.3f %+7.1f%% %9.2f\n", c.Tb - 273.15, c.COP1.Vf,
                c.COP1.value, c.COP1.row.Re, o.COP1.Vf, o.COP1.value, c.COP1.Vf - o.COP1.Vf,
                100(o.COP1.value / c.COP1.value - 1), 100(1 - cop_at_open / c.COP1.value))
    end
    println("ΔVf* = closed − open; ΔCOP1* = open vs closed; loss: closed-cycle COP1 at the open optimum vs its own.")
end
