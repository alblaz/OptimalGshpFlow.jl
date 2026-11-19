#= Optimal brine flow of a ground-source heat pump (evaporator + borehole), with the component API.

Port of `heatpumpborehole_optimalflowrate.jl` (same settings): sweep the brine flow, solve the
evaporator coupled to the borehole by continuation, post-process compressor, condenser and pumps,
and find the flow that maximises the COP. Each point is checked with `validity_issues` (e.g. brine
below its freezing point at low flow).

`study_tb` repeats the sweep over a range of ground temperatures `Tb` (the ground cools down
during the heating season), to see how the optimal flow moves with the operating point and the
brine viscosity.

Needs NonlinearSolve in the environment. Usage:
    include("examples/optimalflow_components.jl")
    rows = study(ClapeyronRefrigerant("R290")); report(rows)
    old  = original(); compare(rows, old)          # original API (spline Refrigerant)
    tb = study_tb(ClapeyronRefrigerant("R290")); report_tb(tb)
=#
using OptimalGshpFlow, NonlinearSolve, CoolProp, LinearAlgebra, Printf

# ---- settings (as in heatpumpborehole_optimalflowrate.jl) ----
const REFRIGERANT = "R290"
const BRINE = "INCOMP::MEA-28%"
const DESIGN = EvaporatorDesign(6000., 273., 273 + 40, REFRIGERANT, BRINE; U = 2000., ΔTs = 3., ΔT_l = 1.5)
const BRINE_PROPS = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), BRINE)   # as the original
const PHX = PlateHeatExchanger(DESIGN.A, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)
const DI = 0.04 * (1 - 1/17)                          # internal pipe diameter
const BOREHOLE = SingleUPipeBorehole(H = 200., D = 4., rpi = DI/2, rb = 0.115, λg = 2.,
                                     pipe_position = ((0.035, 0.0), (-0.035, 0.0)))
const λGROUND = 3.
const TB = 273.15 - 5                                 # undisturbed ground temperature
const ΔTSUP = 3.5
const ηPUMP = 0.3
const TCOND, ηIS = 273.15 + 40, 0.7
const FLOWS = 2:-0.05:0.1                             # × design brine flow, from high to low
const TBS = 273.15 .+ (9:-1:-5)                       # ground temperatures of `study_tb`, warm to cold

# ---- model ----
function heatpump(r, Vf; sf = BRINE_PROPS, Tb = TB, bh = BOREHOLE,
                  refrigerant_pressure_drop = default_refrigerant_pressure_drop(r))
    Vs = Vf * DESIGN.Vs
    evap = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = PHX, hei = DESIGN.hei,
                      Vs = Vs, Vref = DESIGN.Vref, ΔTsup = ΔTSUP, refrigerant_pressure_drop = refrigerant_pressure_drop)
    bhx = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = Vs, λground = λGROUND)
    return ComponentSystem((evaporator = evap, borehole = bhx);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin),
                 T_warm = (:evaporator => :Ts_in,  :borehole => :Tout),
                 T_ground = (:borehole => :Tb,)), fixed = (T_ground = Tb,))
end

# pumping power: borehole U-pipe (2H) and evaporator plates
pump_power(es::EvaporatorSolution, bs::BoreholeSolution, sf) = pumping_power(es, sf, ηPUMP) + pumping_power(bs, sf, ηPUMP)

function performance(es::EvaporatorSolution, bs::BoreholeSolution, r, sf)
    W, Q1 = compressor_condenser_power(TCOND, ηIS, es, r)
    Q2 = es.twophase_region.Q + es.superheat_region.Q
    P  = pump_power(es, bs, sf)
    return (W = W, Q1 = Q1, Q2 = Q2, Ppump = P, COP1 = Q1/(W + P), COP2 = Q2/(W + P))
end

# borehole Reynolds number at the mean brine temperature
function reynolds(bs::BoreholeSolution, sf)
    T = (bs.Tin + bs.Tout)/2
    D = 2bs.borehole.rpi
    return 4sf.ρ(T) * bs.Vs / (π * D * sf.μ(T))
end

# post-processing of the continuation results `res` (ground temperature `Tb`)
function postprocess(r, res, Tb, sf, bh = BOREHOLE)
    rows = []
    for x in res
        x.ok || (@warn "no solution at $(x.λ) × design flow, Tb = $(Tb - 273.15) °C"; continue)
        S = heatpump(r, x.λ; sf = sf, Tb = Tb, bh = bh)
        ph = physical(S, x.u)
        es = evaporator_solution(S.components.evaporator, ph.evaporator, port_values(S, x.u, :evaporator))
        bs = borehole_solution(S.components.borehole, port_values(S, x.u, :borehole))
        push!(rows, (Vf = x.λ, Vs = x.λ * DESIGN.Vs, Te = es.twophase_region.Tcl, Ts_out = ph.nodes.T_cold,
                     Ts_in = ph.nodes.T_warm, Re = reynolds(bs, sf), iterations = x.iterations,
                     performance(es, bs, r, sf)..., issues = unique(i.message for i in validity_issues(S, x.u))))
    end
    return rows
end

first_guess(S, Tb) = initial_guess(S; evaporator = (Te = Tb - 8, Qe = 0.9DESIGN.Qnom, Qsup = 0.1DESIGN.Qnom, Ts_m = Tb - 4),
                                   nodes = (T_cold = Tb - 6, T_warm = Tb - 2))

"""
    study(r; flows = FLOWS, sf = BRINE_PROPS, Tb = TB, bh = BOREHOLE, alg = TrustRegion()) -> Vector

Continuation in the brine flow, from `first(flows)` down; one row per flow with the physical
solution, the performance and the validity issues. `bh`: the borehole, e.g.
`SingleUPipeBorehole(BOREHOLE; transition = SmoothTransition())` for another heat-transfer
model in the laminar-turbulent transition.
"""
function study(r; flows = FLOWS, sf = BRINE_PROPS, Tb = TB, bh = BOREHOLE, alg = TrustRegion())
    u0 = first_guess(heatpump(r, first(flows); sf = sf, Tb = Tb, bh = bh), Tb)
    res = continuation(Vf -> heatpump(r, Vf; sf = sf, Tb = Tb, bh = bh), flows, u0, alg; abstol = 1e-9)
    return postprocess(r, res, Tb, sf, bh)
end

"""
    optimum(rows, k; valid = true) -> NamedTuple

Flow factor `Vf` that maximises `k` (`:COP1` or `:COP2`), refined by a parabola through the best
grid point and its two neighbours (no extra solve); with `valid = true`, among the points without
validity issues only. `at_bound = true` when the best point is the first or last of the grid.
"""
function optimum(rows, k; valid = true)
    pts = valid ? filter(x -> isempty(x.issues), rows) : rows
    sort!(pts; by = x -> x.Vf)
    y = [getfield(x, k) for x in pts]
    i = argmax(y)
    (i == 1 || i == length(pts)) && return (Vf = pts[i].Vf, value = y[i], row = pts[i], at_bound = true)
    x0, x1, x2 = pts[i-1].Vf, pts[i].Vf, pts[i+1].Vf
    y0, y1, y2 = y[i-1], y[i], y[i+1]
    # vertex of the parabola through the three points
    d01, d12 = (y1 - y0)/(x1 - x0), (y2 - y1)/(x2 - x1)
    a = (d12 - d01)/(x2 - x0)
    a < 0 || return (Vf = x1, value = y1, row = pts[i], at_bound = false)    # flat: no vertex
    # parabola y = y0 + d01 (x - x0) + a (x - x0)(x - x1); dy/dx = 0 at xv
    xv = (x0 + x1)/2 - d01/(2a)
    yv = y0 + d01*(xv - x0) + a*(xv - x0)*(xv - x1)
    return (Vf = xv, value = yv, row = pts[i], at_bound = false)
end

"""
    study_tb(r; Tbs = TBS, flows = FLOWS, sf = BRINE_PROPS, bh = BOREHOLE, alg = TrustRegion()) -> Vector

The flow sweep of `study` for each ground temperature in `Tbs` (from warm to cold, as over a
heating season). A continuation in `Tb` at the first flow gives the starting point of each sweep
(warm start from the previous `Tb`). One entry per `Tb`: `(Tb, rows, COP1, COP2)` with the optima
of [`optimum`](@ref).
"""
function study_tb(r; Tbs = TBS, flows = FLOWS, sf = BRINE_PROPS, bh = BOREHOLE, alg = TrustRegion())
    Vf1 = float(first(flows))
    u0 = first_guess(heatpump(r, Vf1; sf = sf, Tb = first(Tbs), bh = bh), first(Tbs))
    starts = continuation(Tb -> heatpump(r, Vf1; sf = sf, Tb = Tb, bh = bh), Tbs, u0, alg; abstol = 1e-9)
    out = []
    for s in starts
        s.ok || (@warn "no solution at Tb = $(s.λ - 273.15) °C, $Vf1 × design flow"; break)
        res = continuation(Vf -> heatpump(r, Vf; sf = sf, Tb = s.λ, bh = bh), flows, s.u, alg; abstol = 1e-9)
        rows = postprocess(r, res, s.λ, sf, bh)
        push!(out, (Tb = s.λ, rows = rows, COP1 = optimum(rows, :COP1), COP2 = optimum(rows, :COP2)))
    end
    return out
end

function report(rows)
    @printf("%6s %8s %8s %8s %8s %8s %8s %7s %7s  %s\n", "Vf", "Vs[l/s]", "Te[°C]", "Tout[°C]", "Q2[W]", "W[W]",
            "Ppump[W]", "COP1", "COP2", "issues")
    for x in rows
        @printf("%6.2f %8.3f %8.2f %8.2f %8.0f %8.0f %8.1f %7.3f %7.3f  %s\n", x.Vf, 1e3x.Vs, x.Te - 273.15,
                x.Ts_out - 273.15, x.Q2, x.W, x.Ppump, x.COP1, x.COP2, join(x.issues, "; "))
    end
    valid = filter(x -> isempty(x.issues), rows)
    for (k, name) in ((:COP1, "heating COP1"), (:COP2, "COP2"))
        best = rows[argmax([getfield(x, k) for x in rows])]
        @printf("max %s = %.3f at %.2f × design flow (%.3f l/s)%s\n", name, getfield(best, k), best.Vf, 1e3best.Vs,
                isempty(best.issues) ? "" : " — has validity issues")
        if !isempty(valid) && !isempty(best.issues)
            bv = valid[argmax([getfield(x, k) for x in valid])]
            @printf("    among points without issues: %.3f at %.2f × design flow\n", getfield(bv, k), bv.Vf)
        end
    end
end

"""
    original(; flows = FLOWS) -> Vector

The original pipeline (`evaporatorborehole_model!`, spline `Refrigerant`, hand-written unknowns and
continuation loop, unscaled TrustRegion), as in `heatpumpborehole_optimalflowrate.jl`, for comparison.
"""
function original(; flows = FLOWS)
    r = Refrigerant(REFRIGERANT)
    p0 = EvaporatorBoreholeParameters(refrigerant = r, secondaryfluid = BRINE_PROPS, plateheatexchanger = PHX,
            hei = DESIGN.hei, Vs = DESIGN.Vs, Vref = DESIGN.Vref, ΔTsup = ΔTSUP, borehole = BOREHOLE,
            λground = λGROUND, Tb = TB)
    u = [to_R(pressure_transform(r), PropsSI("P", "T", DESIGN.Te, "Q", 1., REFRIGERANT)),
         0.9DESIGN.Qnom, 0.1DESIGN.Qnom, sigmoid_01_to_R(0.3), sigmoid_01_to_R(0.3), sigmoid_01_to_R(0.3),
         TB - 6, TB - 4, TB - 2]
    rows = []
    for Vf in flows
        p = EvaporatorBoreholeParameters(p0, Vs = Vf * DESIGN.Vs)
        sol = solve(NonlinearProblem(evaporatorborehole_model!, u, p), TrustRegion(); abstol = 1e-9)
        SciMLBase.successful_retcode(sol) || (@warn "original: no solution at $Vf"; continue)
        u = sol.u
        s = evaporatorborehole_model(u, p)
        push!(rows, (Vf = Vf, Te = s.evaporator.twophase_region.Tcl, Ts_out = s.borehole.Tin,
                     performance(s.evaporator, s.borehole, r, BRINE_PROPS)...))
    end
    return rows
end

function report_tb(tb)
    @printf("%7s | %6s %8s %7s %6s %8s %7s | %6s %7s  %s\n", "Tb[°C]", "Vf*", "Vs*[l/s]", "COP1*", "Re", "Te[°C]",
            "Q2[W]", "Vf*", "COP2*", "notes")
    for x in tb
        o1, o2 = x.COP1, x.COP2
        notes = String[]
        o1.at_bound && push!(notes, "COP1 optimum at the grid bound")
        o2.at_bound && push!(notes, "COP2 optimum at the grid bound")
        length(x.rows) < length(FLOWS) && push!(notes, "$(length(FLOWS) - length(x.rows)) flows unsolved")
        @printf("%7.1f | %6.3f %8.3f %7.3f %6.0f %8.2f %7.0f | %6.3f %7.3f  %s\n", x.Tb - 273.15, o1.Vf,
                1e3o1.Vf * DESIGN.Vs, o1.value, o1.row.Re, o1.row.Te - 273.15, o1.row.Q2, o2.Vf, o2.value,
                join(notes, "; "))
    end
    println("Vf*: optimal flow (× design, parabola through the grid); Re, Te, Q2 at the best grid point.")
end

function compare(new, old)
    worst = (dTe = 0.0, dCOP1 = 0.0, dCOP2 = 0.0)
    for a in new, b in old
        a.Vf ≈ b.Vf || continue
        worst = (dTe = max(worst.dTe, abs(a.Te - b.Te)), dCOP1 = max(worst.dCOP1, abs(a.COP1 - b.COP1)/b.COP1),
                 dCOP2 = max(worst.dCOP2, abs(a.COP2 - b.COP2)/b.COP2))
    end
    return worst
end
