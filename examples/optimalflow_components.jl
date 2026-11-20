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
    c = pressure_drop_tb(ClapeyronRefrigerant("R290")); report_pressure_drop_tb(c)   # with/without refrigerant Δp
    g = pressure_drop_grid(Refrigerant("R290")); report_pressure_drop_grid(g)          # capacity vs speed × flow
    r454c = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
    report_refrigerants_tb("R290" => study_tb(ClapeyronRefrigerant("R290")), "R454C" => study_tb(r454c))  # midpoint
    report_refrigerants_tb((string(c) => study_tb(r454c; design = refrigerant_design(r454c; condenser = c))
                            for c in (:bubble, :midpoint, :dew))...)                                     # conventions
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
"""
    condenser_bubble_T(r, condenser = :midpoint; T = TCOND) -> Tc

Bubble temperature of the condenser (what `compressor_condenser_power` takes: condensing pressure
`p_bubble(Tc)`, liquid at the bubble point) for the condensing temperature `T` and the convention
`condenser` for a blend with glide: `:bubble` (bubble point at `T`), `:midpoint` (mean of bubble
and dew temperature at `T`) or `:dew` (dew point at `T`). `T` for a pure fluid.
"""
condenser_bubble_T(r, condenser = :midpoint; T = TCOND) = T
function condenser_bubble_T(r::ClapeyronRefrigerant, condenser = :midpoint; T = TCOND)
    (OptimalGshpFlow.ispure(r) || condenser === :bubble) && return T
    w = condenser === :midpoint ? 0.5 : condenser === :dew ? 1.0 :
        throw(ArgumentError("condenser: :bubble, :midpoint or :dew"))
    Tb = T                             # fixed point Tb = T − w glide(Tb); the glide varies slowly
    for _ in 1:50
        Tn = T - w * (T_dew(r, OptimalGshpFlow._p_bubble(r, Tb)) - Tb)
        abs(Tn - Tb) < 1e-9 && return Tn
        Tb = Tn
    end
    return Tb
end

"""
    refrigerant_design(r; condenser = :midpoint) -> (; hei, Vref, Tc)

Expansion-valve inlet enthalpy, compressor displacement and condenser bubble temperature `Tc`
([`condenser_bubble_T`](@ref)) for the refrigerant `r`: the design of the example (R290, `DESIGN`)
for pure fluids. For a zeotropic blend the same machine size: the same nominal two-phase duty
`DESIGN.Qnom` at the design dew temperature `DESIGN.Te`, with the liquid at the condenser bubble
point.
"""
refrigerant_design(r; condenser = :midpoint) = (hei = DESIGN.hei, Vref = DESIGN.Vref, Tc = TCOND)
function refrigerant_design(r::ClapeyronRefrigerant; condenser = :midpoint)
    OptimalGshpFlow.ispure(r) && return (hei = DESIGN.hei, Vref = DESIGN.Vref, Tc = TCOND)
    Tc = condenser_bubble_T(r, condenser)
    hei = h_bubble(r, OptimalGshpFlow._p_bubble(r, Tc))
    p = OptimalGshpFlow._p_dew(r, DESIGN.Te)
    return (hei = hei, Vref = DESIGN.Qnom / (ρ_dew(r, p) * (h_dew(r, p) - hei)), Tc = Tc)
end

# `speed`: compressor speed factor (variable-speed compressor): the displaced volume flow `Vref`
function heatpump(r, Vf; sf = BRINE_PROPS, Tb = TB, bh = BOREHOLE, speed = 1.0,
                  refrigerant_pressure_drop = default_refrigerant_pressure_drop(r), design = refrigerant_design(r))
    Vs = Vf * DESIGN.Vs
    evap = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = PHX, hei = design.hei,
                      Vs = Vs, Vref = speed * design.Vref, ΔTsup = ΔTSUP, refrigerant_pressure_drop = refrigerant_pressure_drop)
    bhx = BoreholeHX(borehole = bh, secondaryfluid = sf, Vs = Vs, λground = λGROUND)
    return ComponentSystem((evaporator = evap, borehole = bhx);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin),
                 T_warm = (:evaporator => :Ts_in,  :borehole => :Tout),
                 T_ground = (:borehole => :Tb,)), fixed = (T_ground = Tb,))
end

# pumping power: borehole U-pipe (2H) and evaporator plates
pump_power(es::EvaporatorSolution, bs::BoreholeSolution, sf) = pumping_power(es, sf, ηPUMP) + pumping_power(bs, sf, ηPUMP)

function performance(es::EvaporatorSolution, bs::BoreholeSolution, r, sf; Δp_suction = 0.0, Tc = condenser_bubble_T(r))
    W, Q1 = compressor_condenser_power(Tc, ηIS, es, r; Δp_suction = Δp_suction)
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

# one solution `u` of the system `S`: physical state, refrigerant pressure drops and performance.
# With a refrigerant pressure-drop model, the superheat-zone Δp lowers the compressor suction
# pressure; without one (`nothing`, the original equations) the suction is at the dew pressure.
function evaluate(S, u, r, sf; Tc = condenser_bubble_T(r))
    e = S.components.evaporator
    ph = physical(S, u)
    pv = port_values(S, u, :evaporator)
    es = evaporator_solution(e, ph.evaporator, pv)
    bs = borehole_solution(S.components.borehole, port_values(S, u, :borehole))
    dp = evaporator_pressure_drop(e, ph.evaporator, pv)
    Δp_suction = e.pressure_drop === nothing ? 0.0 : dp.Δp_superheat
    return (ph = ph, es = es, bs = bs, Δp = dp.Δp_twophase, Δp_sup = Δp_suction, Tr_in = dp.Tr_in,
            perf = performance(es, bs, r, sf; Δp_suction = Δp_suction, Tc = Tc))
end

# post-processing of the continuation results `res` (ground temperature `Tb`)
function postprocess(r, res, Tb, sf, bh = BOREHOLE; refrigerant_pressure_drop = default_refrigerant_pressure_drop(r),
                     design = refrigerant_design(r))
    rows = []
    for x in res
        x.ok || (@warn "no solution at $(x.λ) × design flow, Tb = $(Tb - 273.15) °C"; continue)
        S = heatpump(r, x.λ; sf = sf, Tb = Tb, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop, design = design)
        v = evaluate(S, x.u, r, sf; Tc = design.Tc)
        push!(rows, (Vf = x.λ, Vs = x.λ * DESIGN.Vs, Te = v.es.twophase_region.Tcl, Tr_in = v.Tr_in, Δp = v.Δp,
                     Δp_sup = v.Δp_sup, Ts_out = v.ph.nodes.T_cold, Ts_in = v.ph.nodes.T_warm, Re = reynolds(v.bs, sf),
                     iterations = x.iterations, v.perf..., issues = unique(i.message for i in validity_issues(S, x.u))))
    end
    return rows
end

first_guess(S, Tb) = initial_guess(S; evaporator = (Te = Tb - 8, Qe = 0.9DESIGN.Qnom, Qsup = 0.1DESIGN.Qnom, Ts_m = Tb - 4),
                                   nodes = (T_cold = Tb - 6, T_warm = Tb - 2))

"""
    study(r; flows = FLOWS, sf = BRINE_PROPS, Tb = TB, bh = BOREHOLE, alg = TrustRegion(),
          refrigerant_pressure_drop = default_refrigerant_pressure_drop(r)) -> Vector

Continuation in the brine flow, from `first(flows)` down; one row per flow with the physical
solution, the performance and the validity issues. `bh`: the borehole, e.g.
`SingleUPipeBorehole(BOREHOLE; transition = SmoothTransition())` for another heat-transfer
model in the laminar-turbulent transition. `refrigerant_pressure_drop = nothing`: the original
equations (no refrigerant Δp in the evaporator, suction at the dew pressure).
"""
function study(r; flows = FLOWS, sf = BRINE_PROPS, Tb = TB, bh = BOREHOLE, alg = TrustRegion(),
               refrigerant_pressure_drop = default_refrigerant_pressure_drop(r), design = refrigerant_design(r))
    hp(Vf) = heatpump(r, Vf; sf = sf, Tb = Tb, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop, design = design)
    u0 = first_guess(hp(first(flows)), Tb)
    res = continuation(hp, flows, u0, alg; abstol = 1e-9)
    return postprocess(r, res, Tb, sf, bh; refrigerant_pressure_drop = refrigerant_pressure_drop, design = design)
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
    study_tb(r; Tbs = TBS, flows = FLOWS, sf = BRINE_PROPS, bh = BOREHOLE, alg = TrustRegion(),
             refrigerant_pressure_drop = default_refrigerant_pressure_drop(r), design = refrigerant_design(r)) -> Vector

The flow sweep of `study` for each ground temperature in `Tbs` (from warm to cold, as over a
heating season). A continuation in `Tb` at the first flow gives the starting point of each sweep
(warm start from the previous `Tb`). One entry per `Tb`: `(Tb, rows, COP1, COP2)` with the optima
of [`optimum`](@ref).
"""
function study_tb(r; Tbs = TBS, flows = FLOWS, sf = BRINE_PROPS, bh = BOREHOLE, alg = TrustRegion(),
                  refrigerant_pressure_drop = default_refrigerant_pressure_drop(r), design = refrigerant_design(r))
    hp(Vf, Tb) = heatpump(r, Vf; sf = sf, Tb = Tb, bh = bh, refrigerant_pressure_drop = refrigerant_pressure_drop, design = design)
    Vf1 = float(first(flows))
    u0 = first_guess(hp(Vf1, first(Tbs)), first(Tbs))
    starts = continuation(Tb -> hp(Vf1, Tb), Tbs, u0, alg; abstol = 1e-9)
    out = []
    for s in starts
        s.ok || (@warn "no solution at Tb = $(s.λ - 273.15) °C, $Vf1 × design flow"; break)
        res = continuation(Vf -> hp(Vf, s.λ), flows, s.u, alg; abstol = 1e-9)
        rows = postprocess(r, res, s.λ, sf, bh; refrigerant_pressure_drop = refrigerant_pressure_drop, design = design)
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

"""
    pressure_drop_tb(r; Tbs = TBS, kwargs...) -> (without, with)

`study_tb` without the refrigerant pressure drop (`nothing`, the original equations) and with the
default model ([`default_refrigerant_pressure_drop`](@ref), Amalfi 2016 for a pure fluid).
"""
pressure_drop_tb(r; Tbs = TBS, kwargs...) =
    (without = study_tb(r; Tbs = Tbs, refrigerant_pressure_drop = nothing, kwargs...),
     with = study_tb(r; Tbs = Tbs, kwargs...))

"""
    report_pressure_drop_tb(c)

Per ground temperature, at the COP1 optimum of each case: optimal flow, COP1, evaporation (dew)
temperature and two-phase Δp (with) at the best grid point; then the change of COP1 and of the
evaporator duty Q2 (capacity) at a fixed flow (design flow, both cases at the same grid point).
"""
function report_pressure_drop_tb(c)
    @printf("%7s | %6s %7s %8s | %6s %7s %8s %8s %7s | %8s %7s %7s\n", "Tb[°C]", "Vf*", "COP1*", "Te[°C]", "Vf*", "COP1*",
            "Te[°C]", "Δp[kPa]", "ΔTsat", "ΔCOP1*", "ΔCOP1@1", "ΔQ2@1")
    for (a, b) in zip(c.without, c.with)
        oa, ob = a.COP1, b.COP1
        ra = a.rows[findfirst(x -> x.Vf ≈ 1, a.rows)]; rb = b.rows[findfirst(x -> x.Vf ≈ 1, b.rows)]
        @printf("%7.1f | %6.3f %7.3f %8.2f | %6.3f %7.3f %8.2f %8.2f %7.3f | %7.2f%% %6.2f%% %6.2f%%\n", a.Tb - 273.15,
                oa.Vf, oa.value, oa.row.Te - 273.15, ob.Vf, ob.value, ob.row.Te - 273.15, ob.row.Δp/1e3,
                ob.row.Tr_in - ob.row.Te, 100(ob.value/oa.value - 1), 100(rb.COP1/ra.COP1 - 1), 100(rb.Q2/ra.Q2 - 1))
    end
    println("left: without refrigerant Δp; right: with (Amalfi 2016). ΔTsat = T_r,in − Te (two-phase zone).")
end

const PRESSURE_DROP_MODELS = (none = nothing, amalfi = Amalfi2016(), longo = LONGO_FITS.r404a_vaporisation)

"""
    pressure_drop_grid(r; Tbs = (268.15, 282.15), speeds = (0.5, 1.0, 1.5), flows = FLOWS,
                       models = PRESSURE_DROP_MODELS) -> Dict

Capacity and COP over ground temperature × compressor speed × brine flow for each refrigerant
pressure-drop model (`nothing`, Amalfi 2016, Longo-type fit: the upper end of the Δp estimates).
`Dict((model, Tb, speed) => Dict(Vf => (Q2, COP1, Te, Δp)))`; see [`report_pressure_drop_grid`](@ref).
"""
function pressure_drop_grid(r; Tbs = (268.15, 282.15), speeds = (0.5, 1.0, 1.5), flows = FLOWS,
                            models = PRESSURE_DROP_MODELS)
    out = Dict{Tuple{Symbol, Float64, Float64}, Dict{Float64, Any}}()
    for k in keys(models), Tb in Tbs, s in speeds
        hp(Vf) = heatpump(r, Vf; Tb = Tb, speed = s, refrigerant_pressure_drop = models[k])
        res = continuation(hp, flows, first_guess(hp(first(flows)), Tb), TrustRegion(); abstol = 1e-9)
        d = Dict{Float64, Any}()
        for x in res
            x.ok || (@warn "$k: no solution at Tb = $(Tb - 273.15) °C, speed $s, $(x.λ) × design flow"; continue)
            v = evaluate(hp(x.λ), x.u, r, BRINE_PROPS)
            d[round(x.λ, digits = 2)] = (Q2 = v.perf.Q2, COP1 = v.perf.COP1, Te = v.es.twophase_region.Tcl, Δp = v.Δp)
        end
        out[(k, Tb, s)] = d
    end
    return out
end

"""
    report_pressure_drop_grid(g; Vfs = (0.25, 0.5, 1.0, 2.0))

Changes of the evaporator duty Q2 (capacity) and of COP1 with the Amalfi and Longo-type pressure
drops, relative to no refrigerant pressure drop, at the flows `Vfs` of [`pressure_drop_grid`](@ref).
"""
function report_pressure_drop_grid(g; Vfs = (0.25, 0.5, 1.0, 2.0))
    @printf("%5s %5s %5s | %7s %7s | %7s %7s %7s | %7s %7s %7s\n", "Tb", "speed", "Vf", "Q2[W]", "Te[°C]",
            "ΔpA kPa", "dQ2 A%", "dCOP A%", "ΔpL kPa", "dQ2 L%", "dCOP L%")
    for (k, Tb, s) in sort!(collect(filter(k -> k[1] === :none, keys(g))); by = k -> (k[2], k[3])), Vf in Vfs
        a = get(g[(:none, Tb, s)], Vf, nothing); b = get(g[(:amalfi, Tb, s)], Vf, nothing); l = get(g[(:longo, Tb, s)], Vf, nothing)
        (a === nothing || b === nothing || l === nothing) && continue
        @printf("%5.0f %5.1f %5.2f | %7.0f %7.2f | %7.2f %7.2f %7.2f | %7.2f %7.2f %7.2f\n", Tb - 273.15, s, Vf, a.Q2,
                a.Te - 273.15, b.Δp/1e3, 100(b.Q2/a.Q2 - 1), 100(b.COP1/a.COP1 - 1), l.Δp/1e3, 100(l.Q2/a.Q2 - 1),
                100(l.COP1/a.COP1 - 1))
    end
    println("A: Amalfi 2016; L: Longo-type fit (ζ = 1.667). Q2 and Te without refrigerant Δp.")
end

"""
    report_refrigerants_tb(tbs::Pair...)

Side by side, per ground temperature, the COP1 optimum of several `study_tb` results, e.g.
`report_refrigerants_tb("R290" => tb290, "R454C" => tb454)`: optimal flow, COP1, evaporation (dew)
temperature and evaporator duty Q2 at the best grid point.
"""
function report_refrigerants_tb(tbs::Pair...)
    @printf("%7s", "Tb[°C]")
    for (name, _) in tbs
        @printf(" | %-28s", name)
    end
    println()
    @printf("%7s", "")
    for _ in tbs
        @printf(" | %6s %6s %7s %6s", "Vf*", "COP1*", "Te[°C]", "Q2[W]")
    end
    println()
    for i in eachindex(last(first(tbs)))
        @printf("%7.1f", last(first(tbs))[i].Tb - 273.15)
        for (_, tb) in tbs
            o = tb[i].COP1
            @printf(" | %6.3f %6.3f %7.2f %6.0f", o.Vf, o.value, o.row.Te - 273.15, o.row.Q2)
        end
        println()
    end
end
