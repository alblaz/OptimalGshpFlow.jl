#= Overview figures of the closed heat-pump cycle (evaporator + compressor + condenser + borehole).

For each refrigerant, the closed cycle of `examples/cycle.jl` at the design brine flow over the
borehole wall temperatures `TBS` (continuation in Tb, warm to cold), and for a few of them
(`POINT_TBS`) the state of the cycle:
- p-h diagram: saturation dome, vapour isotherms, the cycle (compression as a straight line between
  suction and discharge: the path inside the compressor is not modelled, only its end states);
- temperature profiles along the evaporator (refrigerant, brine) and the condenser (refrigerant,
  water) against the area fraction from the refrigerant inlet (library `evaporator_profiles`,
  `condenser_profiles`);
- energy: heat from the ground (two-phase + superheat zones), electricity (compressor, brine pumps
  of the borehole and of the evaporator), heat to the water (desuperheating + condensation).
Water-side pump of the condenser not modelled (as COP1 of the study: brine pumps only).

Folders (`examples/figures/cycle/`, not tracked by git):
    <refrigerant>/Tb_+08C/  ph.png, profiles.png, energy.png, overview.png   one operating point
    <refrigerant>/vs_Tb.png                                                  heat, power, COP, temperatures vs Tb
    comparison/ph_Tb_+00C.png, comparison/vs_Tb.png                          refrigerants together

Needs NonlinearSolve and CairoMakie in the environment. Usage (one refrigerant per call: a vector of
both has an abstract element type, and the model would be inferred on it):
    include("examples/cycle_figures.jl")
    r290 = ClapeyronRefrigerant("R290")
    r454c = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
    s290 = cycle_sweep(r290); cycle_figures(s290)          # ≈ 15 s after compilation
    s454 = cycle_sweep(r454c); cycle_figures(s454)         # ≈ 1 min
    comparison_figures(s290, s454)
=#
include(joinpath(@__DIR__, "cycle.jl"))
using CairoMakie

const POINT_TBS = 273.15 .+ (8.0, 0.0, -5.0)
const FIGDIR = joinpath(@__DIR__, "figures", "cycle")
const C2K = 273.15
const PTICKS = ([1, 2, 3, 5, 10, 20, 30, 50], ["1", "2", "3", "5", "10", "20", "30", "50"])
const COLORS = (ref = :firebrick, brine = :steelblue, water = :darkorange, dome = :gray35)

# ---- data ----

"""
    cycle_snapshot(S, u; n = 60) -> NamedTuple

Everything the figures need from a closed-cycle solution: the refrigerant states of the p-h diagram,
the temperature profiles of both heat exchangers, the energies and the COP.
"""
function cycle_snapshot(S, u; n = 60)
    (; evaporator, compressor, condenser, borehole) = S.components
    r = evaporator.refrigerant
    ph = physical(S, u)
    nd = ph.nodes
    pv_e, pv_c = port_values(S, u, :evaporator), port_values(S, u, :condenser)
    ze = evaporator_zones(evaporator, ph.evaporator, pv_e)
    es = evaporator_solution(evaporator, ph.evaporator, pv_e)
    bs = borehole_solution(borehole, port_values(S, u, :borehole))
    pe, pc = nd.p_suc, nd.pc
    states = (discharge = (p = pc, h = nd.h_dis, T = vapour_T_ph(r, pc, nd.h_dis)),
              cond_dew = (p = pc, h = h_dew(r, pc), T = T_dew(r, pc)),
              liquid = (p = pc, h = nd.h_liq, T = T_bubble(r, pc)),
              evap_in = (p = pe + ze.twophase.Δp, h = nd.h_liq, T = ze.twophase.Tr_in),
              evap_dew = (p = pe, h = h_dew(r, pe), T = T_dew(r, pe)),
              suction = (p = pe, h = h_pT(r, pe, nd.T_suc; phase = :vapor), T = nd.T_suc))
    W = compressor_power(compressor, port_values(S, u, :compressor))
    P_bh, P_ev = pumping_power(bs, BRINE_PROPS, ηPUMP), pumping_power(es, BRINE_PROPS, ηPUMP)
    E = (Qe = ph.evaporator.Qe, Qsup = ph.evaporator.Qsup, Qd = ph.condenser.Qd, Qc = ph.condenser.Qc,
         W = W, P_borehole = P_bh, P_evaporator = P_ev)
    Q1, Q2, P = E.Qd + E.Qc, E.Qe + E.Qsup, P_bh + P_ev
    return (name = r.name, r = r, Tb = nd.Tb, Vs = evaporator.Vs, states = states,
            nodes = (T_cold = nd.T_cold, T_warm = nd.T_warm, Tw_in = nd.Tw_in, Tw_out = nd.Tw_out),
            evap = evaporator_profiles(evaporator, ph.evaporator, pv_e; n = n),
            cond = condenser_profiles(condenser, ph.condenser, pv_c; n = n),
            energy = E, Q1 = Q1, Q2 = Q2, P = P, COP1 = Q1 / (W + P), COP2 = Q2 / (W + P))
end

"""
    cycle_sweep(r; Tbs = TBS, Vf = 1.0) -> Vector

Closed cycle of `r` at `Vf` × design brine flow for each borehole wall temperature in `Tbs`
(continuation, warm to cold); one [`cycle_snapshot`](@ref) per solved point.
"""
function cycle_sweep(r; Tbs = TBS, Vf = 1.0, design = refrigerant_design(r))
    build(Tb) = cycle_system(r; Vf = Vf, Tb = Tb, design = design)
    Tb1 = first(Tbs)
    res = continuation(build, Tbs, cycle_guess(build(Tb1), Tb1), TrustRegion(); abstol = 1e-10)
    for x in res
        x.ok || @warn "$(r.name): no closed cycle at Tb = $(x.λ - C2K) °C"
    end
    return [cycle_snapshot(build(x.λ), x.u) for x in res if x.ok]
end

tb_label(Tb) = @sprintf("Tb_%+03dC", round(Int, Tb - C2K))
tb_title(st) = @sprintf("%s, Tb = %.0f °C, brine %.2f l/s", st.name, st.Tb - C2K, 1e3st.Vs)
snapshot_at(sweep, Tb) = sweep[argmin([abs(s.Tb - Tb) for s in sweep])]

# ---- panels: each draws into a grid position, for the single figures and the overview ----

# saturation dome and vapour isotherms of the refrigerant
function dome!(ax, r; color = COLORS.dome, isotherms = true, label = nothing)
    # up to 0.75 p_crit: closer, the dew points of a blend (Clapeyron) jump to another root
    pmin, pmax = OptimalGshpFlow._p_dew(r, C2K - 35), 0.75 * OptimalGshpFlow.critical_pressure(r)
    ps = exp.(range(log(pmin), log(pmax), 80))
    lines!(ax, h_bubble.(Ref(r), ps) ./ 1e3, ps ./ 1e5; color = color, label = label)
    lines!(ax, h_dew.(Ref(r), ps) ./ 1e3, ps ./ 1e5; color = color)
    isotherms || return
    for T in C2K .+ (0, 20, 40, 60, 80)
        pT = OptimalGshpFlow._p_dew(r, T)
        pT < pmax || continue                         # near the critical point: not drawn
        pv = exp.(range(log(pmin), log(0.999pT), 40))
        lines!(ax, h_pT.(Ref(r), pv, T; phase = :vapor) ./ 1e3, pv ./ 1e5; color = (color, 0.35), linewidth = 0.8)
        text!(ax, h_pT(r, last(pv), T; phase = :vapor) / 1e3, last(pv) / 1e5; text = @sprintf("%.0f °C", T - C2K),
              fontsize = 9, color = (color, 0.8), align = (:left, :bottom), offset = (2, 1))
    end
end

# the cycle on the p-h diagram: compression (dashed), condenser, valve (dotted), evaporator
function cycle_path!(ax, st; color = COLORS.ref, label = st.name, annotate = true)
    s = st.states
    P(x) = Point2f(x.h / 1e3, x.p / 1e5)
    lines!(ax, [P(s.suction), P(s.discharge)]; color = color, linestyle = :dash)
    lines!(ax, [P(s.discharge), P(s.cond_dew), P(s.liquid)]; color = color, linewidth = 2.5, label = label)
    lines!(ax, [P(s.liquid), P(s.evap_in)]; color = color, linestyle = :dot)
    lines!(ax, [P(s.evap_in), P(s.evap_dew), P(s.suction)]; color = color, linewidth = 2.5)
    pts = (s.discharge, s.cond_dew, s.liquid, s.evap_in, s.evap_dew, s.suction)
    scatter!(ax, [P(x) for x in pts]; color = color, markersize = 7)
    annotate || return
    # number on a fixed side of each point (1–2 and 5–6 lie close together)
    sides = ((:left, :bottom, 5, 3), (:right, :bottom, -5, 3), (:right, :bottom, -5, 3),
             (:right, :top, -5, -3), (:right, :top, -5, -3), (:left, :top, 5, -3))
    for (i, (x, (ha, va, dx, dy))) in enumerate(zip(pts, sides))
        text!(ax, P(x); text = string(i), fontsize = 11, font = :bold, color = color, align = (ha, va), offset = (dx, dy))
    end
    names = ("discharge", "condenser, dew", "condenser outlet", "evaporator inlet", "evaporator, dew", "suction")
    text!(ax, 0.02, 0.98; space = :relative, align = (:left, :top), fontsize = 10, color = color,
          text = join((@sprintf("%d %-17s %6.1f °C %5.2f bar", i, n, x.T - C2K, x.p / 1e5) for (i, (n, x)) in
                       enumerate(zip(names, pts))), "\n"), font = "DejaVu Sans Mono")
end

function panel_ph!(gp, st)
    ax = Axis(gp; xlabel = "specific enthalpy [kJ/kg]", ylabel = "pressure [bar]", yscale = log10, yticks = PTICKS,
              title = "p-h diagram", subtitle = "compression dashed, expansion valve dotted")
    dome!(ax, st.r)
    cycle_path!(ax, st)
    return ax
end

function panel_profiles!(gp, st)
    g = GridLayout(gp)
    ae = Axis(g[1, 1]; xlabel = "area fraction (refrigerant inlet → outlet)", ylabel = "temperature [°C]",
              title = "Evaporator", subtitle = "brine ← counter-flow")
    ac = Axis(g[1, 2]; xlabel = "area fraction (refrigerant inlet → outlet)", title = "Condenser",
              subtitle = "water ← counter-flow")
    e, c = st.evap, st.cond
    lines!(ae, e.a, e.Ts .- C2K; color = COLORS.brine, linewidth = 2, label = "brine")
    lines!(ae, e.a, e.Tr .- C2K; color = COLORS.ref, linewidth = 2, label = "refrigerant")
    vlines!(ae, [e.xs]; color = :gray, linestyle = :dash)
    text!(ae, e.xs / 2, 0.02; text = "two-phase", space = :relative, align = (:center, :bottom), fontsize = 10)
    lines!(ac, c.a, c.Tw .- C2K; color = COLORS.water, linewidth = 2, label = "water")
    lines!(ac, c.a, c.Tr .- C2K; color = COLORS.ref, linewidth = 2, label = "refrigerant")
    vlines!(ac, [c.xd]; color = :gray, linestyle = :dash)
    # zone labels in data coordinates of the area fraction
    text!(ae, (1 + e.xs) / 2, minimum(e.Tr) - C2K; text = "superheat", align = (:center, :bottom), fontsize = 10)
    text!(ac, c.xd / 2, minimum(c.Tw) - C2K; text = "desuper-\nheating", align = (:center, :bottom), fontsize = 10)
    text!(ac, (1 + c.xd) / 2, minimum(c.Tw) - C2K; text = "condensation", align = (:center, :bottom), fontsize = 10)
    axislegend(ae; position = :lt, framevisible = false, labelsize = 10)
    axislegend(ac; position = :rt, framevisible = false, labelsize = 10)
    return g
end

const ENERGY_PARTS = (
    (key = :Qe, bar = 1, label = "ground → two-phase zone", color = :steelblue),
    (key = :Qsup, bar = 1, label = "ground → superheat zone", color = :lightblue),
    (key = :W, bar = 2, label = "compressor", color = :black),
    (key = :P_borehole, bar = 2, label = "pump, borehole", color = :purple),
    (key = :P_evaporator, bar = 2, label = "pump, evaporator", color = :orchid),
    (key = :Qc, bar = 3, label = "condensation → water", color = :darkorange),
    (key = :Qd, bar = 3, label = "desuperheating → water", color = :gold))

function panel_energy!(gp, st)
    E = st.energy
    g = GridLayout(gp)
    ax = Axis(g[1, 1]; ylabel = "power [kW]", xticks = (1:3, ["heat from\nthe ground", "electricity", "heat to\nthe water"]),
              title = "Energy", subtitle = @sprintf("COP1 = Q1/(W + P) = %.2f, COP2 = %.2f", st.COP1, st.COP2))
    xs = [p.bar for p in ENERGY_PARTS]
    ys = [getfield(E, p.key) / 1e3 for p in ENERGY_PARTS]
    grp = collect(1:length(ENERGY_PARTS))
    barplot!(ax, xs, ys; stack = grp, color = [p.color for p in ENERGY_PARTS])
    totals = (st.Q2, E.W + st.P, st.Q1)
    for (i, t) in enumerate(totals)
        text!(ax, i, t / 1e3; text = @sprintf("%.2f kW", t / 1e3), align = (:center, :bottom), offset = (0, 3), fontsize = 11)
    end
    ylims!(ax, 0, 1.15 * maximum(totals) / 1e3)
    Legend(g[1, 2], [PolyElement(color = p.color) for p in ENERGY_PARTS],
           [@sprintf("%s: %.0f W", p.label, getfield(E, p.key)) for p in ENERGY_PARTS]; framevisible = false, labelsize = 10)
    return ax
end

# ---- figures ----

function point_figures(st; dir = joinpath(FIGDIR, st.name, tb_label(st.Tb)))
    mkpath(dir)
    f = Figure(size = (760, 560)); panel_ph!(f[1, 1], st); Label(f[0, :], tb_title(st); fontsize = 14, tellwidth = false)
    save(joinpath(dir, "ph.png"), f; px_per_unit = 2)
    f = Figure(size = (1000, 420)); panel_profiles!(f[1, 1], st); Label(f[0, :], tb_title(st); fontsize = 14, tellwidth = false)
    save(joinpath(dir, "profiles.png"), f; px_per_unit = 2)
    f = Figure(size = (760, 440)); panel_energy!(f[1, 1], st); Label(f[0, :], tb_title(st); fontsize = 14, tellwidth = false)
    save(joinpath(dir, "energy.png"), f; px_per_unit = 2)
    f = Figure(size = (1400, 950))
    Label(f[0, :], tb_title(st); fontsize = 18, font = :bold, tellwidth = false)
    panel_ph!(f[1, 1], st)
    panel_energy!(f[1, 2], st)
    panel_profiles!(f[2, 1:2], st)
    colsize!(f.layout, 1, Relative(0.5))
    save(joinpath(dir, "overview.png"), f; px_per_unit = 2)
    return dir
end

# heat flows, electricity, COP and temperatures over the borehole wall temperature (one or more sweeps)
function vs_tb_figure(sweeps...; title = join((first(s).name for s in sweeps), " vs "))
    f = Figure(size = (1200, 850))
    Label(f[0, :], string(title, ": closed cycle at the design brine flow"); fontsize = 16, font = :bold, tellwidth = false)
    aq = Axis(f[1, 1]; ylabel = "heat [kW]", title = "Heat flows")
    ap = Axis(f[1, 2]; ylabel = "compressor [W]", title = "Electricity")
    ap2 = Axis(f[1, 2]; ylabel = "brine pumps [W]", yaxisposition = :right, ylabelcolor = :purple,
               yticklabelcolor = :purple, ytickcolor = :purple)
    hidespines!(ap2); hidexdecorations!(ap2); hideydecorations!(ap2; label = false, ticklabels = false, ticks = false)
    linkxaxes!(ap, ap2)
    ac = Axis(f[2, 1]; xlabel = "borehole wall temperature Tb [°C]", ylabel = "COP", title = "COP")
    at = Axis(f[2, 2]; xlabel = "borehole wall temperature Tb [°C]", ylabel = "temperature [°C]", title = "Temperatures")
    styles = (:solid, :dash, :dot)
    pel, plab = LineElement[], String[]                             # legend of the two electricity axes
    for (sw, ls) in zip(sweeps, styles)
        tb = [s.Tb - C2K for s in sw]
        nm = first(sw).name
        lines!(aq, tb, [s.Q1 / 1e3 for s in sw]; color = :darkorange, linestyle = ls, label = "Q1 to the water, $nm")
        lines!(aq, tb, [s.Q2 / 1e3 for s in sw]; color = :steelblue, linestyle = ls, label = "Q2 from the ground, $nm")
        lines!(ap, tb, [s.energy.W for s in sw]; color = :black, linestyle = ls)
        lines!(ap2, tb, [s.P for s in sw]; color = :purple, linestyle = ls)
        append!(pel, [LineElement(color = :black, linestyle = ls), LineElement(color = :purple, linestyle = ls)])
        append!(plab, ["compressor, $nm", "brine pumps (right axis), $nm"])
        lines!(ac, tb, [s.COP1 for s in sw]; color = :firebrick, linestyle = ls, label = "COP1, $nm")
        lines!(ac, tb, [s.COP2 for s in sw]; color = :gray40, linestyle = ls, label = "COP2, $nm")
        lines!(at, tb, [s.states.discharge.T - C2K for s in sw]; color = :firebrick, linestyle = ls, label = "discharge, $nm")
        lines!(at, tb, [s.states.cond_dew.T - C2K for s in sw]; color = :darkorange, linestyle = ls, label = "condensation (dew), $nm")
        lines!(at, tb, [s.nodes.Tw_out - C2K for s in sw]; color = :gold, linestyle = ls, label = "water out, $nm")
        lines!(at, tb, [s.nodes.T_warm - C2K for s in sw]; color = :steelblue, linestyle = ls, label = "brine to evaporator, $nm")
        lines!(at, tb, [s.nodes.T_cold - C2K for s in sw]; color = :lightblue, linestyle = ls, label = "brine from evaporator, $nm")
        lines!(at, tb, [s.states.evap_dew.T - C2K for s in sw]; color = :navy, linestyle = ls, label = "evaporation (dew), $nm")
    end
    for ax in (aq, ac)
        axislegend(ax; position = :lt, framevisible = false, labelsize = 9)
    end
    Legend(f[1, 3], pel, plab; framevisible = false, labelsize = 9)
    Legend(f[2, 3], at; framevisible = false, labelsize = 9)
    return f
end

"""
    cycle_figures(sweep; points = POINT_TBS)

Figures of one refrigerant from its [`cycle_sweep`](@ref): `vs_Tb.png` and, for each Tb of `points`,
the folder `Tb_±XXC/` with `ph.png`, `profiles.png`, `energy.png`, `overview.png`.
"""
function cycle_figures(sweep; points = POINT_TBS)
    name = first(sweep).name
    mkpath(joinpath(FIGDIR, name))
    save(joinpath(FIGDIR, name, "vs_Tb.png"), vs_tb_figure(sweep; title = name); px_per_unit = 2)
    dirs = [point_figures(snapshot_at(sweep, Tb)) for Tb in points]
    return [joinpath(FIGDIR, name, "vs_Tb.png"); dirs]
end

"""
    comparison_figures(a, b; Tb = 273.15)

Two refrigerants (their sweeps): p-h diagrams together at `Tb`, and the Tb curves together.
"""
function comparison_figures(a, b; Tb = C2K)
    dir = joinpath(FIGDIR, "comparison")
    mkpath(dir)
    sa, sb = snapshot_at(a, Tb), snapshot_at(b, Tb)
    f = Figure(size = (800, 560))
    ax = Axis(f[1, 1]; xlabel = "specific enthalpy [kJ/kg]", ylabel = "pressure [bar]", yscale = log10, yticks = PTICKS,
              title = @sprintf("p-h diagrams, Tb = %.0f °C, design brine flow", Tb - C2K))
    for (st, col) in ((sa, :firebrick), (sb, :teal))
        dome!(ax, st.r; color = (col, 0.6), isotherms = false)
        cycle_path!(ax, st; color = col, annotate = false,
                    label = @sprintf("%s: COP1 %.2f, Q1 %.2f kW", st.name, st.COP1, st.Q1 / 1e3))
    end
    axislegend(ax; position = :lt, framevisible = false)
    save(joinpath(dir, string("ph_", tb_label(Tb), ".png")), f; px_per_unit = 2)
    save(joinpath(dir, "vs_Tb.png"), vs_tb_figure(a, b); px_per_unit = 2)
    return dir
end
