#= Closed heat-pump cycle (step 9d): R290 and R454C with the condenser on the water side.

Evaporator + compressor + condenser + borehole solved together (`heatpump_cycle`): the condensing
pressure follows from the water side (R290-sized condenser of `examples/condenser.jl`, floor heating
water in at 30 °C, fixed water flow), the evaporator inlet enthalpy from the condenser outlet.
Each point is compared with the open cycle of the optimal-flow study (evaporator + borehole,
`hei` and condensing temperature from a convention, 5h) at the same ground temperature and brine
flow: which convention (`:bubble`, `:midpoint`, `:dew`) the closed cycle supports for R454C, and
the COP difference between the two refrigerants. Machine size as in 5h (same nominal two-phase
duty, `refrigerant_design(r)`, midpoint). No suction pressure drop on either side (compared like
for like).

Usage:
    include("examples/cycle.jl")
    r290 = ClapeyronRefrigerant("R290")
    r454c = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
    rows = [compare_cycle(r290); compare_cycle(r454c)]; report_cycle(rows)
(one refrigerant per call: a vector of both has an abstract element type, and the loop body would
be inferred on the abstract `ClapeyronRefrigerant`, Clapeyron included: no result in 30 min)
=#
include(joinpath(@__DIR__, "condenser.jl"))

const TBS_CYCLE = 273.15 .+ (8.0, 0.0, -5.0)

# closed cycle of the refrigerant `r` (same machine as the open study: `refrigerant_design(r)`)
function cycle_system(r; Vf = 1.0, Tb = TB, Tw_in = TW_IN, design = refrigerant_design(r), sf = BRINE_PROPS)
    Vs = Vf * DESIGN.Vs
    evap = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = PHX, hei = nothing,
                      hei_nominal = design.hei, Vs = Vs, Vref = design.Vref, ΔTsup = ΔTSUP)
    comp = Compressor(refrigerant = r, ηis = ηIS)
    cond = Condenser(refrigerant = r, secondaryfluid = WATER, plateheatexchanger = CPHX, Vs = CONDENSER_DESIGN.Vs,
                     Qscale = CONDENSER_DESIGN.Q1)
    bhx = BoreholeHX(borehole = BOREHOLE, secondaryfluid = sf, Vs = Vs, λground = λGROUND)
    return heatpump_cycle(evaporator = evap, compressor = comp, condenser = cond, borehole = bhx, Tb = Tb, Tw_in = Tw_in)
end

cycle_guess(S, Tb) = heatpump_guess(S; Te = Tb - 8, Tc = TCOND, T_cold = Tb - 6, T_warm = Tb - 2, Tw_out = TW_OUT)

# NonlinearSolve on the closed cycle, results through barriers (`solve` is inferred as `Any`: its
# `u` passed on to the post-processing would infer the model, Clapeyron included, on it)
function cycle_solve(S, u0; alg = TrustRegion(), abstol = 1e-10)
    sol = solve(NonlinearProblem(S, u0), alg; abstol = abstol)
    u = convert(Vector{Float64}, sol.u)::Vector{Float64}
    return (ok = (SciMLBase.successful_retcode(sol)::Bool && all(isfinite, u)), u = u, iterations = Int(sol.stats.nsteps)::Int)
end

# performance of a closed-cycle solution
function cycle_point(S, u)
    (; evaporator, compressor, condenser, borehole) = S.components
    r = evaporator.refrigerant
    ph = physical(S, u); n = ph.nodes
    es = evaporator_solution(evaporator, ph.evaporator, port_values(S, u, :evaporator))
    bs = borehole_solution(borehole, port_values(S, u, :borehole))
    W = compressor_power(compressor, port_values(S, u, :compressor))
    Q1, Q2 = ph.condenser.Qd + ph.condenser.Qc, ph.evaporator.Qe + ph.evaporator.Qsup
    P = pump_power(es, bs, BRINE_PROPS)
    return (Te = T_dew(r, n.p_suc), T_bubble = T_bubble(r, n.pc), T_dew = T_dew(r, n.pc), Tdis = vapour_T_ph(r, n.pc, n.h_dis),
            Tw_out = n.Tw_out, Ts_out = n.T_cold, m = n.m, Q1 = Q1, Q2 = Q2, W = W, Ppump = P, COP1 = Q1 / (W + P),
            issues = String[i.message for i in validity_issues(S, u)])
end

# open cycle (the study) at the same point, condenser convention `conv`
function open_point(r, Tb, conv; Vf = 1.0, u0 = nothing)
    design = refrigerant_design(r; condenser = conv)
    S = heatpump(r, Vf; Tb = Tb, design = design)
    sol = solve(NonlinearProblem(S, u0 === nothing ? first_guess(S, Tb) : u0), TrustRegion(); abstol = 1e-9)
    SciMLBase.successful_retcode(sol) || return nothing
    u = convert(Vector{Float64}, sol.u)::Vector{Float64}   # barrier: `sol.u` is not inferred (blend: Clapeyron on `Any`)
    v = evaluate(S, u, r, BRINE_PROPS; Tc = design.Tc)
    p = performance(v.es, v.bs, r, BRINE_PROPS; Δp_suction = 0.0, Tc = design.Tc)
    return (u = u, Te = v.es.twophase_region.Tcl, Tc = design.Tc, p...)
end

"""
    compare_cycle(r; Tbs = TBS_CYCLE, Vf = 1.0) -> Vector

For the refrigerant `r` and each borehole wall temperature: the closed cycle (`cycle_solve`, warm start
from the previous `Tb`, warmest first) and the open cycle with each condenser convention
(`:bubble`, `:midpoint`, `:dew`; one for a pure fluid).
"""
function compare_cycle(r; Tbs = TBS_CYCLE, Vf = 1.0, log = devnull)
    S = cycle_system(r; Vf = Vf, Tb = first(Tbs))
    u = cycle_guess(S, first(Tbs))
    convs = OptimalGshpFlow.ispure(r) ? (:midpoint,) : (:bubble, :midpoint, :dew)
    rows = []
    for Tb in Tbs
        St = set_fixed(S; Tb = Tb)
        _log(log, r.name, " Tb ", Tb - 273.15, ": closed cycle")
        x = cycle_solve(St, u)
        x.ok || (@warn "$(r.name): closed cycle failed at Tb = $(Tb - 273.15) °C"; break)
        u = x.u
        _log(log, "  newton ", x.iterations, " its; post-processing")
        c = cycle_point(St, u)
        opens = map(convs) do conv
            _log(log, "  open cycle ", conv)
            conv => open_point(r, Tb, conv; Vf = Vf)
        end
        push!(rows, (name = r.name, Tb = Tb, closed = c, open = opens))
    end
    return rows
end

_log(io, args...) = (println(io, args...); flush(io); nothing)

function report_cycle(rows)
    @printf("%-6s %5s | %7s %7s %7s %7s %7s %6s %6s %6s %6s | %s\n", "", "Tb", "Te", "T_bub", "T_dew", "Tdis", "Tw_out",
            "Q1", "W", "Ppump", "COP1", "open cycle: convention Tc → COP1 (Δ vs closed)")
    for x in rows
        c = x.closed
        @printf("%-6s %5.1f | %7.2f %7.2f %7.2f %7.1f %7.2f %6.0f %6.0f %6.1f %6.3f |", x.name, x.Tb - 273.15, c.Te - 273.15,
                c.T_bubble - 273.15, c.T_dew - 273.15, c.Tdis - 273.15, c.Tw_out - 273.15, c.Q1, c.W, c.Ppump, c.COP1)
        for (conv, o) in x.open
            o === nothing ? print(" $conv: failed") :
                @printf(" %s %.2f → %.3f (%+.1f %%)", conv, o.Tc - 273.15, o.COP1, 100(o.COP1 / c.COP1 - 1))
        end
        println(isempty(c.issues) ? "" : "  issues: " * join(c.issues, "; "))
    end
    # the blend against the pure fluid, closed cycle
    for Tb in unique(x.Tb for x in rows)
        rs = filter(x -> x.Tb == Tb, rows)
        length(rs) == 2 || continue
        @printf("Tb %5.1f °C: %s / %s COP1 %+.1f %%, Q1 %+.1f %%\n", Tb - 273.15, rs[2].name, rs[1].name,
                100(rs[2].closed.COP1 / rs[1].closed.COP1 - 1), 100(rs[2].closed.Q1 / rs[1].closed.Q1 - 1))
    end
end
