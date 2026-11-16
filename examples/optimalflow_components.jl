#= Optimal brine flow of a ground-source heat pump (evaporator + borehole), with the component API.

Port of `heatpumpborehole_optimalflowrate.jl` (same settings): sweep the brine flow, solve the
evaporator coupled to the borehole by continuation, post-process compressor, condenser and pumps,
and find the flow that maximises the COP. Each point is checked with `validity_issues` (e.g. brine
below its freezing point at low flow).

Needs NonlinearSolve in the environment. Usage:
    include("examples/optimalflow_components.jl")
    rows = study(ClapeyronRefrigerant("R290")); report(rows)
    old  = original(); compare(rows, old)          # original API (spline Refrigerant)
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

# ---- model ----
function heatpump(r, Vf; sf = BRINE_PROPS, Tb = TB)
    Vs = Vf * DESIGN.Vs
    evap = Evaporator(refrigerant = r, secondaryfluid = sf, plateheatexchanger = PHX, hei = DESIGN.hei,
                      Vs = Vs, Vref = DESIGN.Vref, ΔTsup = ΔTSUP)
    bhx = BoreholeHX(borehole = BOREHOLE, secondaryfluid = sf, Vs = Vs, λground = λGROUND)
    return ComponentSystem((evaporator = evap, borehole = bhx);
        nodes = (T_cold = (:evaporator => :Ts_out, :borehole => :Tin),
                 T_warm = (:evaporator => :Ts_in,  :borehole => :Tout),
                 T_ground = (:borehole => :Tb,)), fixed = (T_ground = Tb,))
end

# pumping power: borehole (U-pipe, 2H) and evaporator plates (as in the original script)
function pump_power(es::EvaporatorSolution, bs::BoreholeSolution, sf)
    Vh = es.superheat_region.Vh
    P_bh = pumping_power(ηPUMP, 2BOREHOLE.H, Vh, 2BOREHOLE.rpi, sf, (bs.Tin + bs.Tout)/2)
    Tmean = (es.superheat_region.Thh + es.twophase_region.Thl)/2
    P_ev = pumping_power(pressure_drop_evap(Vh, PHX, sf, Tmean), Vh, ηPUMP)
    return P_bh + P_ev
end

function performance(es::EvaporatorSolution, bs::BoreholeSolution, r, sf)
    W  = compressor_power(TCOND, ηIS, es, r)
    Q1 = condenser_power(TCOND, ηIS, es, r)
    Q2 = es.twophase_region.Q + es.superheat_region.Q
    P  = pump_power(es, bs, sf)
    return (W = W, Q1 = Q1, Q2 = Q2, Ppump = P, COP1 = Q1/(W + P), COP2 = Q2/(W + P))
end

"""
    study(r; flows = FLOWS, sf = BRINE_PROPS, alg = TrustRegion()) -> Vector

Continuation in the brine flow, from `first(flows)` down; one row per flow with the physical
solution, the performance and the validity issues.
"""
function study(r; flows = FLOWS, sf = BRINE_PROPS, alg = TrustRegion())
    S1 = heatpump(r, first(flows); sf = sf)
    u0 = initial_guess(S1; evaporator = (Te = TB - 8, Qe = 0.9DESIGN.Qnom, Qsup = 0.1DESIGN.Qnom, Ts_m = TB - 4),
                       nodes = (T_cold = TB - 6, T_warm = TB - 2))
    res = continuation(Vf -> heatpump(r, Vf; sf = sf), flows, u0, alg; abstol = 1e-9)
    rows = []
    for x in res
        x.ok || (@warn "no solution at $(x.λ) × design flow"; continue)
        S = heatpump(r, x.λ; sf = sf)
        ph = physical(S, x.u)
        es = evaporator_solution(S.components.evaporator, ph.evaporator, port_values(S, x.u, :evaporator))
        bs = borehole_solution(S.components.borehole, port_values(S, x.u, :borehole))
        push!(rows, (Vf = x.λ, Vs = x.λ * DESIGN.Vs, Te = es.twophase_region.Tcl, Ts_out = ph.nodes.T_cold,
                     Ts_in = ph.nodes.T_warm, iterations = x.iterations, performance(es, bs, r, sf)...,
                     issues = unique(i.message for i in validity_issues(S, x.u))))
    end
    return rows
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

function compare(new, old)
    worst = (dTe = 0.0, dCOP1 = 0.0, dCOP2 = 0.0)
    for a in new, b in old
        a.Vf ≈ b.Vf || continue
        worst = (dTe = max(worst.dTe, abs(a.Te - b.Te)), dCOP1 = max(worst.dCOP1, abs(a.COP1 - b.COP1)/b.COP1),
                 dCOP2 = max(worst.dCOP2, abs(a.COP2 - b.COP2)/b.COP2))
    end
    return worst
end
