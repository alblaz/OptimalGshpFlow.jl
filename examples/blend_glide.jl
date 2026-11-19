#= Evaporator with a zeotropic blend: temperature glide in the two-phase zone (ROADMAP step 5e,
docs/dev/glide.md).

R454C (R32/R1234yf 21.5/78.5 % by mass, ≈ 8 K glide between bubble and dew point). The refrigerant
enters the two-phase zone at T_r,in = T(pe, hei) (p-h flash), below the dew temperature Te at the
end of the zone, and the model takes its temperature linear along the zone (Qiao 2025: pseudo
capacity C_r = ṁ (h_dew − hei)/(Te − T_r,in), φ = C_s/C_r in the effectiveness).

- `glide_sweep(r)`: standalone evaporator (brine inlet fixed) over the brine flow: glide, φ (crosses
  1: balanced counter-flow), duty, temperatures.
- `reference_error(S, u)`: the two-phase zone integrated along the area with the exact T(h) of the
  blend (flash at each step, RK4) from the model's brine outlet and refrigerant inlet: the
  mismatch at the end of the zone (brine temperature vs Ts_m, enthalpy vs h_dew) measures the error
  of the linear temperature profile.

Needs NonlinearSolve in the environment. Usage:
    include("examples/blend_glide.jl")
    r = ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785], name = "R454C")
    rows = glide_sweep(r); report_glide(rows)
    c = pressure_drop_effect(r); report_pressure_drop_effect(c)   # with/without refrigerant Δp
=#
using OptimalGshpFlow, NonlinearSolve, Printf
const Clapeyron = OptimalGshpFlow.Clapeyron

const BRINE = "INCOMP::MEA-28%"
const SF = SecondaryFluid(range(-15. + 273.15, 30 + 273.15, length = 40), BRINE)
const DESIGN = EvaporatorDesign(6000., 273., 273 + 40, "R290", BRINE; U = 2000., ΔTs = 3., ΔT_l = 1.5)
const PHX = PlateHeatExchanger(DESIGN.A, 0.5, 0.15; b = 0.002, ϕ = 1.189, β = 30.)
const TCOND, TS_IN, ΔTSUP = 273.15 + 40, 273.15 + 5, 3.5

# inlet enthalpy: saturated liquid at the condensing (bubble) pressure
inlet_enthalpy(r) = h_bubble(r, Clapeyron.bubble_pressure(r.eos, TCOND, r.z)[1])

# standalone evaporator, brine inlet fixed at TS_IN, brine flow Vf × design
function evaporator(r, Vf; hei = inlet_enthalpy(r), refrigerant_pressure_drop = default_refrigerant_pressure_drop(r))
    e = Evaporator(; refrigerant = r, secondaryfluid = SF, plateheatexchanger = PHX, hei = hei,
                   Vs = Vf * DESIGN.Vs, Vref = DESIGN.Vref, ΔTsup = ΔTSUP,
                   refrigerant_pressure_drop = refrigerant_pressure_drop)
    return ComponentSystem((evaporator = e,); nodes = (T_out = (:evaporator => :Ts_out,), T_in = (:evaporator => :Ts_in,)),
                           fixed = (T_in = TS_IN,))
end

"""
    glide_sweep(r; flows = 2.0:-0.1:0.3, refrigerant_pressure_drop = default_refrigerant_pressure_drop(r)) -> Vector

Continuation in the brine flow (× design) from `first(flows)`; one row per flow: zone quantities
(`Te − T_r,in`, `φ`, `γ`, `NTU`, `Δp`), duties, brine temperatures and validity issues.
`refrigerant_pressure_drop = nothing`: glide only.
"""
function glide_sweep(r; flows = 2.0:-0.1:0.3, refrigerant_pressure_drop = default_refrigerant_pressure_drop(r))
    hei = inlet_enthalpy(r)
    hp(Vf) = evaporator(r, Vf; hei = hei, refrigerant_pressure_drop = refrigerant_pressure_drop)
    S0 = hp(first(flows))
    u0 = initial_guess(S0; evaporator = (Te = TS_IN - 10, Qe = 3000., Qsup = 300., Ts_m = TS_IN - 1),
                       nodes = (T_out = TS_IN - 3,))
    res = continuation(hp, flows, u0, TrustRegion(); abstol = 1e-9)
    for x in res
        x.ok || @warn "no solution at $(x.λ) × design flow"
    end
    # concretely typed rows (a comprehension, not push! into []): the report calls model functions
    # on the fields; on a Vector{Any} they would be inferred on abstract types (≈ 5 min, PR #32)
    return [sweep_row(hp(x.λ), x) for x in res if x.ok]
end

function sweep_row(S, x)
    ph = physical(S, x.u); pv = port_values(S, x.u, :evaporator)
    z = evaporator_zones(S.components.evaporator, ph.evaporator, pv).twophase
    return (Vf = x.λ, S = S, u = x.u, Te = z.Te, Tr_in = z.Tr_in, φ = z.φ, γ = z.γ, Δp = z.Δp, NTU = z.NTU, Qe = ph.evaporator.Qe,
            Qsup = ph.evaporator.Qsup, xs = ph.evaporator.xs, Ts_out = pv.Ts_out, Ts_m = ph.evaporator.Ts_m,
            issues = String[i.message for i in validity_issues(S, x.u)] |> unique)
end

"""
    reference_error(S, u; n = 40) -> (; ΔTs, Δh, ΔQ)

Two-phase zone of the solution `u` of the system `S` integrated along the area (RK4, `n` steps)
with the exact refrigerant temperature T(p, h) (flash; pressure linear along the zone from
`pe + Δp` to `pe`, as in the model) and the model's U value and capacities,
from the refrigerant inlet (enthalpy `hei`, brine leaving at `Ts_out`). At the end of the zone the
model has the brine at `Ts_m` and the refrigerant at `h_dew`: `ΔTs` and `Δh` are the mismatches,
`ΔQ = ṁ Δh` the duty error of the linear profile [W] (positive: the exact profile absorbs more).
"""
function reference_error(S, u; n = 40)
    e = S.components.evaporator
    r = e.refrigerant
    ph = physical(S, u); pv = port_values(S, u, :evaporator)
    x = ph.evaporator
    z = evaporator_zones(e, x, pv).twophase
    pe = x.pe
    m = ρ_pT(r, pe, z.Te + e.ΔTsup; phase = :vapor) * e.Vref
    UA = z.NTU * z.Cs                                   # U A_zone
    # state (Ts, h) along s = A/A_zone from the refrigerant inlet; the brine flows towards s = 0
    f(s, y) = (q = UA * (y[1] - T_ph(r, pe + z.Δp * (1 - s), y[2])); (q / z.Cs, q / m))
    y = (pv.Ts_out, e.hei); ds = 1/n; s = 0.0
    for _ in 1:n
        k1 = f(s, y)
        k2 = f(s + ds/2, y .+ ds/2 .* k1)
        k3 = f(s + ds/2, y .+ ds/2 .* k2)
        k4 = f(s + ds, y .+ ds .* k3)
        y = y .+ ds/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); s += ds
    end
    Δh = y[2] - h_dew(r, pe)
    return (ΔTs = y[1] - x.Ts_m, Δh = Δh, ΔQ = m * Δh)
end

function report_glide(rows; reference = true)
    @printf("%5s %8s %8s %7s %6s %6s %7s %7s %8s %8s %7s", "Vf", "Te[°C]", "Tr,in", "glide", "φ", "NTU", "Qe[W]",
            "Qsup[W]", "Ts,out", "Ts,m", "xs")
    reference && @printf(" %8s %8s", "ΔTs[K]", "ΔQe[%]")
    println("  issues")
    for x in rows
        @printf("%5.2f %8.2f %8.2f %7.2f %6.3f %6.3f %7.0f %7.0f %8.2f %8.2f %7.3f", x.Vf, x.Te - 273.15, x.Tr_in - 273.15,
                x.Te - x.Tr_in, x.φ, x.NTU, x.Qe, x.Qsup, x.Ts_out - 273.15, x.Ts_m - 273.15, x.xs)
        if reference
            d = reference_error(x.S, x.u)
            @printf(" %8.3f %8.2f", d.ΔTs, 100d.ΔQ / x.Qe)
        end
        println("  ", join(x.issues, "; "))
    end
    println("glide = Te − T_r,in (two-phase zone); φ = C_brine/C_r; ΔTs, ΔQe: exact T(h) profile vs the linear one.")
end

"""
    pressure_drop_effect(r; flows = 2.0:-0.1:0.3) -> (without, with)

`glide_sweep` without and with the refrigerant pressure drop (default model, Amalfi 2016);
see [`report_pressure_drop_effect`](@ref).
"""
pressure_drop_effect(r; flows = 2.0:-0.1:0.3) =
    (without = glide_sweep(r; flows = flows, refrigerant_pressure_drop = nothing), with = glide_sweep(r; flows = flows))

function report_pressure_drop_effect(c)
    @printf("%5s | %8s %8s %7s | %8s %8s %7s %7s %6s %6s | %7s %7s\n", "Vf", "Te", "Tr,in", "Qe", "Te", "Tr,in", "Qe",
            "Δp kPa", "γ", "φ", "ΔTe[K]", "ΔQe[%]")
    for (a, b) in zip(c.without, c.with)
        a.Vf ≈ b.Vf || continue
        @printf("%5.2f | %8.2f %8.2f %7.0f | %8.2f %8.2f %7.0f %7.2f %6.3f %6.3f | %7.3f %7.2f\n", a.Vf, a.Te - 273.15,
                a.Tr_in - 273.15, a.Qe, b.Te - 273.15, b.Tr_in - 273.15, b.Qe, b.Δp/1e3, b.γ, b.φ, b.Te - a.Te, 100(b.Qe/a.Qe - 1))
    end
    println("left: glide only; right: glide and refrigerant pressure drop (Amalfi 2016). Temperatures in °C.")
end
