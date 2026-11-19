#=
Compressor and condenser (post-processing of an evaporator solution).

Subcritical cycle: the condensing pressure is the bubble pressure at `Tc`. A transcritical high
side (R744 gas cooler) is not covered.
=#

# condensing (bubble) pressure at T
_p_bubble(r::Refrigerant, T) = _coolprop_value(r.coolprop, COOLPROP_QT, COOLPROP_P, 0.0, T)
_p_bubble(r::ClapeyronRefrigerant, T) = ispure(r) ? p_sat(r, T) : Clapeyron.bubble_pressure(r.eos, T, r.z)[1]

"""
    compression(refrigerant, pe, Tsuc, Δp, pc, ηis) -> (; h_in, h_out, ρ_in)

Compression of the vapour leaving the evaporator at pressure `pe` and temperature `Tsuc`, after an
adiabatic suction pressure drop `Δp` (isenthalpic: the compressor inlet is at `pe − Δp` with the
evaporator outlet enthalpy), to the discharge pressure `pc` with isentropic efficiency `ηis`.
Returns the inlet enthalpy `h_in` (= evaporator outlet), the outlet enthalpy `h_out`, and the
density `ρ_in` at the evaporator outlet: the refrigerant volume flow `Vref` of the evaporator is
referred to it, so `ρ_in Vref` is the evaporator's mass flow (B11).
"""
function compression(r::Refrigerant, pe, Tsuc, Δp, pc, ηis)
    # CoolProp (HEOS, as the original PropsSI calls) through the cached state: one flash for the
    # evaporator outlet, one for the suction state if Δp ≠ 0, one for the isentropic outlet
    h_in, s_e, ρ_in = _coolprop_values(r.coolprop, COOLPROP_PT, pe, Tsuc, (COOLPROP_H, COOLPROP_S, COOLPROP_D))
    s_in = iszero(Δp) ? s_e : _coolprop_value(r.coolprop, COOLPROP_HP, COOLPROP_S, h_in, pe - Δp)
    h_is = _coolprop_value(r.coolprop, COOLPROP_PS, COOLPROP_H, pc, s_in)
    return (h_in = h_in, h_out = h_in + (h_is - h_in)/ηis, ρ_in = ρ_in)
end

function compression(r::ClapeyronRefrigerant, pe, Tsuc, Δp, pc, ηis)
    h_in = h_pT(r, pe, Tsuc; phase = :vapor)
    ρ_in = ρ_pT(r, pe, Tsuc; phase = :vapor)
    s_in = if iszero(Δp)
        s_pT(r, pe, Tsuc; phase = :vapor)
    else
        # suction temperature at constant enthalpy (Joule-Thomson, ≈ 0.03 K per 10 kPa for R290):
        # Newton in T at p = pe − Δp from Tsuc; two steps reach machine precision
        ps, T = pe - Δp, Tsuc
        for _ in 1:2
            T -= (h_pT(r, ps, T; phase = :vapor) - h_in) / cp_pT(r, ps, T; phase = :vapor)
        end
        s_pT(r, ps, T; phase = :vapor)
    end
    h_is = h_ps(r, pc, s_in; phase = :vapor)
    return (h_in = h_in, h_out = h_in + (h_is - h_in)/ηis, ρ_in = ρ_in)
end

# saturated-liquid enthalpy at the condenser outlet (no subcooling)
_h_condenser_outlet(r::Refrigerant, Tc, pc) = _coolprop_value(r.coolprop, COOLPROP_QT, COOLPROP_H, 0.0, Tc)
_h_condenser_outlet(r::ClapeyronRefrigerant, Tc, pc) = h_bubble(r, pc)

function _compression(Tc, ηis, es::EvaporatorSolution, r, Δp_suction)
    pc = _p_bubble(r, Tc)                                  # condensation pressure
    pe = _p_dew(r, es.twophase_region.Tcl)                  # evaporator outlet pressure
    return pc, compression(r, pe, es.superheat_region.Tch, Δp_suction, pc, ηis)
end

"""
    compressor_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::AbstractRefrigerant; ξ = 0.07)

Power delivered by the compressor [W].

- `Tc`: condensing temperature [K];
- `ηis`: isentropic efficiency of the compressor;
- `evaporator_solution`: results of the evaporator (refrigerant volume flow, evaporation and
  suction temperatures);
- `ξ`: fraction of heat lost through the compressor envelope, default 7 % (lower for large
  hermetic compressors, larger for non-hermetic ones);
- `Δp_suction`: pressure drop from the end of the two-phase zone to the compressor inlet [Pa]
  (e.g. the superheat zone, [`evaporator_pressure_drop`](@ref)); default 0. Adiabatic: the
  compressor inlet has the evaporator outlet enthalpy at the lower pressure; the mass flow is the
  evaporator's (`Vref` referred to the evaporator outlet), so here the pressure drop raises the
  work per kg but not the capacity (a lower suction density would need the suction pressure in
  the evaporator residual: not modelled).

Spline `Refrigerant`: CoolProp (HEOS, cached state); `ClapeyronRefrigerant`: the property
interface. For both powers use [`compressor_condenser_power`](@ref) (one compression).
"""
function compressor_power(Tc, ηis, es::EvaporatorSolution, r::AbstractRefrigerant; ξ = 0.07, Δp_suction = 0.0)
    _, c = _compression(Tc, ηis, es, r, Δp_suction)
    return c.ρ_in * es.twophase_region.Vc * (c.h_out - c.h_in) / (1 - ξ)
end

"""
    condenser_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::AbstractRefrigerant)

Power delivered at the condenser [W]: refrigerant mass flow × (compressor outlet enthalpy −
saturated-liquid enthalpy at `Tc`; no subcooling). Arguments as [`compressor_power`](@ref).
"""
function condenser_power(Tc, ηis, es::EvaporatorSolution, r::AbstractRefrigerant; Δp_suction = 0.0)
    pc, c = _compression(Tc, ηis, es, r, Δp_suction)
    return c.ρ_in * es.twophase_region.Vc * (c.h_out - _h_condenser_outlet(r, Tc, pc))
end

"""
    compressor_condenser_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant; ξ = 0.07) -> (; W, Q1)

[`compressor_power`](@ref) `W` and [`condenser_power`](@ref) `Q1` from one compression
calculation (half the cost of calling both).
"""
function compressor_condenser_power(Tc, ηis, es::EvaporatorSolution, r::AbstractRefrigerant; ξ = 0.07, Δp_suction = 0.0)
    pc, c = _compression(Tc, ηis, es, r, Δp_suction)
    m = c.ρ_in * es.twophase_region.Vc
    return (W = m * (c.h_out - c.h_in) / (1 - ξ), Q1 = m * (c.h_out - _h_condenser_outlet(r, Tc, pc)))
end
