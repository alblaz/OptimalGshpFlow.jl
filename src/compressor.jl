#=
Compressor and condenser (post-processing of an evaporator solution).

Subcritical cycle: the condensing pressure is the bubble pressure at `Tc`. A transcritical high
side (R744 gas cooler) is not covered.
=#

# condensing (bubble) pressure at T
_p_bubble(r::Refrigerant, T) = PropsSI("P", "T", T, "Q", 0, r.name)
_p_bubble(r::ClapeyronRefrigerant, T) = ispure(r) ? p_sat(r, T) : Clapeyron.bubble_pressure(r.eos, T, r.z)[1]

"""
    compression(refrigerant, pe, Tsuc, pc, ηis) -> (; h_in, h_out, ρ_in)

Compressor inlet enthalpy, outlet enthalpy (isentropic efficiency `ηis`) and inlet density, for
suction at pressure `pe` and temperature `Tsuc` and discharge at pressure `pc`.
"""
function compression(r::Refrigerant, pe, Tsuc, pc, ηis)
    # CoolProp by name, as in the original code
    h_in  = PropsSI("H", "T", Tsuc, "P", pe, r.name)
    s_in  = PropsSI("S", "T", Tsuc, "P", pe, r.name)
    h_is  = PropsSI("H", "S", s_in, "P", pc, r.name)
    ρ_in  = PropsSI("D", "T", Tsuc, "P", pe, r.name)
    return (h_in = h_in, h_out = h_in + (h_is - h_in)/ηis, ρ_in = ρ_in)
end

function compression(r::ClapeyronRefrigerant, pe, Tsuc, pc, ηis)
    h_in = h_pT(r, pe, Tsuc; phase = :vapor)
    s_in = s_pT(r, pe, Tsuc; phase = :vapor)
    h_is = h_ps(r, pc, s_in; phase = :vapor)
    ρ_in = ρ_pT(r, pe, Tsuc; phase = :vapor)
    return (h_in = h_in, h_out = h_in + (h_is - h_in)/ηis, ρ_in = ρ_in)
end

# saturated-liquid enthalpy at the condenser outlet (no subcooling)
_h_condenser_outlet(r::Refrigerant, Tc, pc) = PropsSI("H", "T", Tc, "Q", 0, r.name)
_h_condenser_outlet(r::ClapeyronRefrigerant, Tc, pc) = h_bubble(r, pc)

function _compression(Tc, ηis, es::EvaporatorSolution, r)
    pc = _p_bubble(r, Tc)                                  # condensation pressure
    pe = _p_dew(r, es.twophase_region.Tcl)                 # evaporation pressure
    return pc, compression(r, pe, es.superheat_region.Tch, pc, ηis)
end

"""
    compressor_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::AbstractRefrigerant; ξ = 0.07)

Power delivered by the compressor [W].

- `Tc`: condensing temperature [K];
- `ηis`: isentropic efficiency of the compressor;
- `evaporator_solution`: results of the evaporator (refrigerant volume flow, evaporation and
  suction temperatures);
- `ξ`: fraction of heat lost through the compressor envelope, default 7 % (lower for large
  hermetic compressors, larger for non-hermetic ones).

Spline `Refrigerant`: CoolProp by fluid name (original code); `ClapeyronRefrigerant`: the
property interface.
"""
function compressor_power(Tc, ηis, es::EvaporatorSolution, r::AbstractRefrigerant; ξ = 0.07)
    _, c = _compression(Tc, ηis, es, r)
    return c.ρ_in * es.twophase_region.Vc * (c.h_out - c.h_in) / (1 - ξ)
end

"""
    condenser_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::AbstractRefrigerant)

Power delivered at the condenser [W]: refrigerant mass flow × (compressor outlet enthalpy −
saturated-liquid enthalpy at `Tc`; no subcooling). Arguments as [`compressor_power`](@ref).
"""
function condenser_power(Tc, ηis, es::EvaporatorSolution, r::AbstractRefrigerant)
    pc, c = _compression(Tc, ηis, es, r)
    return c.ρ_in * es.twophase_region.Vc * (c.h_out - _h_condenser_outlet(r, Tc, pc))
end
