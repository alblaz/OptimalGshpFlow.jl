#=
(T, v) state layer for `ClapeyronRefrigerant`.

With temperature T [K] and molar volume v [m³/mol] as the state, every property is an explicit
evaluation of the Helmholtz energy (no inner iteration, exact automatic differentiation). A model
that carries (T, v) among its unknowns closes the state with the pressure residual
`pressure_residual(r, T, v, p) = p_Tv(r, T, v)/p - 1`; the outer Newton iteration then acts as the
warm start of the volume solve.

The volume is molar (as in Clapeyron); all outputs are per kg, as in the rest of the package.
=#

"""
    volume_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown)

Molar volume [m³/mol] at `(p, T)` (inner volume solve). Use it to build initial guesses for the
(T, v) unknowns.
"""
volume_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = Clapeyron.volume(r.eos, p, T, r.z; phase = phase)

"Pressure [Pa] at temperature `T` [K] and molar volume `v` [m³/mol]."
p_Tv(r::ClapeyronRefrigerant, T, v) = Clapeyron.pressure(r.eos, v, T, r.z)
"Specific enthalpy [J/kg] at `(T, v)`."
h_Tv(r::ClapeyronRefrigerant, T, v) = Clapeyron.VT_enthalpy(r.eos, v, T, r.z) / r.M
"Specific entropy [J/kg/K] at `(T, v)`."
s_Tv(r::ClapeyronRefrigerant, T, v) = Clapeyron.VT_entropy(r.eos, v, T, r.z) / r.M
"Isobaric specific heat [J/kg/K] at `(T, v)`."
cp_Tv(r::ClapeyronRefrigerant, T, v) = Clapeyron.VT_isobaric_heat_capacity(r.eos, v, T, r.z) / r.M
"Density [kg/m³] at molar volume `v` [m³/mol]."
ρ_Tv(r::ClapeyronRefrigerant, T, v) = r.M / v

"Viscosity [Pa s] at `(T, v)`."
μ_Tv(r::ClapeyronRefrigerant, T, v) = _transport_Tv(r, :viscosity, T, v)
"Thermal conductivity [W/m/K] at `(T, v)`."
λ_Tv(r::ClapeyronRefrigerant, T, v) = _transport_Tv(r, :conductivity, T, v)

# EntropyScaling evaluates the residual entropy of its own equation of state at (T, v)
function _transport_Tv(r::ClapeyronRefrigerant{<:Any, EntropyScalingTransport}, prop, T, v)
    if prop === :viscosity
        return EntropyScaling.VT_viscosity(r.transport_model, v, T, r.z)
    else
        return EntropyScaling.VT_thermal_conductivity(r.transport_model, v, T, r.z)
    end
end

# CoolProp with (mass density, T) inputs: no flash, derivatives by central differences
_transport_Tv(r::ClapeyronRefrigerant{<:Any, CoolPropTransport}, prop, T, v) =
    _coolprop(r.transport_model, COOLPROP_DT, _coolprop_key(prop), r.M / v, T, 1e-6, 1e-6)

"""
    pressure_residual(r, T, v, p)

Scaled residual `p_Tv(r, T, v)/p - 1` that closes a (T, v) state at pressure `p`.
"""
pressure_residual(r::ClapeyronRefrigerant, T, v, p) = p_Tv(r, T, v) / p - 1

"""
    volume_from_R(r, u) = r.vmin * (1 + exp(u))
    volume_to_R(r, v)   = log(v / r.vmin - 1)

Map between the molar volume, bounded below by the co-volume `r.vmin`, and an unbounded unknown
on a log scale (like the pressure). Liquid states map to `u ≈ 0`, vapour at low reduced
pressure to `u ≈ 4–6`.
"""
volume_from_R(r::ClapeyronRefrigerant, u) = r.vmin * (1 + exp(u))
volume_to_R(r::ClapeyronRefrigerant, v) = log(v / r.vmin - 1)
