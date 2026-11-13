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

########################
# saturation at T      #
########################

"""
    saturation_T(r::ClapeyronRefrigerant, T) -> (p, v_l, v_v)

Saturation pressure [Pa] and molar volumes [m³/mol] of saturated liquid and vapour of a pure
fluid at temperature `T` [K].

When the CoolProp fluid file provides superancillary equations (Bell et al., Chebyshev expansions
in T accurate to ≈ 1e-12, `r.superancillary == true`; true for all common refrigerants) this is
an explicit evaluation (≈ 0.2 µs, exact AD) with no inner iteration. Otherwise it falls back to
Clapeyron's saturation solver.

Valid for `r.Ttriple < T < r.Tcrit`; the superancillaries return `NaN` outside. Combine with the
(T, v) functions, e.g. the dew-point enthalpy is `h_Tv(r, T, v_v)`. Mixtures (dew and bubble
points with glide) are not covered (see the roadmap, glide branch).
"""
function saturation_T(r::ClapeyronRefrigerant, T)
    ispure(r) || throw(ArgumentError("saturation_T is implemented for pure fluids only"))
    return r.superancillary ? _saturation_superanc(r.eos.pures[1], T) : _saturation_solve(r.eos.pures[1], T)
end

"Saturation pressure [Pa] of a pure fluid at temperature `T` [K]; see [`saturation_T`](@ref)."
p_sat(r::ClapeyronRefrigerant, T) = saturation_T(r, T)[1]

# Clapeyron loads the superancillaries of a CoolProp fluid file as the ancillary (initial-guess)
# equations, available through x0_psat and x0_sat_pure
function _saturation_superanc(pure, T)
    vl, vv = Clapeyron.x0_sat_pure(pure, T)
    return Clapeyron.x0_psat(pure, T), vl, vv
end
_saturation_solve(pure, T) = Clapeyron.saturation_pressure(pure, T)

# the ancillaries are superancillaries if they satisfy the phase equilibrium to near round-off
function _has_superancillary(pure, Tc)
    for T in (0.6Tc, 0.8Tc, 0.95Tc)
        p, vl, vv = _saturation_superanc(pure, T)
        all(isfinite, (p, vl, vv)) || return false
        gl = Clapeyron.VT_gibbs_free_energy(pure, vl, T)
        gv = Clapeyron.VT_gibbs_free_energy(pure, vv, T)
        pl, pv = Clapeyron.pressure(pure, vl, T), Clapeyron.pressure(pure, vv, T)
        (abs(pl / p - 1) < 1e-9 && abs(pv / p - 1) < 1e-9 && abs(gl - gv) / (Clapeyron.R̄ * T) < 1e-9) || return false
    end
    return true
end
