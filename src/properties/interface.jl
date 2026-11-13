#=
Refrigerant property interface.

The models only call the functions below, so any refrigerant backend that implements them can be
used: the spline-based `Refrigerant` (CoolProp sampling) or `ClapeyronRefrigerant` (Helmholtz
equations of state evaluated on the fly).

All quantities are SI and per kg: p [Pa], T [K], h [J/kg], s [J/kg/K], ρ [kg/m³], cp [J/kg/K],
μ [Pa s], λ [W/m/K]. Implementations must be generic in the number type so that
ForwardDiff.Dual numbers pass through.

Required by the evaporator models
---------------------------------
- `critical_pressure(r)`, `molar_mass(r)` [kg/mol]
- `pressure_bounds(r)` (not exported): `(pmin, pmax)` [Pa], admissible evaporation pressures;
  the models map the pressure unknown with `pressure_transform(r) = LogBox(pressure_bounds(r)...)`
- saturated vapour (dew point) at p: `T_dew`, `h_dew`, `ρ_dew`
- superheated vapour at p and superheat ΔT above the dew point: `h_superheated`,
  `ρ_superheated`, `cp_superheated`, `μ_superheated`, `λ_superheated`

Optional (backend dependent)
----------------------------
- bubble point: `T_bubble`, `h_bubble`, `ρ_bubble`; `s_dew`, `cp_dew`, `μ_dew`, `λ_dew`
- single phase at (p, T) with a phase hint: `h_pT`, `s_pT`, `ρ_pT`, `cp_pT`, `μ_pT`, `λ_pT`
- isentropic: `T_ps`, `h_ps`
=#

"""
    AbstractRefrigerant

Supertype of the refrigerant property backends ([`Refrigerant`](@ref),
[`ClapeyronRefrigerant`](@ref)). See `src/properties/interface.jl` for the functions a backend
implements.
"""
abstract type AbstractRefrigerant end

"Critical pressure [Pa] (of the mixture at its composition, for a mixture)."
critical_pressure(r::AbstractRefrigerant) = r.pcrit

"Molar mass [kg/mol]."
molar_mass(r::AbstractRefrigerant) = r.M

"""
    pressure_bounds(r) -> (pmin, pmax)

Interval [Pa] of admissible evaporation pressures, `0 < pmin < pmax`. For the spline
`Refrigerant` it is the sampled pressure range; for `ClapeyronRefrigerant` it is (triple-point
pressure, critical pressure).

Not exported: map the pressure with [`pressure_transform`](@ref) (`to_R`/`from_R`) and read the
interval with `bounds(pressure_transform(r))`. A backend extends
`OptimalGshpFlow.pressure_bounds`.
"""
function pressure_bounds end

"""
    pressure_transform(r) = LogBox(pressure_bounds(r)...)

Map of the evaporation-pressure unknown: `pe = from_R(pressure_transform(r), u)`,
`u = to_R(pressure_transform(r), pe)` (sigmoid in `log p`, see [`LogBox`](@ref)).
"""
pressure_transform(r::AbstractRefrigerant) = LogBox(pressure_bounds(r)...)

# saturation (dew = saturated vapour, bubble = saturated liquid)
function T_dew end
function T_bubble end
function h_dew end
function h_bubble end
function ρ_dew end
function ρ_bubble end
function s_dew end
function cp_dew end
function μ_dew end
function λ_dew end

# single phase at (p, T); keyword `phase = :unknown | :vapor | :liquid`
function h_pT end
function s_pT end
function ρ_pT end
function cp_pT end
function μ_pT end
function λ_pT end

"""
    h_superheated(r, p, ΔT)

Specific enthalpy [J/kg] of vapour at pressure `p` and `ΔT` kelvin above the dew point
temperature, i.e. at `(p, T_dew(r, p) + ΔT)`. Likewise `s_superheated`, `ρ_superheated`,
`cp_superheated`, `μ_superheated`, `λ_superheated`.

The default method evaluates the `(p, T)` function with `phase = :vapor`; the spline backend
interpolates directly in `(log10(p), ΔT)`.
"""
h_superheated(r::AbstractRefrigerant, p, ΔT) = h_pT(r, p, T_dew(r, p) + ΔT; phase = :vapor)
s_superheated(r::AbstractRefrigerant, p, ΔT) = s_pT(r, p, T_dew(r, p) + ΔT; phase = :vapor)
ρ_superheated(r::AbstractRefrigerant, p, ΔT) = ρ_pT(r, p, T_dew(r, p) + ΔT; phase = :vapor)
cp_superheated(r::AbstractRefrigerant, p, ΔT) = cp_pT(r, p, T_dew(r, p) + ΔT; phase = :vapor)
μ_superheated(r::AbstractRefrigerant, p, ΔT) = μ_pT(r, p, T_dew(r, p) + ΔT; phase = :vapor)
λ_superheated(r::AbstractRefrigerant, p, ΔT) = λ_pT(r, p, T_dew(r, p) + ΔT; phase = :vapor)
