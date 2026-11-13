#=
How the evaporator models obtain the refrigerant state (feature/property-interface, sub-step 4).

The heat balances (`evaporator_residuals!`) are the same for every formulation; only the first
unknown and, for `ExplicitTv`, extra (T, v) unknowns differ:

| formulation               | u[1]          | extra unknowns (after the base ones) | inner iterations        |
|---------------------------|---------------|--------------------------------------|-------------------------|
| `NestedPressure()`        | pressure      | none                                 | dew point, (p,T) volume |
| `SaturationTemperature()` | dew temp. Te  | none                                 | 2 (p,T) volume solves   |
| `ExplicitTv()`            | dew temp. Te  | volumes at Te+ΔTsup and Te+ΔTsup/2   | none                    |

`NestedPressure` works with every `AbstractRefrigerant`; the other two need a pure
`ClapeyronRefrigerant` (explicit saturation from the superancillaries, see `saturation_T`).
=#

"""
    PropertyFormulation

How the evaporator models obtain the refrigerant state: [`NestedPressure`](@ref) (default),
[`SaturationTemperature`](@ref) or [`ExplicitTv`](@ref). Pass it as last argument, e.g.
`evaporator_model!(du, u, p, ExplicitTv())`.
"""
abstract type PropertyFormulation end

"""
    NestedPressure()

Unknown `u[1]`: evaporation pressure (sigmoid on [`pressure_bounds`](@ref)). Properties come from
the refrigerant interface at (p, ΔT), with inner iterations for the Clapeyron backend. The
original formulation, valid for every refrigerant backend.
"""
struct NestedPressure <: PropertyFormulation end

"""
    SaturationTemperature()

Unknown `u[1]`: dew (evaporation) temperature `Te` (sigmoid on `(Ttriple, Tcrit)`). Saturation is
explicit ([`saturation_T`](@ref)); the superheated states need one (p, T) volume solve each
(outlet at `Te + ΔTsup`, mean at `Te + ΔTsup/2`). Pure `ClapeyronRefrigerant` only.
"""
struct SaturationTemperature <: PropertyFormulation end

"""
    ExplicitTv()

As [`SaturationTemperature`](@ref), plus two unknowns appended after the base ones: the molar
volumes (mapped with [`volume_to_R`](@ref)) of the superheated vapour at the outlet
`Te + ΔTsup` and at the mean superheat `Te + ΔTsup/2`, closed by the residuals
`pressure_residual(r, T, v, pe)`. All properties are explicit; no inner iteration.
Pure `ClapeyronRefrigerant` only.
"""
struct ExplicitTv <: PropertyFormulation end

"Number of unknowns a formulation appends to the base unknowns of the model."
n_extra_unknowns(::PropertyFormulation) = 0
n_extra_unknowns(::ExplicitTv) = 2

"Interval `(Tmin, Tmax)` [K] onto which the evaporation-temperature unknown is mapped."
temperature_bounds(r::ClapeyronRefrigerant) = (r.Ttriple, r.Tcrit)

"""
    evaporator_refrigerant_state!(du, formulation, refrigerant, u, nbase, ΔTsup)

Refrigerant state used by `evaporator_residuals!`, as a named tuple
`(pe, Te, hes, hsup, ρsup, cpm, μm, λm)`: pressure, dew temperature, dew-point enthalpy, outlet
enthalpy and density (at `ΔTsup` superheat), and cp, μ, λ at the mean superheat `ΔTsup/2`.
`nbase` is the number of base unknowns of the model; a formulation with extra unknowns writes
their residuals into `du[nbase+1:end]`.
"""
function evaporator_refrigerant_state!(du, ::NestedPressure, r::AbstractRefrigerant, u, nbase, ΔTsup)
    pe = sigmoid_R_to_ab(u[1], pressure_bounds(r)...)
    return (pe = pe, Te = T_dew(r, pe), hes = h_dew(r, pe),
            hsup = h_superheated(r, pe, ΔTsup), ρsup = ρ_superheated(r, pe, ΔTsup),
            cpm = cp_superheated(r, pe, ΔTsup/2), μm = μ_superheated(r, pe, ΔTsup/2), λm = λ_superheated(r, pe, ΔTsup/2))
end

function evaporator_refrigerant_state!(du, ::SaturationTemperature, r::ClapeyronRefrigerant, u, nbase, ΔTsup)
    Te = sigmoid_R_to_ab(u[1], temperature_bounds(r)...)
    pe, _, vv = saturation_T(r, Te)
    Tsup, Tm = Te + ΔTsup, Te + ΔTsup/2
    vsup = volume_pT(r, pe, Tsup; phase = :vapor)
    vm   = volume_pT(r, pe, Tm; phase = :vapor)
    return _superheat_state(r, pe, Te, vv, Tsup, vsup, Tm, vm)
end

function evaporator_refrigerant_state!(du, ::ExplicitTv, r::ClapeyronRefrigerant, u, nbase, ΔTsup)
    Te = sigmoid_R_to_ab(u[1], temperature_bounds(r)...)
    pe, _, vv = saturation_T(r, Te)
    Tsup, Tm = Te + ΔTsup, Te + ΔTsup/2
    vsup = volume_from_R(r, u[nbase + 1])
    vm   = volume_from_R(r, u[nbase + 2])
    du[nbase + 1] = pressure_residual(r, Tsup, vsup, pe)
    du[nbase + 2] = pressure_residual(r, Tm, vm, pe)
    return _superheat_state(r, pe, Te, vv, Tsup, vsup, Tm, vm)
end

_superheat_state(r, pe, Te, vv, Tsup, vsup, Tm, vm) =
    (pe = pe, Te = Te, hes = h_Tv(r, Te, vv), hsup = h_Tv(r, Tsup, vsup), ρsup = ρ_Tv(r, Tsup, vsup),
     cpm = cp_Tv(r, Tm, vm), μm = μ_Tv(r, Tm, vm), λm = λ_Tv(r, Tm, vm))

########################
# unknowns             #
########################

# dew pressure at temperature T, for building initial guesses
_p_dew(r::Refrigerant, T) = PropsSI("P", "T", T, "Q", 1, r.name)
_p_dew(r::ClapeyronRefrigerant, T) = ispure(r) ? p_sat(r, T) : Clapeyron.dew_pressure(r.eos, T, r.z)[1]

"""
    evaporator_unknowns(formulation, refrigerant, ΔTsup; Te, Qe, Qsup, Ts_l, Ts_m,
                        ϵe = 0.5, ϵsup = 0.5, xs = 0.5, Ts_h = nothing)

Vector of unknowns for [`evaporator_model!`](@ref) (or, with `Ts_h`, for
[`evaporatorborehole_model!`](@ref)) in the given formulation, from physical values: dew
temperature `Te`, heat flows, effectivenesses, moving-boundary position and secondary-fluid
temperatures.
"""
function evaporator_unknowns(formulation::PropertyFormulation, r, ΔTsup; Te, Qe, Qsup, Ts_l, Ts_m,
                             ϵe = 0.5, ϵsup = 0.5, xs = 0.5, Ts_h = nothing)
    u1 = formulation isa NestedPressure ? sigmoid_ab_to_R(_p_dew(r, Te), pressure_bounds(r)...) :
                                          sigmoid_ab_to_R(Te, temperature_bounds(r)...)
    u = [u1, Qe, Qsup, sigmoid_01_to_R(ϵe), sigmoid_01_to_R(ϵsup), sigmoid_01_to_R(xs), Ts_l, Ts_m]
    Ts_h === nothing || push!(u, Ts_h)
    if formulation isa ExplicitTv
        pe = _p_dew(r, Te)
        push!(u, volume_to_R(r, volume_pT(r, pe, Te + ΔTsup; phase = :vapor)),
                 volume_to_R(r, volume_pT(r, pe, Te + ΔTsup/2; phase = :vapor)))
    end
    return u
end

"""
    nested_unknowns(formulation, refrigerant, u, nbase)

Convert the unknowns `u` of a formulation to those of `NestedPressure()` (drop the extra
unknowns, map `u[1]` to the pressure), e.g. to post-process with [`evaporator_model`](@ref).
`nbase` is 8 for the evaporator and 9 for the evaporator coupled to the borehole.
"""
nested_unknowns(::NestedPressure, r, u, nbase) = u[1:nbase]
function nested_unknowns(::Union{SaturationTemperature, ExplicitTv}, r::ClapeyronRefrigerant, u, nbase)
    Te = sigmoid_R_to_ab(u[1], temperature_bounds(r)...)
    v = u[1:nbase]
    v[1] = sigmoid_ab_to_R(p_sat(r, Te), pressure_bounds(r)...)
    return v
end
