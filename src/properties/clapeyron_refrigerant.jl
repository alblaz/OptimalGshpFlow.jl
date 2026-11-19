#=
Refrigerant properties computed on the fly with Clapeyron.jl (thermodynamics) and
EntropyScaling.jl or CoolProp (transport), as an alternative to the spline-based `Refrigerant`.

All functions take and return SI, mass-based quantities: p [Pa], T [K], h [J/kg], s [J/kg/K],
ρ [kg/m³], cp [J/kg/K], μ [Pa s], λ [W/m/K]. They are generic in the number type so that
ForwardDiff.Dual numbers pass through (AD-ready for the nonlinear solvers).

Clapeyron reads the reference Helmholtz equations of state from the CoolProp fluid library
(loaded through the CoolProp dependency), so thermodynamic properties of pure fluids reproduce
CoolProp to round-off.
=#

########################
# transport backends   #
########################

"""
    TransportBackend

How viscosity and thermal conductivity are computed for a [`ClapeyronRefrigerant`](@ref).

- `CoolPropTransport()` (default for pure fluids): CoolProp reference correlations. Pure fluids
  only. Derivatives for AD are computed by central finite differences.
- `EntropyScalingTransport()` (default for mixtures): EntropyScaling.jl `RefpropRES` model. Works
  for mixtures and with automatic differentiation. Accuracy vs CoolProp (pure fluids, vapour,
  reduced pressure 0.02–0.7): viscosity within ≈ 5–14 %, thermal conductivity within ≈ 6–17 %
  (R290: +17 %); liquid within ≈ 4 %.
"""
abstract type TransportBackend end
struct EntropyScalingTransport <: TransportBackend end
struct CoolPropTransport <: TransportBackend end

# names used by EntropyScaling for refrigerants known by their ASHRAE number
const ENTROPYSCALING_NAMES = Dict(
    "R290" => "propane", "R600" => "butane", "R600a" => "isobutane", "R1270" => "propene",
    "R717" => "ammonia", "R744" => "carbon dioxide", "R718" => "water", "R50" => "methane",
    "R170" => "ethane")

########################
# refrigerant type     #
########################

"Composition of a [`ClapeyronRefrigerant`](@ref), a type parameter: `Pure` or `Mixture`."
abstract type Composition end
struct Pure <: Composition end
struct Mixture <: Composition end

"""
    ClapeyronRefrigerant(components; mass_fractions = nothing, mole_fractions = nothing,
                         transport = nothing, name = ...)

Refrigerant whose properties are computed on the fly with Clapeyron.jl (multi-parameter
Helmholtz equations of state, `MultiFluid`).

- `components`: a fluid name (`"R290"`, `"R32"`, ...) or a vector of names for a mixture.
- `mass_fractions` / `mole_fractions`: composition of a mixture (one of them is required for
  mixtures). Refrigerant blends are usually specified by mass, e.g. R454C:
  `ClapeyronRefrigerant(["R32", "R1234yf"]; mass_fractions = [0.215, 0.785])`.
- `transport`: [`TransportBackend`](@ref). Default (`nothing`): `CoolPropTransport()` for a pure
  fluid, `EntropyScalingTransport()` for a mixture.

Performance: pass `phase = :vapor` (or `:liquid`) to the `(p, T)` functions whenever the phase is
known. With the default `phase = :unknown` both volume roots are computed and compared, which is
≈ 5–20× slower.

Fields: `eos` (Clapeyron model), `transport_model`, `z` (mole fractions), `M` (molar mass, kg/mol),
`Tcrit`, `pcrit` (critical point of the fluid or of the mixture at composition `z`), `vmin`
(co-volume, the lower bound of the molar volume [m³/mol], see [`volume_from_R`](@ref)),
`Ttriple`, `ptriple` (triple point; for a mixture `Ttriple = NaN` and `ptriple` is the lowest
triple-point pressure of the components), `superancillary` (pure fluid whose saturation curve is an explicit function of T, see
[`saturation_T`](@ref)).
"""
struct ClapeyronRefrigerant{E, TB <: TransportBackend, TM, C <: Composition} <: AbstractRefrigerant
    name::String
    components::Vector{String}
    eos::E
    transport::TB
    transport_model::TM
    z::Vector{Float64}
    M::Float64
    Tcrit::Float64
    pcrit::Float64
    vmin::Float64
    Ttriple::Float64
    ptriple::Float64
    superancillary::Bool
end

function ClapeyronRefrigerant(components::AbstractVector{<:AbstractString};
                              mass_fractions = nothing, mole_fractions = nothing,
                              transport::Union{Nothing, TransportBackend} = nothing,
                              name = join(components, "&"))
    comps = String.(components)
    if length(comps) == 1
        return _clapeyron_refrigerant(comps, something(transport, CoolPropTransport()), Pure(), mass_fractions,
                                      mole_fractions, String(name))
    else
        return _clapeyron_refrigerant(comps, something(transport, EntropyScalingTransport()), Mixture(),
                                      mass_fractions, mole_fractions, String(name))
    end
end
# a single fluid: purity and the default transport are known statically, so the type is inferred
ClapeyronRefrigerant(component::AbstractString; transport::TransportBackend = CoolPropTransport(),
                     name = component, mass_fractions = nothing, mole_fractions = nothing) =
    _clapeyron_refrigerant([String(component)], transport, Pure(), mass_fractions, mole_fractions, String(name))

# function barrier: specialized on the transport backend and the composition, so that the type of
# the result (transport model, composition) is inferred
function _clapeyron_refrigerant(comps::Vector{String}, transport::TransportBackend, composition::C, mass_fractions,
                                mole_fractions, name::String) where {C <: Composition}
    n = length(comps)
    eos = Clapeyron.MultiFluid(comps)
    Mi = Clapeyron.mw(eos) ./ 1000                     # kg/mol

    z = if n == 1
        [1.0]
    elseif mole_fractions !== nothing
        mole_fractions ./ sum(mole_fractions)
    elseif mass_fractions !== nothing
        x = mass_fractions ./ Mi
        x ./ sum(x)
    else
        throw(ArgumentError("a mixture needs `mass_fractions` or `mole_fractions`"))
    end
    length(z) == n || throw(ArgumentError("composition has $(length(z)) entries for $n components"))

    M = sum(z .* Mi)
    Tc, pc, _ = n == 1 ? Clapeyron.crit_pure(eos) : Clapeyron.crit_mix(eos, z)

    transport_model = if transport isa EntropyScalingTransport
        EntropyScaling.RefpropRES([get(ENTROPYSCALING_NAMES, c, c) for c in comps])
    else
        n == 1 || throw(ArgumentError("CoolPropTransport supports pure fluids only"))
        CoolPropState(comps[1])
    end
    vmin = Clapeyron.lb_volume(eos, Tc, z)   # independent of T for MultiFluid
    Ttriple = n == 1 ? Float64(eos.pures[1].properties.Ttp) : NaN
    ptriple = minimum(Float64(pure.properties.ptp) for pure in eos.pures)
    superanc = n == 1 && _has_superancillary(eos.pures[1], Ttriple, Tc)

    return ClapeyronRefrigerant{typeof(eos), typeof(transport), typeof(transport_model), C}(name, comps, eos,
               transport, transport_model, collect(float.(z)), M, Tc, pc, vmin, Ttriple, ptriple, superanc)
end

"`true` for a pure fluid: from the type, so that branches on it are resolved at compile time."
ispure(::ClapeyronRefrigerant{<:Any, <:Any, <:Any, Pure}) = true
ispure(::ClapeyronRefrigerant{<:Any, <:Any, <:Any, Mixture}) = false

pressure_bounds(r::ClapeyronRefrigerant) = (r.ptriple, r.pcrit)

Base.show(io::IO, r::ClapeyronRefrigerant) =
    print(io, "ClapeyronRefrigerant(", r.name, ", z = ", r.z, ", transport = ", nameof(typeof(r.transport)), ")")

########################
# saturation           #
########################

"""
    saturation_dew(r, p)    -> (T, v)
    saturation_bubble(r, p) -> (T, v)

Dew (saturated vapour) and bubble (saturated liquid) point at pressure `p`: temperature [K] and
molar volume [m³/mol] of the phase with the overall composition. For a pure fluid dew and
bubble temperatures coincide; for a zeotropic mixture `T_dew - T_bubble` is the temperature glide.
"""
function saturation_dew(r::ClapeyronRefrigerant, p)
    if ispure(r)
        T, vl, vv = _saturation_temperature_pure(r, p)
    else
        T, vl, vv, _ = Clapeyron.dew_temperature(r.eos, p, r.z)
    end
    return T, vv
end

function saturation_bubble(r::ClapeyronRefrigerant, p)
    if ispure(r)
        T, vl, vv = _saturation_temperature_pure(r, p)
    else
        T, vl, vv, _ = Clapeyron.bubble_temperature(r.eos, p, r.z)
    end
    return T, vl
end

# Pure fluid at pressure p: invert the superancillary p_sat(T) when available (explicit and
# exact, see saturation_T), otherwise Clapeyron's solver. Clapeyron 0.6.30
# saturation_temperature returns NaN for R744 at every pressure.
function _saturation_temperature_pure(r::ClapeyronRefrigerant, p)
    r.superancillary || return Clapeyron.saturation_temperature(r.eos, p)
    T = _T_sat_superanc(r, p)
    _, vl, vv = saturation_T(r, T)
    return T, vl, vv
end

_primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))
_primal(x) = x

# Newton on ln p_sat(T) = ln p with Float64 values, then one Newton step with the (possibly Dual)
# pressure: the step propagates the exact derivative dT/dp = 1/(p d ln p_sat/dT).
# NaN outside (p_triple, p_crit).
function _T_sat_superanc(r::ClapeyronRefrigerant, p)
    lnp = log(_primal(p))
    lnpsat(T) = log(p_sat(r, T))
    Ttp, Tc = r.Ttriple, r.Tcrit
    lnptp, lnpc = lnpsat(Ttp + 1e-6 * (Tc - Ttp)), log(r.pcrit)
    lnptp < lnp < lnpc || return oftype(float(_primal(p)), NaN) * one(p)
    # initial guess: ln p linear in 1/T between the triple and critical points
    T = 1 / (1/Ttp + (lnp - lnptp) / (lnpc - lnptp) * (1/Tc - 1/Ttp))
    for _ in 1:30
        ΔT = (lnpsat(T) - lnp) / ForwardDiff.derivative(lnpsat, T)
        T = clamp(T - ΔT, Ttp + 1e-9, Tc - 1e-9)
        abs(ΔT) < 1e-11 * T && break
    end
    return T - (lnpsat(T) - log(p)) / ForwardDiff.derivative(lnpsat, T)
end

"Dew point temperature [K] at pressure `p` [Pa]."
T_dew(r::ClapeyronRefrigerant, p) = saturation_dew(r, p)[1]
"Bubble point temperature [K] at pressure `p` [Pa]."
T_bubble(r::ClapeyronRefrigerant, p) = saturation_bubble(r, p)[1]

"Specific enthalpy [J/kg] of saturated vapour (dew point) at pressure `p`."
h_dew(r::ClapeyronRefrigerant, p) = ((T, v) = saturation_dew(r, p); Clapeyron.VT_enthalpy(r.eos, v, T, r.z) / r.M)
"Specific enthalpy [J/kg] of saturated liquid (bubble point) at pressure `p`."
h_bubble(r::ClapeyronRefrigerant, p) = ((T, v) = saturation_bubble(r, p); Clapeyron.VT_enthalpy(r.eos, v, T, r.z) / r.M)

"Density [kg/m³] of saturated vapour (dew point) at pressure `p`."
ρ_dew(r::ClapeyronRefrigerant, p) = r.M / saturation_dew(r, p)[2]
"Density [kg/m³] of saturated liquid (bubble point) at pressure `p`."
ρ_bubble(r::ClapeyronRefrigerant, p) = r.M / saturation_bubble(r, p)[2]

"Specific entropy [J/kg/K] of saturated vapour (dew point) at pressure `p`."
s_dew(r::ClapeyronRefrigerant, p) = ((T, v) = saturation_dew(r, p); Clapeyron.VT_entropy(r.eos, v, T, r.z) / r.M)

"Isobaric specific heat [J/kg/K] of saturated vapour (dew point) at pressure `p`."
cp_dew(r::ClapeyronRefrigerant, p) = ((T, v) = saturation_dew(r, p); Clapeyron.VT_isobaric_heat_capacity(r.eos, v, T, r.z) / r.M)

"Viscosity [Pa s] of saturated vapour (dew point) at pressure `p`."
μ_dew(r::ClapeyronRefrigerant, p) = _transport_dew(r, :viscosity, p)
"Thermal conductivity [W/m/K] of saturated vapour (dew point) at pressure `p`."
λ_dew(r::ClapeyronRefrigerant, p) = _transport_dew(r, :conductivity, p)

"""
    cp_twophase(r, p)

Pseudo two-phase specific heat `(h_dew - h_bubble)/(T_dew - T_bubble)` [J/kg/K] of a zeotropic
mixture at pressure `p` (Qiao 2025, eq. 3). Infinite for a pure fluid (no glide).
"""
function cp_twophase(r::ClapeyronRefrigerant, p)
    ispure(r) && return Inf
    return (h_dew(r, p) - h_bubble(r, p)) / (T_dew(r, p) - T_bubble(r, p))
end

"""
    T_ph(r, p, h)

Temperature [K] at pressure `p` and specific enthalpy `h` [J/kg] (Clapeyron p-h flash): in the
two-phase region of a zeotropic mixture it lies between `T_bubble(r, p)` and `T_dew(r, p)`. Used
for the refrigerant temperature at the evaporator inlet (glide). An iterative flash: ≈ 0.5–0.7 ms
for a binary mixture, about as much as all the saturation calls of an evaporator residual.
Clapeyron's flash is not inferable: the result is asserted to the promoted type of `p` and `h`
(an inference barrier; without it everything downstream, e.g. the evaporator residual, infers `Any`).
"""
T_ph(r::ClapeyronRefrigerant, p, h) = Clapeyron.PH.temperature(r.eos, p, h * r.M, r.z)::float(promote_type(typeof(p), typeof(h)))

########################
# single phase (p, T)  #
########################

"Specific enthalpy [J/kg] at `(p, T)`. `phase` (`:unknown`, `:vapor`, `:liquid`) selects the volume root."
h_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = Clapeyron.enthalpy(r.eos, p, T, r.z; phase = phase) / r.M
"Specific entropy [J/kg/K] at `(p, T)`."
s_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = Clapeyron.entropy(r.eos, p, T, r.z; phase = phase) / r.M
"Density [kg/m³] at `(p, T)`."
ρ_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = Clapeyron.mass_density(r.eos, p, T, r.z; phase = phase)
"Isobaric specific heat [J/kg/K] at `(p, T)`."
cp_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = Clapeyron.isobaric_heat_capacity(r.eos, p, T, r.z; phase = phase) / r.M
"Viscosity [Pa s] at `(p, T)`."
μ_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = _transport(r, :viscosity, p, T, phase)
"Thermal conductivity [W/m/K] at `(p, T)`."
λ_pT(r::ClapeyronRefrigerant, p, T; phase = :unknown) = _transport(r, :conductivity, p, T, phase)

########################
# isentropic           #
########################

"""
    T_ps(r, p, s; phase = :unknown, T0 = nothing)

Temperature [K] at pressure `p` and specific entropy `s` [J/kg/K]. `T0`: optional initial guess.

With `phase = :vapor`, a vapour state below the critical pressure (`s ≥ s_dew(p)`, e.g. an
isentropic compressor outlet) is solved by Newton in `T` at constant `p` (`ds = cp dT / T`), which
reuses the single-phase functions of the model; any other case goes to Clapeyron's p-s flash.
"""
function T_ps(r::ClapeyronRefrigerant, p, s; phase = :unknown, T0 = nothing)
    if phase === :vapor
        T = _T_ps_vapor(r, p, s, T0)
        T === nothing || return T
    end
    # Clapeyron's flash behind an inference barrier: inferring it costs ≈ 16 s of compilation, paid
    # only when it is actually called; its result (inferred `Any`) is converted to a concrete type
    Tf = Base.inferencebarrier(_T_ps_flash)(r, p, s, phase, T0)
    TT = promote_type(typeof(p), typeof(s), Float64)
    return convert(TT, Tf)::TT
end

_T_ps_flash(r, p, s, phase, T0) =
    T0 === nothing ? Clapeyron.PS.temperature(r.eos, p, s * r.M, r.z; phase = phase) :
                     Clapeyron.PS.temperature(r.eos, p, s * r.M, r.z; phase = phase, T0 = T0)

# Newton in T at constant p for a vapour state; `nothing` if the state is not a subcritical vapour
# or if Newton does not converge (the caller then uses Clapeyron's flash). Iterates on the primal
# values; the final step, with the Duals of `p` and `s`, gives the exact derivatives (implicit
# function theorem: ∂T/∂s = T/cp, ∂T/∂p = -T/cp ∂s/∂p).
function _T_ps_vapor(r::ClapeyronRefrigerant, p, s, T0)
    pp, ss = _primal(p), _primal(s)
    pp < r.pcrit || return nothing
    Td = T_dew(r, pp)
    (isfinite(Td) && ss >= s_dew(r, pp)) || return nothing
    newton_step(p, s, T) = (s_pT(r, p, T; phase = :vapor) - s) * T / cp_pT(r, p, T; phase = :vapor)
    # initial guess: constant cp from the dew point
    T = T0 === nothing ? Td * exp((ss - s_dew(r, pp)) / cp_pT(r, pp, Td; phase = :vapor)) : float(_primal(T0))
    for _ in 1:50
        ΔT = newton_step(pp, ss, T)
        T = max(T - ΔT, Td)
        abs(ΔT) < 1e-10 * T && return T - newton_step(p, s, T)
    end
    return nothing
end
"Specific enthalpy [J/kg] at pressure `p` and specific entropy `s` (e.g. isentropic compressor outlet: `phase = :vapor`)."
h_ps(r::ClapeyronRefrigerant, p, s; phase = :unknown, T0 = nothing) = h_pT(r, p, T_ps(r, p, s; phase = phase, T0 = T0); phase = phase)

########################
# transport internals  #
########################

function _transport(r::ClapeyronRefrigerant{<:Any, EntropyScalingTransport}, prop, p, T, phase)
    if prop === :viscosity
        return EntropyScaling.viscosity(r.transport_model, p, T, r.z; phase = phase)
    else
        return EntropyScaling.thermal_conductivity(r.transport_model, p, T, r.z; phase = phase)
    end
end

function _transport(r::ClapeyronRefrigerant{<:Any, CoolPropTransport}, prop, p, T, phase)
    return _coolprop(r.transport_model, COOLPROP_PT, _coolprop_key(prop), p, T, 1e-6, 1e-4)
end

_transport_dew(r::ClapeyronRefrigerant{<:Any, EntropyScalingTransport}, prop, p) =
    _transport(r, prop, p, T_dew(r, p), :vapor)

# saturated vapour with CoolProp: (p, Q = 1) avoids evaluating exactly on the phase boundary
_transport_dew(r::ClapeyronRefrigerant{<:Any, CoolPropTransport}, prop, p) =
    _coolprop(r.transport_model, COOLPROP_PQ, _coolprop_key(prop), p, 1.0, 1e-6, 0.0)

########################
# CoolProp low level   #
########################

_coolprop_key(prop) = prop === :viscosity ? COOLPROP_V : COOLPROP_L

# surface tension: CoolProp (pure fluids with CoolProp transport)
function surface_tension(r::ClapeyronRefrigerant, p)
    cs = r.transport_model
    cs isa CoolPropState || throw(ArgumentError("surface_tension needs a pure fluid with CoolPropTransport"))
    return _coolprop(cs, COOLPROP_PQ, COOLPROP_I, p, 0.0, 1e-6, 0.0)
end
