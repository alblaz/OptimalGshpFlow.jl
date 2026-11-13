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
(co-volume, the lower bound of the molar volume [m³/mol], see [`volume_from_R`](@ref)).
"""
struct ClapeyronRefrigerant{E, TB <: TransportBackend, TM} <: AbstractRefrigerant
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
end

function ClapeyronRefrigerant(components::AbstractVector{<:AbstractString};
                              mass_fractions = nothing, mole_fractions = nothing,
                              transport::Union{Nothing, TransportBackend} = nothing,
                              name = join(components, "&"))
    comps = String.(components)
    n = length(comps)
    transport === nothing && (transport = n == 1 ? CoolPropTransport() : EntropyScalingTransport())
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

    return ClapeyronRefrigerant(name, comps, eos, transport, transport_model, collect(float.(z)), M, Tc, pc, vmin)
end
ClapeyronRefrigerant(component::AbstractString; kwargs...) =
    ClapeyronRefrigerant([component]; name = component, kwargs...)

ispure(r::ClapeyronRefrigerant) = length(r.components) == 1

# no lower bound: the saturation solver works down to the triple point
pressure_bounds(r::ClapeyronRefrigerant) = (zero(r.pcrit), r.pcrit)

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
        T, vl, vv = Clapeyron.saturation_temperature(r.eos, p)
    else
        T, vl, vv, _ = Clapeyron.dew_temperature(r.eos, p, r.z)
    end
    return T, vv
end

function saturation_bubble(r::ClapeyronRefrigerant, p)
    if ispure(r)
        T, vl, vv = Clapeyron.saturation_temperature(r.eos, p)
    else
        T, vl, vv, _ = Clapeyron.bubble_temperature(r.eos, p, r.z)
    end
    return T, vl
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

"Temperature [K] at pressure `p` and specific entropy `s` [J/kg/K] (p-s flash). `T0`: optional initial guess."
T_ps(r::ClapeyronRefrigerant, p, s; phase = :unknown, T0 = nothing) =
    T0 === nothing ? Clapeyron.PS.temperature(r.eos, p, s * r.M, r.z; phase = phase) :
                     Clapeyron.PS.temperature(r.eos, p, s * r.M, r.z; phase = phase, T0 = T0)
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

"""
    CoolPropState(fluid)

Cached CoolProp `AbstractState` (HEOS backend) of a pure fluid, used for transport properties.
Updating it and reading outputs costs ≈ 3 µs, against ≈ 100 µs for one `PropsSI` call, with
identical results. A lock serialises access, so a refrigerant can be shared between threads.
"""
mutable struct CoolPropState
    fluid::String
    handle::Clong
    lock::ReentrantLock
    function CoolPropState(fluid::AbstractString)
        cs = new(String(fluid), CoolProp.AbstractState_factory("HEOS", fluid), ReentrantLock())
        finalizer(x -> CoolProp.AbstractState_free(x.handle), cs)
        return cs
    end
end
Base.show(io::IO, cs::CoolPropState) = print(io, "CoolPropState(", cs.fluid, ")")

const COOLPROP_PT = Ref{Clong}(-1)
const COOLPROP_PQ = Ref{Clong}(-1)
const COOLPROP_DT = Ref{Clong}(-1)
const COOLPROP_V  = Ref{Clong}(-1)
const COOLPROP_L  = Ref{Clong}(-1)

function _init_coolprop_indices()
    COOLPROP_PT[] = CoolProp.get_input_pair_index("PT_INPUTS")
    COOLPROP_PQ[] = CoolProp.get_input_pair_index("PQ_INPUTS")
    COOLPROP_DT[] = CoolProp.get_input_pair_index("DmassT_INPUTS")
    COOLPROP_V[]  = CoolProp.get_param_index("V")
    COOLPROP_L[]  = CoolProp.get_param_index("CONDUCTIVITY")
    return nothing
end

_coolprop_key(prop) = prop === :viscosity ? COOLPROP_V : COOLPROP_L

function _coolprop_value(cs::CoolPropState, pair, key, x1, x2)
    lock(cs.lock) do
        CoolProp.AbstractState_update(cs.handle, pair[], x1, x2)
        CoolProp.AbstractState_keyed_output(cs.handle, key[])
    end
end

# CoolProp output for the input pair (x1, x2); for ForwardDiff.Dual inputs the derivatives are
# central differences with relative steps rel1, rel2
function _coolprop(cs::CoolPropState, pair, key, x1, x2, rel1, rel2)
    (x1 isa ForwardDiff.Dual || x2 isa ForwardDiff.Dual) || return _coolprop_value(cs, pair, key, x1, x2)
    a, b = ForwardDiff.value(x1), ForwardDiff.value(x2)
    f = _coolprop_value(cs, pair, key, a, b)
    D = x1 isa ForwardDiff.Dual ? typeof(x1) : typeof(x2)
    ∂ = zero(ForwardDiff.partials(one(D)))
    if x1 isa ForwardDiff.Dual
        h = rel1 * a
        ∂ += (_coolprop_value(cs, pair, key, a + h, b) - _coolprop_value(cs, pair, key, a - h, b)) / (2h) * ForwardDiff.partials(x1)
    end
    if x2 isa ForwardDiff.Dual
        h = rel2 * b
        ∂ += (_coolprop_value(cs, pair, key, a, b + h) - _coolprop_value(cs, pair, key, a, b - h)) / (2h) * ForwardDiff.partials(x2)
    end
    return D(f, ∂)
end
