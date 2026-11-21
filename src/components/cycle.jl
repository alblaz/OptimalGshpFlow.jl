#=
Closed cycle (step 9d): port transforms shared by the components of the refrigerant loop, and the
compressor as a component. The loop: evaporator (outlet `p_suc`, `T_suc`, mass flow `m`) →
compressor (discharge enthalpy `h_dis` at the condensing pressure `pc`) → condenser (outlet
`h_out` = saturated liquid at `pc`) → expansion valve (isenthalpic: the condenser outlet node is
the evaporator inlet `hei`) → evaporator.
=#

"""
    CyclePorts(refrigerant; Te = 273.15)

Transforms of the refrigerant-loop ports, the same in every component built for `refrigerant` (all
ports of a node must have the same transform): pressure `p` ([`pressure_transform`](@ref)),
enthalpy `h` (`Linear` from the saturated liquid at `Te`, in units of the latent heat), mass flow
`m` (`Above(0, 1)`: log scale, independent of the machine size). `pscale` (dew pressure at `Te`)
and `hscale` (latent heat at `Te`) scale the residuals of these ports.
"""
struct CyclePorts{T <: Real}
    p::LogBox{T}
    h::Linear{T}
    m::Above{T}
    pscale::T
    hscale::T
end

function CyclePorts(r::AbstractRefrigerant; Te = 273.15)
    p0 = _p_dew(r, Te)
    hl, Δh = h_bubble(r, p0), h_dew(r, p0) - h_bubble(r, p0)
    hl, Δh, p0 = Float64(hl), Float64(Δh), Float64(p0)
    return CyclePorts(pressure_transform(r), Linear(hl, Δh), Above(0.0, 1.0), p0, Δh)
end

"""
    vapour_h_ps(r, p, s) -> h

Enthalpy of superheated vapour at pressure `p` and entropy `s` (≥ the dew entropy), generic in the
number type: [`h_ps`](@ref) (`phase = :vapor`) for a `ClapeyronRefrigerant`, Newton in T on
`s_pT` (`∂s/∂T = cp/T`) from the dew point otherwise.
"""
vapour_h_ps(r::ClapeyronRefrigerant, p, s) = h_ps(r, p, s; phase = :vapor)
function vapour_h_ps(r::AbstractRefrigerant, p, s)
    Td = T_dew(r, p)
    T = Td * exp((s - s_dew(r, p)) / cp_dew(r, p))
    for _ in 1:4
        T -= (s_pT(r, p, T; phase = :vapor) - s) * T / cp_pT(r, p, T; phase = :vapor)
        T = max(T, Td)
    end
    return h_pT(r, p, T; phase = :vapor)
end

"""
    compression_state(r, pe, Tsuc, Δp, pc, ηis) -> (; h_in, h_out)

As [`compression`](@ref), on the property interface for every refrigerant (generic in the number
type, so the compressor can sit inside a residual; the spline `Refrigerant` goes through its
superheat table instead of CoolProp, see [`vapour_h_ps`](@ref)).
"""
function compression_state(r::AbstractRefrigerant, pe, Tsuc, Δp, pc, ηis)
    h_in = h_pT(r, pe, Tsuc; phase = :vapor)
    # adiabatic suction pressure drop: isenthalpic to pe − Δp (Newton in T, two steps)
    ps, T = pe - Δp, Tsuc
    for _ in 1:2
        T -= (h_pT(r, ps, T; phase = :vapor) - h_in) / cp_pT(r, ps, T; phase = :vapor)
    end
    h_is = vapour_h_ps(r, pc, s_pT(r, ps, T; phase = :vapor))
    return (h_in = h_in, h_out = h_in + (h_is - h_in) / ηis)
end

"""
    Compressor(; refrigerant, ηis, ξ = 0.07, Δp_suction = 0.0)

Compressor of a closed cycle, isentropic efficiency `ηis`. The displacement is the `Vref` of the
[`Evaporator`](@ref) (the mass flow `m` is its `ρ_sup Vref`, see B11); `ξ` is the fraction of the
shaft power lost through the envelope (for the power only); `Δp_suction` an adiabatic suction
pressure drop (see [`compressor_power`](@ref)). Volumetric and inverter efficiencies: not yet.

- Ports: `p_suc`, `T_suc` (evaporator outlet), `m`, `pc` (condensing pressure) and `h_dis`
  (discharge enthalpy). Transforms from [`CyclePorts`](@ref).
- No internal unknowns; one residual `h_dis − h_out(p_suc, T_suc, pc)`, scaled by the latent heat.

Power: [`compressor_power`](@ref)`(c::Compressor, p)`.
"""
struct Compressor{R <: AbstractRefrigerant, T <: Real} <: AbstractComponent
    refrigerant::R
    ηis::T
    ξ::T
    Δp_suction::T
    cycleports::CyclePorts{T}
end

function Compressor(; refrigerant, ηis, ξ = 0.07, Δp_suction = 0.0)
    cp = CyclePorts(refrigerant)
    T = typeof(cp.pscale)
    return Compressor(refrigerant, T(ηis), T(ξ), T(Δp_suction), cp)
end

function ports(c::Compressor)
    cp = c.cycleports
    return (p_suc = cp.p, T_suc = Free(), m = cp.m, pc = cp.p, h_dis = cp.h)
end
residual_scales(c::Compressor) = (c.cycleports.hscale,)

_compression_state(c::Compressor, p) = compression_state(c.refrigerant, p.p_suc, p.T_suc, c.Δp_suction, p.pc, c.ηis)

function residuals!(r, c::Compressor, x, p)
    r[1] = p.h_dis - _compression_state(c, p).h_out
    return nothing
end

"""
    compressor_power(c::Compressor, p) -> W

Electric power of the compressor [W] at the port values `p` (e.g. from [`port_values`](@ref)):
`m (h_out − h_in)/(1 − ξ)`.
"""
function compressor_power(c::Compressor, p)
    s = _compression_state(c, p)
    return p.m * (s.h_out - s.h_in) / (1 - c.ξ)
end
