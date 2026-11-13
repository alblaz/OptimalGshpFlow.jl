#=
Evaporator and borehole as components (feature/component-core).
=#

"""
    default_formulation(refrigerant)

`SaturationTemperature()` for a pure `ClapeyronRefrigerant` with explicit saturation
(superancillaries), `NestedPressure()` otherwise.
"""
default_formulation(::AbstractRefrigerant) = NestedPressure()
default_formulation(r::ClapeyronRefrigerant) = ispure(r) && r.superancillary ? SaturationTemperature() : NestedPressure()

"""
    nominal_duty(refrigerant, Vref, hei; Te = 273.15)

Evaporator duty [W] used to scale heat flows and heat balances: the two-phase duty
`ρ_dew Vref (h_dew - hei)` at the dew temperature `Te`.
"""
function nominal_duty(r, Vref, hei; Te = 273.15)
    p = _p_dew(r, Te)
    return abs(ρ_dew(r, p) * Vref * (h_dew(r, p) - hei))
end

"""
    Evaporator(; refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup,
               formulation = default_formulation(refrigerant), Qscale = nominal_duty(...))

Two-zone (two-phase + superheat) moving-boundary plate evaporator, ε-NTU in each zone; the same
equations as [`evaporator_model!`](@ref).

- Ports: `Ts_in` (secondary-fluid inlet, the warm side), `Ts_out` (outlet).
- Unknowns: `pe` (`NestedPressure`) or `Te` (otherwise), `Qe`, `Qsup` (`Linear(0, Qscale)`),
  `ϵe`, `ϵsup`, `xs` (`Unit()`), `Ts_m` (secondary fluid at the moving boundary), and for
  `ExplicitTv` the molar volumes `vsup`, `vm`.
- Residuals: heat balances scaled by `Qscale`, the ε-NTU relations and the volume residuals as
  they are (dimensionless).

Guess: `guess(evap; Te, Qe, Qsup, Ts_m, ϵe = 0.5, ϵsup = 0.5, xs = 0.5)`.
"""
struct Evaporator{R <: AbstractRefrigerant, F <: PropertyFormulation, P <: PlateHeatExchanger, T <: Real} <: AbstractComponent
    refrigerant::R
    secondaryfluid::SecondaryFluid
    plateheatexchanger::P
    hei::T
    Vs::T
    Vref::T
    ΔTsup::T
    formulation::F
    Qscale::T
end

function Evaporator(; refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup,
                    formulation::PropertyFormulation = default_formulation(refrigerant),
                    Qscale = nominal_duty(refrigerant, Vref, hei))
    hei, Vs, Vref, ΔTsup, Qscale = promote(float.((hei, Vs, Vref, ΔTsup, Qscale))...)
    return Evaporator(refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, formulation, Qscale)
end

_first_name(::NestedPressure) = :pe
_first_name(::PropertyFormulation) = :Te
_extra_names(::PropertyFormulation) = ()
_extra_names(::ExplicitTv) = (:vsup, :vm)

function unknowns(c::Evaporator)
    F, r, Q = c.formulation, c.refrigerant, Linear(0, c.Qscale)
    names = (_first_name(F), :Qe, :Qsup, :ϵe, :ϵsup, :xs, :Ts_m, _extra_names(F)...)
    return NamedTuple{names}((first_transform(F, r), Q, Q, Unit(), Unit(), Unit(), Free(), extra_transforms(F, r)...))
end

ports(::Evaporator) = (Ts_in = Free(), Ts_out = Free())

function residual_scales(c::Evaporator)
    Q = c.Qscale
    return (Q, Q, Q, 1.0, Q, Q, Q, 1.0, map(_ -> 1.0, _extra_names(c.formulation))...)
end

function residuals!(r, c::Evaporator, x, p)
    F = c.formulation
    vols = map(k -> x[k], _extra_names(F))
    rs = refrigerant_state!(view(r, 9:8+length(vols)), F, c.refrigerant, x[1], vols, c.ΔTsup)
    xb = (Qe = x.Qe, Qsup = x.Qsup, ϵe = x.ϵe, ϵsup = x.ϵsup, xs = x.xs, Ts_l = p.Ts_out, Ts_m = x.Ts_m)
    evaporator_balances!(r, xb, p.Ts_in, rs, c.refrigerant, c.secondaryfluid, c.plateheatexchanger,
                         c.hei, c.Vs, c.Vref, c.ΔTsup)
    return nothing
end

function guess(c::Evaporator; Te, Qe, Qsup, Ts_m, ϵe = 0.5, ϵsup = 0.5, xs = 0.5)
    F, r = c.formulation, c.refrigerant
    x1 = F isa NestedPressure ? _p_dew(r, Te) : Te
    extra = if F isa ExplicitTv
        pe = _p_dew(r, Te)
        (volume_pT(r, pe, Te + c.ΔTsup; phase = :vapor), volume_pT(r, pe, Te + c.ΔTsup/2; phase = :vapor))
    else
        ()
    end
    return NamedTuple{keys(unknowns(c))}((x1, Qe, Qsup, ϵe, ϵsup, xs, Ts_m, extra...))
end

"""
    EvaporatorParameters(c::Evaporator, Ts_h)

Parameters of the original model ([`evaporator_model!`](@ref)) for the same evaporator, with
secondary-fluid inlet temperature `Ts_h`.
"""
EvaporatorParameters(c::Evaporator, Ts_h) =
    EvaporatorParameters(refrigerant = c.refrigerant, secondaryfluid = c.secondaryfluid,
                         plateheatexchanger = c.plateheatexchanger, hei = c.hei, Vs = c.Vs, Vref = c.Vref,
                         ΔTsup = c.ΔTsup, Ts_h = Ts_h)

"""
    evaporator_solution(c::Evaporator, x, p) -> EvaporatorSolution

Post-processing of a solution: `x` physical unknowns and `p` port values of the evaporator (e.g.
from [`physical`](@ref) and [`port_values`](@ref)).
"""
function evaporator_solution(c::Evaporator, x, p)
    r = c.refrigerant
    pe = c.formulation isa NestedPressure ? x[1] : p_sat(r, x[1])
    u = [to_R(pressure_transform(r), pe), x.Qe, x.Qsup, to_R(Unit(), x.ϵe), to_R(Unit(), x.ϵsup),
         to_R(Unit(), x.xs), p.Ts_out, x.Ts_m]
    return evaporator_model(u, EvaporatorParameters(c, p.Ts_in))
end

"""
    BoreholeHX(; borehole, secondaryfluid, Vs, λground, Tb)

Borehole heat exchanger (Hellström, uniform borehole wall temperature `Tb`), as in
[`borehole_model!`](@ref).

- Ports: `Tin`, `Tout` (secondary fluid).
- No internal unknowns; one residual `kin Tin - Tout + (1 - kin) Tb` [K], with the coefficient
  and the fluid density at the mean temperature `(Tin + Tout)/2`.
"""
struct BoreholeHX{B <: Borehole, T <: Real} <: AbstractComponent
    borehole::B
    secondaryfluid::SecondaryFluid
    Vs::T
    λground::T
    Tb::T
end
function BoreholeHX(; borehole, secondaryfluid, Vs, λground, Tb)
    Vs, λground, Tb = promote(float.((Vs, λground, Tb))...)
    return BoreholeHX(borehole, secondaryfluid, Vs, λground, Tb)
end

ports(::BoreholeHX) = (Tin = Free(), Tout = Free())
residual_scales(::BoreholeHX) = (1.0,)

function residuals!(r, c::BoreholeHX, x, p)
    Tm = (p.Tin + p.Tout)/2
    sf = c.secondaryfluid
    kin = uniform_Tb_coeffs(c.borehole, c.λground, sf.ρ(Tm) * c.Vs, Tm, sf)
    r[1] = kin * p.Tin - p.Tout + (1 - kin) * c.Tb
    return nothing
end

"Post-processing: `BoreholeSolution` at the port values `p`."
borehole_solution(c::BoreholeHX, p) = BoreholeSolution(c.borehole, c.λground, p.Tin, p.Tout, c.Tb, c.Vs, c.secondaryfluid)
