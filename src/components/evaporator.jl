#=
Evaporator and borehole as components (feature/component-core).
=#

"""
    default_formulation(refrigerant)

`SaturationTemperature()` for a pure `ClapeyronRefrigerant` (explicit saturation from the
superancillaries when available, Clapeyron's saturation solver otherwise), `NestedPressure()`
for a mixture and for the spline `Refrigerant`. Determined by the refrigerant type only, so that
an `Evaporator` built with the default is inferred.
"""
default_formulation(::AbstractRefrigerant) = NestedPressure()
default_formulation(r::ClapeyronRefrigerant) = ispure(r) ? SaturationTemperature() : NestedPressure()

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
    default_Te_min(refrigerant)

Default lower bound of the evaporation temperature of an [`Evaporator`](@ref): the dew
temperature at a reduced pressure of 0.01, or just above the lowest pressure of the backend if
that is higher (R290: −60 °C; R744: −56.5 °C, its triple point). It keeps the solver away from
the zero-flow root at the triple point, where the vapour density, hence the refrigerant mass
flow and all heat flows, vanish.
"""
default_Te_min(r::AbstractRefrigerant) = T_dew(r, max(0.01 * critical_pressure(r), 1.001 * first(pressure_bounds(r))))

"""
    Evaporator(; refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup,
               formulation = default_formulation(refrigerant), Qscale = nominal_duty(...),
               Te_min = default_Te_min(refrigerant),
               refrigerant_pressure_drop = default_refrigerant_pressure_drop(refrigerant))

Two-zone (two-phase + superheat) moving-boundary plate evaporator, ε-NTU in each zone.

Two-phase zone: `Qe = ε C_s (T_s,m − T_r,in)` with the effectiveness of Qiao (2025,
[`effectiveness_counterflow`](@ref)): the refrigerant pressure drop over the zone
(`refrigerant_pressure_drop`, a [`TwoPhaseFriction`](@ref) model, by default [`Amalfi2016`](@ref)
for pure fluids; see [`twophase_pressure_drop`](@ref)) makes it enter at
`T_r,in = T_dew(pe + Δp)`, warmer than `Te = T_dew(pe)` at the end of the zone.
`refrigerant_pressure_drop = nothing` gives `T_r,in = Te`, `ε = 1 − e^{−NTU}`: the equations of
[`evaporator_model!`](@ref). The set of unknowns is the same in both cases. The inlet quality is
taken at `pe`; the superheat zone has no pressure drop inside the model (see
[`evaporator_pressure_drop`](@ref) for the compressor suction pressure).

- Ports: `Ts_in` (secondary-fluid inlet, the warm side), `Ts_out` (outlet).
- Unknowns: `pe` (`NestedPressure`, `LogBox` from the dew pressure at `Te_min` to the upper
  pressure bound) or `Te` (otherwise, `Box(Te_min, Tcrit)`), `Qe`, `Qsup` (`Linear(0, Qscale)`),
  `ϵsup`, `xs` (`Unit()`), `Ts_m` (secondary fluid at the moving boundary), and for
  `ExplicitTv` the molar volumes `vsup`, `vm`. (The two-phase effectiveness is explicit, not an
  unknown.)
- Residuals: heat balances scaled by `Qscale`, the superheat ε-NTU relation and the volume
  residuals as they are (dimensionless).

Guess: `guess(evap; Te, Qe, Qsup, Ts_m, ϵsup = 0.5, xs = 0.5)`.
"""
struct Evaporator{R <: AbstractRefrigerant, F <: PropertyFormulation, P <: PlateHeatExchanger, T <: Real,
                  X <: VariableTransform, S <: SecondaryFluid, D <: Union{Nothing, TwoPhaseFriction}} <: AbstractComponent
    refrigerant::R
    secondaryfluid::S
    plateheatexchanger::P
    hei::T
    Vs::T
    Vref::T
    ΔTsup::T
    formulation::F
    Qscale::T
    Te_min::T
    x1transform::X          # map of the first unknown (pe or Te), from Te_min
    pressure_drop::D        # refrigerant two-phase pressure drop model (nothing: none)
end

"""
    default_refrigerant_pressure_drop(refrigerant)

[`Amalfi2016`](@ref) for pure fluids; `nothing` (no refrigerant pressure drop) for mixtures, whose
two-phase inlet temperature also needs the glide (not implemented yet). Determined by the
refrigerant type, so that an `Evaporator` built with the default is inferred.
"""
default_refrigerant_pressure_drop(::AbstractRefrigerant) = Amalfi2016()
default_refrigerant_pressure_drop(r::ClapeyronRefrigerant) = ispure(r) ? Amalfi2016() : nothing

function Evaporator(; refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup,
                    formulation::PropertyFormulation = default_formulation(refrigerant),
                    Qscale = nominal_duty(refrigerant, Vref, hei), Te_min = default_Te_min(refrigerant),
                    refrigerant_pressure_drop::Union{Nothing, TwoPhaseFriction} = default_refrigerant_pressure_drop(refrigerant))
    refrigerant_pressure_drop === nothing || refrigerant isa ClapeyronRefrigerant && !ispure(refrigerant) &&
        throw(ArgumentError("refrigerant pressure drop: pure fluids only so far (mixtures need the glide)"))
    hei, Vs, Vref, ΔTsup, Qscale, Te_min = promote(float.((hei, Vs, Vref, ΔTsup, Qscale, Te_min))...)
    x1 = _first_transform(formulation, refrigerant, Te_min)
    return Evaporator(refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, formulation, Qscale,
                      Te_min, x1, refrigerant_pressure_drop)
end

function _first_transform(::NestedPressure, r, Te_min)
    pmin, pmax = pressure_bounds(r)
    return LogBox(max(pmin, _p_dew(r, Te_min)), pmax)
end
function _first_transform(F::PropertyFormulation, r, Te_min)
    lo, hi = bounds(first_transform(F, r))
    return Box(max(lo, Te_min), hi)
end

_first_name(::NestedPressure) = :pe
_first_name(::PropertyFormulation) = :Te
_extra_names(::PropertyFormulation) = ()
_extra_names(::ExplicitTv) = (:vsup, :vm)

function unknowns(c::Evaporator)
    F, r, Q = c.formulation, c.refrigerant, Linear(0, c.Qscale)
    names = (_first_name(F), :Qe, :Qsup, :ϵsup, :xs, :Ts_m, _extra_names(F)...)
    return NamedTuple{names}((c.x1transform, Q, Q, Unit(), Unit(), Free(), extra_transforms(F, r)...))
end

ports(::Evaporator) = (Ts_in = Free(), Ts_out = Free())

function residual_scales(c::Evaporator)
    Q = c.Qscale
    return (Q, Q, Q, Q, Q, Q, 1.0, map(_ -> 1.0, _extra_names(c.formulation))...)
end

function residuals!(r, c::Evaporator, x, p)
    F = c.formulation
    vols = map(k -> x[k], _extra_names(F))
    rs = refrigerant_state!(view(r, 8:7+length(vols)), F, c.refrigerant, x[1], vols, c.ΔTsup)
    _evaporator_balances!(r, c, x, p.Ts_in, p.Ts_out, rs)
    return nothing
end

# Refrigerant temperature at the inlet of the two-phase zone and the pressure drop over it, for
# the refrigerant mass flow m and the zone fraction xs (rs: refrigerant state at the end of the zone)
_twophase_inlet(::Nothing, c::Evaporator, rs, m, xs) = (rs.Te, zero(rs.Te))
function _twophase_inlet(model::TwoPhaseFriction, c::Evaporator, rs, m, xs)
    r, pe, phx = c.refrigerant, rs.pe, c.plateheatexchanger
    hl = h_bubble(r, pe)
    x_in = (c.hei - hl) / (rs.hes - hl)
    Δp = twophase_pressure_drop(model, r, phx, m, pe, x_in, one(x_in); L = xs * phx.L).total
    return (T_dew(r, pe + Δp), Δp)
end

# two-phase zone: capacity of the secondary fluid, NTU, effectiveness and refrigerant inlet temperature
function _twophase_zone(c::Evaporator, Qe, xs, Ts_l, Ts_m, rs)
    sf, phx = c.secondaryfluid, c.plateheatexchanger
    Tsm = (Ts_l + Ts_m)/2
    U, _, _ = evaporator_uvalue_twophase(Tsm, c.Vs, Qe, rs.pe, xs, sf, c.refrigerant, phx)
    Cs = sf.ρ(Tsm) * sf.cp(Tsm) * c.Vs
    Tr_in, Δp = _twophase_inlet(c.pressure_drop, c, rs, rs.ρsup * c.Vref, xs)
    NTU = U * phx.A * xs / Cs
    γ = (rs.Te - Tr_in) / (Tr_in - Ts_m)
    return (U = U, Cs = Cs, NTU = NTU, γ = γ, ϵ = effectiveness_counterflow(NTU, γ), Tr_in = Tr_in, Δp = Δp)
end

# superheat zone: U value, capacities of the secondary fluid and of the refrigerant (mean over the
# superheat), NTU referred to the smaller capacity, capacity ratio
function _superheat_zone(c::Evaporator, xs, Ts_h, Ts_m, rs)
    sf, phx = c.secondaryfluid, c.plateheatexchanger
    Tsm = (Ts_h + Ts_m)/2
    U, _, _ = evaporator_uvalue_superheat(Tsm, c.Vs, c.Vref, rs, sf, phx)
    Cs = sf.ρ(Tsm) * sf.cp(Tsm) * c.Vs
    Cref = rs.ρsup * c.Vref * (rs.hsup - rs.hes) / c.ΔTsup
    Cmin = min(Cref, Cs)
    return (U = U, Cs = Cs, Cref = Cref, Cmin = Cmin, Cr = Cmin / max(Cref, Cs), NTU = U * phx.A * (1 - xs) / Cmin)
end

function _evaporator_balances!(du, c::Evaporator, x, Ts_h, Ts_l, rs)
    (; Qe, Qsup, ϵsup, xs, Ts_m) = x
    ΔTsup = c.ΔTsup
    (; Te, hes, ρsup) = rs
    m = ρsup * c.Vref

    # two-phase (counter-flow, refrigerant temperature changing with the pressure drop)
    z = _twophase_zone(c, Qe, xs, Ts_l, Ts_m, rs)

    # superheat (as evaporator_balances!)
    h = _superheat_zone(c, xs, Ts_h, Ts_m, rs)

    du[1] = z.Cs * (Ts_m - Ts_l) - Qe                 # secondary fluid, two-phase zone
    du[2] = m * (hes - c.hei) - Qe                     # refrigerant, two-phase zone
    du[3] = z.ϵ * z.Cs * (Ts_m - z.Tr_in) - Qe         # ε-NTU (Qiao)
    du[4] = h.Cs * (Ts_h - Ts_m) - Qsup                # secondary fluid, superheat
    du[5] = h.Cref * ΔTsup - Qsup                      # refrigerant, superheat
    du[6] = ϵsup * h.Cmin * (Ts_h - Te) - Qsup         # ε-NTU
    du[7] = (1 - exp(-h.NTU * (1 - h.Cr))) / (1 - h.Cr * exp(-h.NTU * (1 - h.Cr))) - ϵsup
    return nothing
end

function guess(c::Evaporator; Te, Qe, Qsup, Ts_m, ϵsup = 0.5, xs = 0.5)
    F, r = c.formulation, c.refrigerant
    x1 = F isa NestedPressure ? _p_dew(r, Te) : Te
    extra = if F isa ExplicitTv
        pe = _p_dew(r, Te)
        (volume_pT(r, pe, Te + c.ΔTsup; phase = :vapor), volume_pT(r, pe, Te + c.ΔTsup/2; phase = :vapor))
    else
        ()
    end
    return NamedTuple{keys(unknowns(c))}((x1, Qe, Qsup, ϵsup, xs, Ts_m, extra...))
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
    pe = _pe(c, x)
    # the original solution types: two-phase effectiveness in the classical definition, i.e.
    # referred to the refrigerant temperature Te at the end of the zone
    Tsm = (p.Ts_out + x.Ts_m)/2
    ϵe = x.Qe / (c.secondaryfluid.ρ(Tsm) * c.secondaryfluid.cp(Tsm) * c.Vs * (x.Ts_m - T_dew(r, pe)))
    u = [to_R(pressure_transform(r), pe), x.Qe, x.Qsup, to_R(Unit(), ϵe), to_R(Unit(), x.ϵsup),
         to_R(Unit(), x.xs), p.Ts_out, x.Ts_m]
    return evaporator_model(u, EvaporatorParameters(c, p.Ts_in))
end

_pe(c::Evaporator, x) = c.formulation isa NestedPressure ? x[1] : p_sat(c.refrigerant, x[1])

"""
    evaporator_zones(c::Evaporator, x, p) -> (; twophase, superheat)

Zone quantities of a solution (`x` physical unknowns, `p` port values), as used in the residuals:
- `twophase`: U value `U`, brine capacity `Cs`, `NTU = U A xs / Cs`, `γ`, refrigerant temperatures
  `Tr_in` (inlet) and `Te` (end of the zone), pressure drop `Δp`;
- `superheat`: `U`, brine and refrigerant capacities `Cs`, `Cref` (mean over the superheat),
  `NTU_s = U A (1 − xs)/Cs` and `NTU_r = U A (1 − xs)/Cref`.
Constant within each zone, so the temperature profiles follow in closed form (e.g.
`examples/evaporator_profiles.jl`). The refrigerant state is evaluated from the pressure
(`NestedPressure`) whatever the formulation.
"""
function evaporator_zones(c::Evaporator, x, p)
    pe = _pe(c, x)
    rs = refrigerant_state!(nothing, NestedPressure(), c.refrigerant, pe, (), c.ΔTsup)
    z = _twophase_zone(c, x.Qe, x.xs, p.Ts_out, x.Ts_m, rs)
    h = _superheat_zone(c, x.xs, p.Ts_in, x.Ts_m, rs)
    Asup = c.plateheatexchanger.A * (1 - x.xs)
    return (twophase = (U = z.U, Cs = z.Cs, NTU = z.NTU, γ = z.γ, Tr_in = z.Tr_in, Te = rs.Te, Δp = z.Δp),
            superheat = (U = h.U, Cs = h.Cs, Cref = h.Cref, NTU_s = h.U * Asup / h.Cs, NTU_r = h.U * Asup / h.Cref))
end

"""
    evaporator_pressure_drop(c::Evaporator, x, p) -> NamedTuple

Refrigerant side of a solution (`x` physical unknowns, `p` port values): two-phase pressure drop
`Δp_twophase` [Pa], refrigerant temperature at the inlet of the two-phase zone `Tr_in` and at its
end `Te`, `γ` and `NTU` of the zone, and the single-phase pressure drop of the superheat zone
`Δp_superheat` (Martin; for the compressor suction pressure, `compressor_power(…; Δp_suction)`).
"""
function evaporator_pressure_drop(c::Evaporator, x, p)
    r, phx = c.refrigerant, c.plateheatexchanger
    pe = _pe(c, x)
    Te = T_dew(r, pe)
    Tsup = Te + c.ΔTsup
    ρsup = ρ_pT(r, pe, Tsup; phase = :vapor)
    rs = (pe = pe, Te = Te, hes = h_dew(r, pe), ρsup = ρsup)
    z = _twophase_zone(c, x.Qe, x.xs, p.Ts_out, x.Ts_m, rs)
    m = ρsup * c.Vref
    Tm = Te + c.ΔTsup/2
    Δp_sup = vapour_pressure_drop(phx, m, ρ_pT(r, pe, Tm; phase = :vapor), μ_pT(r, pe, Tm; phase = :vapor), (1 - x.xs) * phx.L)
    return (Δp_twophase = z.Δp, Tr_in = z.Tr_in, Te = Te, γ = z.γ, NTU = z.NTU, Δp_superheat = Δp_sup)
end

"""
    BoreholeHX(; borehole, secondaryfluid, Vs, λground)

Borehole heat exchanger (Hellström, uniform borehole wall temperature `Tb`), as in
[`borehole_model!`](@ref).

- Ports: `Tin`, `Tout` (secondary fluid) and `Tb`, the ground (borehole wall) temperature:
  usually a fixed node, e.g. `nodes = (…, T_ground = (:borehole => :Tb,)), fixed = (T_ground = 281.15,)` (8 °C).
  Change it with [`set_fixed`](@ref) (same system type, nothing recompiled), or connect it to a
  component that computes it (e.g. a dynamic ground model).
- No internal unknowns; one residual `kin Tin - Tout + (1 - kin) Tb` [K], with the coefficient
  and the fluid density at the mean temperature `(Tin + Tout)/2`.
"""
struct BoreholeHX{B <: Borehole, T <: Real, S <: SecondaryFluid} <: AbstractComponent
    borehole::B
    secondaryfluid::S
    Vs::T
    λground::T
end
function BoreholeHX(; borehole, secondaryfluid, Vs, λground, Tb = nothing)
    Tb === nothing || throw(ArgumentError("BoreholeHX: `Tb` is a port; connect it to a node, e.g. " *
        "`nodes = (…, T_ground = (:borehole => :Tb,)), fixed = (T_ground = $Tb,)`"))
    Vs, λground = promote(float(Vs), float(λground))
    return BoreholeHX(borehole, secondaryfluid, Vs, λground)
end

ports(::BoreholeHX) = (Tin = Free(), Tout = Free(), Tb = Free())
residual_scales(::BoreholeHX) = (1.0,)

function residuals!(r, c::BoreholeHX, x, p)
    Tm = (p.Tin + p.Tout)/2
    sf = c.secondaryfluid
    kin = uniform_Tb_coeffs(c.borehole, c.λground, sf.ρ(Tm) * c.Vs, Tm, sf)
    r[1] = kin * p.Tin - p.Tout + (1 - kin) * p.Tb
    return nothing
end

"Post-processing: `BoreholeSolution` at the port values `p`."
borehole_solution(c::BoreholeHX, p) = BoreholeSolution(c.borehole, c.λground, p.Tin, p.Tout, p.Tb, c.Vs, c.secondaryfluid)

"""
    validity(c::Evaporator, x, p; ΔTpinch = 0.1, xs_tol = 1e-3)

Issues reported for an evaporator solution:
- secondary fluid (`Ts_in`, `Ts_m`, `Ts_out`) below its freezing point or outside its property
  table;
- evaporation temperature at its lower bound `Te_min` (within 0.01 K): not a physical solution;
- pinch (within `ΔTpinch`): `Ts_out - Te` (two-phase zone) or `Ts_in - (Te + ΔTsup)` (superheat);
- a zone almost vanished: `xs < xs_tol` or `xs > 1 - xs_tol`.
"""
function validity(c::Evaporator, x, p; ΔTpinch = 0.1, xs_tol = 1e-3)
    issues = NamedTuple[]
    _secondary_issues!(issues, c.secondaryfluid, (Ts_in = p.Ts_in, Ts_m = x.Ts_m, Ts_out = p.Ts_out))
    Te = c.formulation isa NestedPressure ? T_dew(c.refrigerant, x.pe) : x.Te
    Te < c.Te_min + 0.01 &&
        push!(issues, _issue(:Te, Te, c.Te_min, "evaporation temperature at its lower bound Te_min"))
    p.Ts_out - Te < ΔTpinch &&
        push!(issues, _issue(:Te, Te, p.Ts_out - ΔTpinch, "two-phase pinch: Ts_out - Te below ΔTpinch"))
    p.Ts_in - (Te + c.ΔTsup) < ΔTpinch &&
        push!(issues, _issue(:Te, Te, p.Ts_in - c.ΔTsup - ΔTpinch, "superheat pinch: Ts_in - (Te + ΔTsup) below ΔTpinch"))
    x.xs < xs_tol && push!(issues, _issue(:xs, x.xs, xs_tol, "two-phase zone almost vanished"))
    x.xs > 1 - xs_tol && push!(issues, _issue(:xs, x.xs, 1 - xs_tol, "superheat zone almost vanished"))
    return issues
end

"Issues reported for a borehole: secondary fluid (`Tin`, `Tout`) below freezing or outside its table."
validity(c::BoreholeHX, x, p) = _secondary_issues!(NamedTuple[], c.secondaryfluid, (Tin = p.Tin, Tout = p.Tout))
