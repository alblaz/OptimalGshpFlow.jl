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
(`refrigerant_pressure_drop`, a [`TwoPhaseFriction`](@ref) model, by default [`Amalfi2016`](@ref);
see [`twophase_pressure_drop`](@ref)) makes it enter at
`T_r,in = T_dew(pe + Δp)`, warmer than `Te = T_dew(pe)` at the end of the zone.
`refrigerant_pressure_drop = nothing` gives `T_r,in = Te`, `ε = 1 − e^{−NTU}`: the equations of
[`evaporator_model!`](@ref). The set of unknowns is the same in both cases.
Zeotropic blend (`ClapeyronRefrigerant` mixture): the refrigerant enters at `T_r,in = T_ph(p_in, hei)`
(flash at the inlet enthalpy and pressure `p_in = pe + Δp`), below the dew temperature, and its
temperature changes linearly along the zone: the glide `T_dew(p_in) − T_r,in` with the pseudo
capacity `C_r = ṁ (h_dew − hei)/glide` (`φ = C_s/C_r`), and the saturation shift
`Te − T_dew(p_in)` of the pressure drop (`γ`). The surface tension of a blend is a mole-fraction
average of its components ([`surface_tension`](@ref)). The boiling heat-transfer coefficient
(Cooper with the mixture properties) is reduced by the mass-transfer resistance of the blend
([`mixture_boiling_factor`](@ref), Thome & Shakir 1987). The inlet quality is
taken at `pe`; the superheat zone has no pressure drop inside the model (see
[`evaporator_pressure_drop`](@ref) for the compressor suction pressure).

- Ports: `Ts_in` (secondary-fluid inlet, the warm side), `Ts_out` (outlet). In a closed cycle
  (`hei = nothing`, give `hei_nominal` for the scales and the guess) also `hei` (refrigerant inlet
  enthalpy, from the condenser outlet through the expansion valve) and the outlet ports `p_suc`
  (`pe`), `T_suc` (`Te + ΔTsup`) and `m` (`ρ_sup Vref`), each with its residual; transforms from
  [`CyclePorts`](@ref).
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
                  X <: VariableTransform, S <: SecondaryFluid, D <: Union{Nothing, TwoPhaseFriction},
                  H <: Union{Nothing, T}, K <: Union{Nothing, CyclePorts{T}}} <: AbstractComponent
    refrigerant::R
    secondaryfluid::S
    plateheatexchanger::P
    hei::H                  # refrigerant inlet enthalpy; nothing: a port (closed cycle)
    Vs::T
    Vref::T
    ΔTsup::T
    formulation::F
    Qscale::T
    Te_min::T
    x1transform::X          # map of the first unknown (pe or Te), from Te_min
    pressure_drop::D        # refrigerant two-phase pressure drop model (nothing: none)
    cycleports::K           # transforms of the cycle ports (hei = nothing), else nothing
end


"""
    default_refrigerant_pressure_drop(refrigerant)

[`Amalfi2016`](@ref) for every refrigerant (pure fluids and blends; `nothing` in `Evaporator`
for no refrigerant pressure drop).
"""
default_refrigerant_pressure_drop(::AbstractRefrigerant) = Amalfi2016()

function Evaporator(; refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup,
                    formulation::PropertyFormulation = default_formulation(refrigerant),
                    hei_nominal = hei, Qscale = nothing, Te_min = default_Te_min(refrigerant),
                    refrigerant_pressure_drop::Union{Nothing, TwoPhaseFriction} = default_refrigerant_pressure_drop(refrigerant))
    hei_nominal === nothing && throw(ArgumentError("Evaporator: with `hei = nothing` (closed cycle) give `hei_nominal`"))
    Qscale = something(Qscale, nominal_duty(refrigerant, Vref, hei_nominal))
    Vs, Vref, ΔTsup, Qscale, Te_min = promote(float.((Vs, Vref, ΔTsup, Qscale, Te_min))...)
    x1 = _first_transform(formulation, refrigerant, Te_min)
    hei, cp = hei === nothing ? (nothing, CyclePorts(refrigerant)) : (convert(typeof(Vs), hei), nothing)
    return Evaporator(refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, formulation, Qscale,
                      Te_min, x1, refrigerant_pressure_drop, cp)
end

# inlet enthalpy: the parameter, or the port of a cycle evaporator
_hei(hei::Real, p) = hei
_hei(::Nothing, p) = p.hei
_hei(c::Evaporator, p) = _hei(c.hei, p)

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

# a closed-cycle evaporator (cycleports, hei = nothing): inlet enthalpy and outlet state as ports
# (dispatch on the field: a method on Evaporator{…, Nothing} is not taken as more specific)
ports(c::Evaporator) = _evaporator_ports(c.cycleports)
_evaporator_ports(::Nothing) = (Ts_in = Free(), Ts_out = Free())
_evaporator_ports(cp::CyclePorts) = (Ts_in = Free(), Ts_out = Free(), hei = cp.h, p_suc = cp.p, T_suc = Free(), m = cp.m)

function residual_scales(c::Evaporator)
    Q = c.Qscale
    return (Q, Q, Q, Q, Q, Q, 1.0, map(_ -> 1.0, _extra_names(c.formulation))..., _cycle_scales(c)...)
end
_cycle_scales(c::Evaporator) = _cycle_scales(c.cycleports, c)
_cycle_scales(::Nothing, c) = ()
_cycle_scales(cp::CyclePorts, c) = (cp.pscale, 1.0, c.Qscale / cp.hscale)

function residuals!(r, c::Evaporator, x, p)
    F = c.formulation
    vols = map(k -> x[k], _extra_names(F))
    rs = refrigerant_state!(view(r, 8:7+length(vols)), F, c.refrigerant, x[1], vols, c.ΔTsup)
    _evaporator_balances!(r, c, x, p.Ts_in, p.Ts_out, rs, _hei(c, p))
    _cycle_residuals!(r, 7 + length(vols), c.cycleports, c, p, rs)
    return nothing
end

_cycle_residuals!(r, n, ::Nothing, c, p, rs) = nothing
function _cycle_residuals!(r, n, ::CyclePorts, c, p, rs)
    r[n+1] = p.p_suc - rs.pe                           # outlet pressure (end of the two-phase zone)
    r[n+2] = p.T_suc - (rs.Te + c.ΔTsup)               # outlet temperature
    r[n+3] = p.m - rs.ρsup * c.Vref                    # mass flow (displacement at the outlet density)
    return nothing
end

# Refrigerant temperature at the inlet of the two-phase zone and the pressure drop over it, for
# the refrigerant mass flow m and the zone fraction xs (rs: refrigerant state at the end of the zone)
_twophase_inlet(::Nothing, c::Evaporator, rs, m, xs, hei) = (rs.Te, zero(rs.Te))
function _twophase_inlet(model::TwoPhaseFriction, c::Evaporator, rs, m, xs, hei)
    r, pe, phx = c.refrigerant, rs.pe, c.plateheatexchanger
    hl = h_bubble(r, pe)
    x_in = (hei - hl) / (rs.hes - hl)
    Δp = twophase_pressure_drop(model, r, phx, m, pe, x_in, one(x_in); L = xs * phx.L).total
    return (T_dew(r, pe + Δp), Δp)
end

# refrigerant temperature at the inlet of the two-phase zone at the pressure pe, from the inlet
# enthalpy (flash): below Te = T_dew(pe) for a zeotropic blend (glide), Te for a pure fluid
_glide_inlet_T(::AbstractRefrigerant, pe, hei, Te) = Te
_glide_inlet_T(r::ClapeyronRefrigerant, pe, hei, Te) = ispure(r) ? Te : T_ph(r, pe, hei)

# two-phase zone: capacity of the secondary fluid, NTU, effectiveness and refrigerant inlet
# temperature. Refrigerant temperature linear along the zone (Qiao 2025): it rises by the glide
# Te − T_in (blends; capacity C_r = ṁ (h_dew − hei)/(Te − T_in), φ = Cs/C_r) and falls by the
# saturation shift of the pressure drop, Te − T_dew(pe + Δp) < 0 (γ).
function _twophase_zone(c::Evaporator, Qe, xs, Ts_l, Ts_m, rs, hei)
    sf, phx = c.secondaryfluid, c.plateheatexchanger
    Tsm = (Ts_l + Ts_m)/2
    U, _, _ = evaporator_uvalue_twophase(Tsm, c.Vs, Qe, rs.pe, xs, sf, c.refrigerant, phx)
    Cs = sf.ρ(Tsm) * sf.cp(Tsm) * c.Vs
    m = rs.ρsup * c.Vref
    # pressure drop: dew temperature at the inlet pressure p_in = pe + Δp (Tsat_in); blends enter
    # below it by the glide, T_r,in = T(p_in, hei). Te − T_r,in = (Te − Tsat_in) + glide: the
    # saturation shift goes into γ, the glide into φ.
    Tsat_in, Δp = _twophase_inlet(c.pressure_drop, c, rs, m, xs, hei)
    Tr_in = _glide_inlet_T(c.refrigerant, rs.pe + Δp, hei, Tsat_in)
    glide = Tsat_in - Tr_in
    φ = Cs * glide / (m * (rs.hes - hei))
    NTU = U * phx.A * xs / Cs
    γ = (rs.Te - Tsat_in) / (Tr_in - Ts_m)
    return (U = U, Cs = Cs, NTU = NTU, γ = γ, φ = φ, ϵ = effectiveness_counterflow(NTU, γ, φ), Tr_in = Tr_in, Δp = Δp)
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

function _evaporator_balances!(du, c::Evaporator, x, Ts_h, Ts_l, rs, hei)
    (; Qe, Qsup, ϵsup, xs, Ts_m) = x
    ΔTsup = c.ΔTsup
    (; Te, hes, ρsup) = rs
    m = ρsup * c.Vref

    # two-phase (counter-flow, refrigerant temperature changing with the pressure drop)
    z = _twophase_zone(c, Qe, xs, Ts_l, Ts_m, rs, hei)

    # superheat (as evaporator_balances!)
    h = _superheat_zone(c, xs, Ts_h, Ts_m, rs)

    du[1] = z.Cs * (Ts_m - Ts_l) - Qe                 # secondary fluid, two-phase zone
    du[2] = m * (hes - hei) - Qe                     # refrigerant, two-phase zone
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
    EvaporatorParameters(c::Evaporator, Ts_h, hei = c.hei)

Parameters of the original model ([`evaporator_model!`](@ref)) for the same evaporator, with
secondary-fluid inlet temperature `Ts_h` (and inlet enthalpy `hei`: needed for a cycle evaporator).
"""
EvaporatorParameters(c::Evaporator, Ts_h, hei = c.hei) =
    EvaporatorParameters(refrigerant = c.refrigerant, secondaryfluid = c.secondaryfluid,
                         plateheatexchanger = c.plateheatexchanger, hei = hei, Vs = c.Vs, Vref = c.Vref,
                         ΔTsup = c.ΔTsup, Ts_h = Ts_h)

"""
    evaporator_solution(c::Evaporator, x, p) -> EvaporatorSolution

Post-processing of a solution: `x` physical unknowns and `p` port values of the evaporator (e.g.
from [`physical`](@ref) and [`port_values`](@ref)).
"""
function evaporator_solution(c::Evaporator, x, p)
    r = c.refrigerant
    pe = _pe(c, x)
    # the original solution types: two-phase effectiveness referred to a lower bound of the
    # refrigerant temperature (so ≤ 1): Te at the end of the zone, or T(pe, hei) of a blend
    Tsm = (p.Ts_out + x.Ts_m)/2
    Te = T_dew(r, pe)
    hei = _hei(c, p)
    ϵe = x.Qe / (c.secondaryfluid.ρ(Tsm) * c.secondaryfluid.cp(Tsm) * c.Vs * (x.Ts_m - _glide_inlet_T(r, pe, hei, Te)))
    u = [to_R(pressure_transform(r), pe), x.Qe, x.Qsup, to_R(Unit(), ϵe), to_R(Unit(), x.ϵsup),
         to_R(Unit(), x.xs), p.Ts_out, x.Ts_m]
    return evaporator_model(u, EvaporatorParameters(c, p.Ts_in, hei))
end

_pe(c::Evaporator, x) = c.formulation isa NestedPressure ? x[1] : p_sat(c.refrigerant, x[1])

"""
    evaporator_zones(c::Evaporator, x, p) -> (; twophase, superheat)

Zone quantities of a solution (`x` physical unknowns, `p` port values), as used in the residuals:
- `twophase`: U value `U`, brine capacity `Cs`, `NTU = U A xs / Cs`, `γ` (pressure drop), `φ`
  (glide: brine capacity over the refrigerant pseudo capacity, 0 for a pure fluid), refrigerant
  temperatures `Tr_in` (inlet) and `Te` (end of the zone), pressure drop `Δp`;
- `superheat`: `U`, brine and refrigerant capacities `Cs`, `Cref` (mean over the superheat),
  `NTU_s = U A (1 − xs)/Cs` and `NTU_r = U A (1 − xs)/Cref`.
Constant within each zone, so the temperature profiles follow in closed form (e.g.
`examples/evaporator_profiles.jl`). The refrigerant state is evaluated from the pressure
(`NestedPressure`) whatever the formulation.
"""
function evaporator_zones(c::Evaporator, x, p)
    pe = _pe(c, x)
    rs = refrigerant_state!(nothing, NestedPressure(), c.refrigerant, pe, (), c.ΔTsup)
    z = _twophase_zone(c, x.Qe, x.xs, p.Ts_out, x.Ts_m, rs, _hei(c, p))
    h = _superheat_zone(c, x.xs, p.Ts_in, x.Ts_m, rs)
    Asup = c.plateheatexchanger.A * (1 - x.xs)
    return (twophase = (U = z.U, Cs = z.Cs, NTU = z.NTU, γ = z.γ, φ = z.φ, Tr_in = z.Tr_in, Te = rs.Te, Δp = z.Δp),
            superheat = (U = h.U, Cs = h.Cs, Cref = h.Cref, NTU_s = h.U * Asup / h.Cs, NTU_r = h.U * Asup / h.Cref))
end

"""
    evaporator_profiles(c::Evaporator, x, p; n = 100) -> NamedTuple

Temperature profiles along the evaporator from a solution (`x` physical unknowns, `p` port values):
area fraction `a` from the refrigerant inlet (0) to the outlet (1), refrigerant `Tr` and secondary
fluid `Ts`, `n` points per zone, the zone boundary `xs`, and `end_error`: the largest difference [K]
between the far ends of the integrated profiles and the solution (`Ts_m`, `Te + ΔTsup`, `Ts_in`).
- Two-phase zone (0 … xs), as the model (Qiao 2025): the refrigerant temperature changes by the
  pressure drop linearly in the area (ΔT_p = Te − T_sat,in) and by the glide linearly in the heat
  (φ = C_s glide / Qe); with w = Ts − Ts_out, w′ = NTU [(1 − φ) w + D₀ − ΔT_p s], D₀ = Ts_out − T_r,in
  (closed form, finite at φ = 1), Tr = T_r,in + ΔT_p s + φ w.
- Superheat zone (xs … 1): counter-flow single-phase streams with the zone's U value and capacities
  (`evaporator_zones`): dTr/dt = NTU_r D, dTs/dt = NTU_s D, from Te and Ts_m.
"""
function evaporator_profiles(c::Evaporator, x, p; n::Integer = 100)
    zn = evaporator_zones(c, x, p)
    xs, Te, Tr_in, N, φ = x.xs, zn.twophase.Te, zn.twophase.Tr_in, zn.twophase.NTU, zn.twophase.φ
    glide = φ * x.Qe / zn.twophase.Cs
    ΔTp = Te - Tr_in - glide                       # saturation shift of the pressure drop
    s = range(0, 1, n)
    # two-phase zone: w = NTU [D₀ E₁(s) − ΔT_p E₂(s)], E₁ = ∫₀ˢ e^{βt}, E₂ = ∫₀ˢ E₁, β = (1 − φ) NTU
    β = (1 - φ) * N
    E1(t) = abs(β * t) < 1e-6 ? t + β * t^2 / 2 : expm1(β * t) / β
    E2(t) = abs(β * t) < 1e-6 ? t^2 / 2 + β * t^3 / 6 : (expm1(β * t) / β - t) / β
    D0 = p.Ts_out - Tr_in
    w = [N * (D0 * E1(t) - ΔTp * E2(t)) for t in s]
    Ts2 = p.Ts_out .+ w
    Tr2 = Tr_in .+ ΔTp .* s .+ φ .* w
    # superheat zone
    (; NTU_s, NTU_r) = zn.superheat
    k = NTU_s - NTU_r
    D1 = x.Ts_m - Te
    I(t) = abs(k) < 1e-12 ? D1 * t : D1 * expm1(k * t) / k                  # ∫₀ᵗ D
    Tr3 = Te .+ NTU_r .* I.(s)
    Ts3 = x.Ts_m .+ NTU_s .* I.(s)
    end_error = maximum(abs, (Ts2[end] - x.Ts_m, Tr2[end] - Te, Tr3[end] - (Te + c.ΔTsup), Ts3[end] - p.Ts_in))
    return (a = [xs .* s; xs .+ (1 - xs) .* s], Tr = [Tr2; Tr3], Ts = [Ts2; Ts3], xs = xs, end_error = end_error)
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
    z = _twophase_zone(c, x.Qe, x.xs, p.Ts_out, x.Ts_m, rs, _hei(c, p))
    m = ρsup * c.Vref
    Tm = Te + c.ΔTsup/2
    Δp_sup = vapour_pressure_drop(phx, m, ρ_pT(r, pe, Tm; phase = :vapor), μ_pT(r, pe, Tm; phase = :vapor), (1 - x.xs) * phx.L)
    return (Δp_twophase = z.Δp, Tr_in = z.Tr_in, Te = Te, γ = z.γ, NTU = z.NTU, Δp_superheat = Δp_sup)
end

"""
    BoreholeHX(; borehole, secondaryfluid, Vs, λground)

Borehole heat exchanger (Hellström, uniform borehole wall temperature `Tb`), as in
[`borehole_model!`](@ref).

- Ports: `Tin`, `Tout` (secondary fluid) and `Tb`, the borehole wall temperature (not the undisturbed ground `Tg`, see [`InfiniteLineSource`](@ref)):
  usually a fixed node, e.g. `nodes = (…, Tb = (:borehole => :Tb,)), fixed = (Tb = 278.15,)`.
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
        "`nodes = (…, Tb = (:borehole => :Tb,)), fixed = (Tb = $Tb,)`"))
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
- pinch (within `ΔTpinch`): `Ts_out - T_r,in` (two-phase zone, brine outlet against the refrigerant
  inlet: below `Te` for a blend, above it with a pressure drop) or `Ts_in - (Te + ΔTsup)` (superheat);
- a zone almost vanished: `xs < xs_tol` or `xs > 1 - xs_tol`.
"""
function validity(c::Evaporator, x, p; ΔTpinch = 0.1, xs_tol = 1e-3)
    issues = NamedTuple[]
    _secondary_issues!(issues, c.secondaryfluid, (Ts_in = p.Ts_in, Ts_m = x.Ts_m, Ts_out = p.Ts_out))
    Te = c.formulation isa NestedPressure ? T_dew(c.refrigerant, x.pe) : x.Te
    Te < c.Te_min + 0.01 &&
        push!(issues, _issue(:Te, Te, c.Te_min, "evaporation temperature at its lower bound Te_min"))
    Tr_in = evaporator_zones(c, x, p).twophase.Tr_in
    p.Ts_out - Tr_in < ΔTpinch &&
        push!(issues, _issue(:Te, Te, Te - (Tr_in - p.Ts_out) - ΔTpinch, "two-phase pinch: Ts_out - T_r,in below ΔTpinch"))
    p.Ts_in - (Te + c.ΔTsup) < ΔTpinch &&
        push!(issues, _issue(:Te, Te, p.Ts_in - c.ΔTsup - ΔTpinch, "superheat pinch: Ts_in - (Te + ΔTsup) below ΔTpinch"))
    x.xs < xs_tol && push!(issues, _issue(:xs, x.xs, xs_tol, "two-phase zone almost vanished"))
    x.xs > 1 - xs_tol && push!(issues, _issue(:xs, x.xs, 1 - xs_tol, "superheat zone almost vanished"))
    return issues
end

"Issues reported for a borehole: secondary fluid (`Tin`, `Tout`) below freezing or outside its table."
validity(c::BoreholeHX, x, p) = _secondary_issues!(NamedTuple[], c.secondaryfluid, (Tin = p.Tin, Tout = p.Tout))
