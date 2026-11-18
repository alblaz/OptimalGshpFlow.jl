#= Refrigerant-side pressure drop in a brazed plate heat exchanger (ROADMAP step 5, design note
`docs/dev/glide.md`): two-phase (boiling) and single-phase vapour zones.

Total two-phase pressure drop (Amalfi et al. 2016, part 1, eq. 8): friction + acceleration
(momentum) + gravity + ports/manifolds, the last three from the homogeneous model. Both fluids
flow in channels of the same geometry (stacked brazed plates): `PlateHeatExchanger` serves both
sides; `Nc = (np − 1)/2` channels per side, each of cross-section `W b`.
=#

const GRAVITY = 9.81

"""
    TwoPhaseFriction

Frictional pressure-drop model of a boiling refrigerant in a plate heat exchanger:
[`Amalfi2016`](@ref) (general, multi-fluid) or [`KineticEnergyFriction`](@ref) (Longo-type
fits, linear in the kinetic energy per unit volume).
"""
abstract type TwoPhaseFriction end

"""
    Amalfi2016()

Amalfi, Vakili-Farahani & Thome (2016), Int. J. Refrig. 61, 185–203 (part 2, eqs. 7–8) with the
Fanning form of part 1 (eq. 18):

`Δp_f = 2 f_tp L G²/(d_h ρ_m)`, `f_tp = C 15.698 We_m^−0.475 Bd^0.255 (ρ_l/ρ_v)^−0.571`,
`C = 2.125 (β/70°)^9.993 + 0.955`, `We_m = G² d_h/(ρ_m σ)`, `Bd = (ρ_l − ρ_v) g d_h²/σ`,
`d_h = 2b/ϕ`, `ρ_m` homogeneous at the mean quality. Database: 1513 points, mean absolute error
21.5 %; refrigerants at saturation temperatures −3…35 °C.
"""
struct Amalfi2016 <: TwoPhaseFriction end

"""
    KineticEnergyFriction(ζ; L_ref = 0.278, dh_ref = 0.004)

Frictional pressure drop linear in the kinetic energy per unit volume (Longo and co-workers):
`Δp_f = ζ KE/V`, `KE/V = G²/(2ρ_m)`, `ρ_m` homogeneous at the mean quality, `ζ` in Pa per
J m⁻³ (1000 × the constants of the papers, which give `Δp_f` in kPa) for a plate of flow length
`L_ref` and hydraulic diameter `dh_ref = 2b`. For another plate (and a zone of length `L`) it is
scaled as an equivalent constant friction factor: `ζ L/(2b) ⋅ dh_ref/L_ref` (an assumption).

Fits on the same BPHE (L = 278 mm, W = 72 mm, chevron 65°, b = 2 mm): [`LONGO_FITS`](@ref).
"""
struct KineticEnergyFriction{T <: Real} <: TwoPhaseFriction
    ζ::T
    L_ref::T
    dh_ref::T
end
KineticEnergyFriction(ζ; L_ref = 0.278, dh_ref = 0.004) = KineticEnergyFriction(promote(float(ζ), L_ref, dh_ref)...)

"""
    LONGO_FITS

Longo-type kinetic-energy fits on the same BPHE (`ζ` in Pa per J m⁻³):
- `hfc_vaporisation`: R134a, R410A, R236fa vaporisation, Longo & Gasparella 2007 (ζ = 1490);
- `r404a_vaporisation`: R404A vaporisation at 20 °C, Longo et al. 2016 (ζ = 1667);
- `hc_condensation`: R600a, R290, R1270 condensation, Longo 2010 (ζ = 1900).
The R290 vaporisation fit (Longo 2012, ASME J. Heat Transfer 134, 101801) is not included.
"""
const LONGO_FITS = (hfc_vaporisation = KineticEnergyFriction(1490.),
                    r404a_vaporisation = KineticEnergyFriction(1667.),
                    hc_condensation = KineticEnergyFriction(1900.))

"Homogeneous two-phase density at the vapour quality `x`."
homogeneous_density(ρl, ρv, x) = inv(x/ρv + (1 - x)/ρl)

"Mass flux [kg/(m² s)] of the mass flow `m` in the channels of one side of `phx`."
channel_mass_flux(phx::PlateHeatExchanger, m) = m / (phx.W * phx.b * (phx.np - 1)/2)

"""
    friction_pressure_drop(model::TwoPhaseFriction, phx, G, xm, ρl, ρv, σ, L)

Two-phase frictional pressure drop [Pa] over the flow length `L` at mass flux `G`, mean quality
`xm`, saturated densities `ρl`, `ρv` and surface tension `σ`.
"""
function friction_pressure_drop(::Amalfi2016, phx::PlateHeatExchanger, G, xm, ρl, ρv, σ, L)
    dh = 2phx.b / phx.ϕ
    ρm = homogeneous_density(ρl, ρv, xm)
    We = G^2 * dh / (ρm * σ)
    Bd = (ρl - ρv) * GRAVITY * dh^2 / σ
    C = 2.125 * (phx.β / 70)^9.993 + 0.955
    f = C * 15.698 * We^(-0.475) * Bd^0.255 * (ρl / ρv)^(-0.571)
    return 2f * L * G^2 / (dh * ρm)
end

function friction_pressure_drop(m::KineticEnergyFriction, phx::PlateHeatExchanger, G, xm, ρl, ρv, σ, L)
    ρm = homogeneous_density(ρl, ρv, xm)
    ζ = m.ζ * (L / (2phx.b)) * (m.dh_ref / m.L_ref)
    return ζ * G^2 / (2ρm)
end

"Acceleration (momentum) pressure drop [Pa], homogeneous model, quality `x_in` → `x_out`."
acceleration_pressure_drop(G, x_in, x_out, ρl, ρv) = G^2 * (x_out - x_in) * (1/ρv - 1/ρl)

"Gravity pressure drop [Pa] over the height `L` for upward flow, homogeneous density `ρm`."
gravity_pressure_drop(ρm, L) = ρm * GRAVITY * L

"""
    port_pressure_drop(m, ρ_in, ρ_out; D_port = 0.04, N_pass = 1)

Inlet and outlet ports and manifolds (Shah & Focke 1988, as in Amalfi et al. 2016 part 1, eq. 11):
`0.75 N_pass [G_p²/(2ρ_in) + G_p²/(2ρ_out)]`, `G_p` the mass flux in a port of diameter `D_port`.
"""
function port_pressure_drop(m, ρ_in, ρ_out; D_port = 0.04, N_pass = 1)
    Gp = m / (π * D_port^2 / 4)
    return 0.75 * N_pass * (Gp^2 / (2ρ_in) + Gp^2 / (2ρ_out))
end

"""
    vapour_pressure_drop(phx, m, ρ, μ, L)

Single-phase (superheated vapour) frictional pressure drop [Pa] over the length `L`: Martin's
friction factor ([`f_Martin`](@ref), as on the brine side), `d_h = 2b/ϕ`.
"""
function vapour_pressure_drop(phx::PlateHeatExchanger, m, ρ, μ, L)
    dh = 2phx.b / phx.ϕ
    G = channel_mass_flux(phx, m)
    Re = G * dh / μ
    return f_Martin(Re, phx.β) * L / dh * G^2 / (2ρ)
end

"""
    twophase_pressure_drop(model, r::AbstractRefrigerant, phx, m, p, x_in, x_out; L = phx.L, gravity = true)

Pressure drop [Pa] of the refrigerant mass flow `m` boiling from quality `x_in` to `x_out` over the
flow length `L` at the pressure `p` (properties at saturation): `(; friction, acceleration,
gravity, total)`; `gravity = false` for horizontal flow (upward flow otherwise). Ports are not
included (see [`port_pressure_drop`](@ref)).
"""
function twophase_pressure_drop(model::TwoPhaseFriction, r::AbstractRefrigerant, phx::PlateHeatExchanger, m, p,
                                x_in, x_out; L = phx.L, gravity::Bool = true)
    ρl, ρv, σ = ρ_bubble(r, p), ρ_dew(r, p), surface_tension(r, p)
    G = channel_mass_flux(phx, m)
    xm = (x_in + x_out) / 2
    fr = friction_pressure_drop(model, phx, G, xm, ρl, ρv, σ, L)
    ac = acceleration_pressure_drop(G, x_in, x_out, ρl, ρv)
    gr = gravity ? gravity_pressure_drop(homogeneous_density(ρl, ρv, xm), L) : zero(fr)
    return (friction = fr, acceleration = ac, gravity = gr, total = fr + ac + gr)
end
