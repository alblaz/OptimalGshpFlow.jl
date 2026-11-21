#=
Refrigerant-side heat transfer during condensation in a brazed plate heat exchanger (step 9b).

Longo 2010 (Int. J. Refrig. 33, 944, `assets/docs/`): R290, R600a and R1270 condensing in a BPHE
(65° chevron, b = 2 mm). Two regimes, with the transition at an equivalent Reynolds number of
about 1600–1700 (G ≈ 15–18 kg/(m² s)):
- gravity controlled: Nusselt (1916) film condensation on a vertical surface (14.8 % mean
  deviation from the data);
- forced convection: Akers et al. (1959), local in the vapour quality (21 %).
Longo refers h and q to the projected plate area (and multiplies by the enlargement factor Φ);
here both are referred to the actual area `A` of the `PlateHeatExchanger`, as in the evaporator,
so Φ drops out: h_projected = Φ h, q_projected = Φ q.
Default (`LongoJH`): Longo's single curve of the heat transfer factor J_H(Re_eq) (Fig. 6) instead
of the two correlations, which jump at the transition (Akers over-predicts just above it).
=#

"""
    CondensationModel

Refrigerant-side condensation heat-transfer model: [`LongoJH`](@ref) (default), [`Longo2010`](@ref).
"""
abstract type CondensationModel end

"""
    LongoJH(; J_gravity = 55.8, Re_t = 1648.0, n = 0.535, δ = 0.1, Φ = 1.24)

Condensation in a brazed plate heat exchanger from Longo's heat transfer factor
`J_H = Nu/Pr_l^⅓` (`Nu = h d_h/λ_l`, `d_h = 2b`) against the equivalent Reynolds number at the mean
quality (Longo 2010, Fig. 6: one curve for R600a, R290 and R1270 at 25–40 °C, Re_eq ≈ 700–3500;
in agreement with his HFC data): constant `J_gravity` in the gravity-controlled region, the power
law `J_gravity (Re_eq/Re_t)^n` above the transition `Re_t`, joined by a C¹ fillet in log-log
coordinates over `ln(Re_eq/Re_t) ∈ [−δ, δ]` (at most `n δ/4` = 1.3 % above the corner).

Read off Fig. 6 (pixel trace, 2026-11-21): plateau 55.8 for Re_eq 1000–1600, 76.9 at 3000, 82.2
at 3400; the transition line at 1648. The power law through the two upper anchors (exponent 0.535)
meets the plateau at the transition. The exponent is above the "0.3–0.4" of Longo's text, but
consistent with his h ∝ G^0.65–0.7 above the transition (35–40 % for +60 % of G, Fig. 4).
Longo's h is per projected area: divided by his enlargement factor `Φ` for the actual area.

Unlike [`Longo2010`](@ref) (Nusselt below the transition, Akers above, joined over 1500–1800), h is
flat in the gravity region, as Longo's data, and does not jump to Akers (which over-predicts the
data by 16–42 % just above the transition): that jump made h fall by 20 % for −4 % of the mass flow
in the closed cycle (condenser at Re_eq 1580–1840).
"""
struct LongoJH <: CondensationModel
    J_gravity::Float64
    Re_t::Float64
    n::Float64
    δ::Float64
    Φ::Float64
end
LongoJH(; J_gravity = 55.8, Re_t = 1648.0, n = 0.535, δ = 0.1, Φ = 1.24) =
    LongoJH(float(J_gravity), float(Re_t), float(n), float(δ), float(Φ))

"""
    longo_jh(model::LongoJH, Re_eq)

Heat transfer factor `J_H = Nu/Pr_l^⅓` of [`LongoJH`](@ref) at the equivalent Reynolds number `Re_eq`.
"""
function longo_jh(model::LongoJH, Re)
    s = log(Re / model.Re_t)
    δ, n = model.δ, model.n
    y = s <= -δ ? zero(s) : s >= δ ? n * s : n * (s + δ)^2 / (4δ)
    return model.J_gravity * exp(y)
end

"""
    Longo2010(; Re_low = 1500.0, Re_high = 1800.0)

Condensation in a brazed plate heat exchanger after Longo 2010: Nusselt film condensation
([`nusselt_film_htc`](@ref), gravity controlled) below the transition and Akers forced-convection
condensation ([`akers_htc`](@ref)) above it, with Longo's criterion: the equivalent Reynolds number
at the mean quality, transition at 1600–1700 for all his fluids. The two are joined by a C¹ blend
(cubic Hermite weight) between `Re_low` and `Re_high`, as the solver needs a differentiable model.
(A smooth maximum of the two is not used: below the transition it takes Akers, which over-predicts
there, and adds the two where they coincide, up to +40 % against Longo's R290 data.)
"""
struct Longo2010 <: CondensationModel
    Re_low::Float64
    Re_high::Float64
end
Longo2010(; Re_low = 1500.0, Re_high = 1800.0) = Longo2010(float(Re_low), float(Re_high))

"""
    nusselt_film_htc(λl, ρl, μl, Δh_lv, L, q)

Nusselt (1916) laminar film condensation on a vertical surface of height `L` [m], mean
heat-transfer coefficient [W/(m² K)] written with the heat flux `q` [W/m²] instead of the wall
superheat ΔT = q/h: `h = 0.943 [λ³ρ²gΔh_lv/(μ L ΔT)]^¼` becomes `h = [0.943⁴ λ³ρ²gΔh_lv/(μ L q)]^⅓`
(liquid conductivity `λl`, density `ρl`, viscosity `μl`; latent heat `Δh_lv`). Falls with the
heat flux (thicker film).
"""
nusselt_film_htc(λl, ρl, μl, Δh_lv, L, q) = (0.943^4 * λl^3 * ρl^2 * GRAVITY * Δh_lv / (μl * L * q))^(1/3)

"""
    akers_htc(G, x_in, x_out, ρl, ρv, μl, λl, cpl, dh)

Akers et al. (1959) forced-convection condensation, `h = 5.03 (λl/dh) Re_eq^⅓ Pr_l^⅓` with
`Re_eq = G [(1 − x) + x (ρl/ρv)^½] dh/μl`, averaged over the vapour quality from `x_in` to
`x_out` (uniform heat flux: x linear in the area; Re_eq linear in x, closed form).
"""
function akers_htc(G, x_in, x_out, ρl, ρv, μl, λl, cpl, dh)
    a = G * dh / μl                           # Re_eq = a + b x
    b = a * (sqrt(ρl / ρv) - 1)
    Δx = x_in - x_out
    Re13 = abs(b * Δx) < 1e-8 * a ? cbrt(a + b * (x_in + x_out) / 2) :
           3 * ((a + b * x_in)^(4/3) - (a + b * x_out)^(4/3)) / (4b * Δx)
    return 5.03 * λl / dh * Re13 * cbrt(cpl * μl / λl)
end

"Equivalent Reynolds number of Akers/Longo at the vapour quality `x`."
reynolds_equivalent(G, x, ρl, ρv, μl, dh) = G * ((1 - x) + x * sqrt(ρl / ρv)) * dh / μl

"""
    condensation_htc(model::Longo2010, r, phx, m, p, x_in, x_out, q; L = phx.L) -> NamedTuple

Mean refrigerant-side heat-transfer coefficient [W/(m² K), actual area] of the refrigerant mass
flow `m` condensing at the pressure `p` from the quality `x_in` to `x_out` with the heat flux `q`
[W/m², actual area] over the flow length `L` of the zone, in the plate exchanger `phx` (`d_h = 2b`,
as Longo). Returns `(; h, gravity, forced, Re_eq)` with `Re_eq` at the mean quality.

Zeotropic blend: Silver–Bell–Ghaly, `1/h = 1/h_film + Z_g/h_g` with `Z_g = x cp_g dT/dh` (the
sensible cooling of the vapour as a fraction of the total heat, `dT/dh = (T_dew − T_bubble)/Δh_lv`)
and `h_g` the vapour flowing alone (Wanniarachchi), at the mean quality.
"""
function condensation_htc(model::Longo2010, r::AbstractRefrigerant, phx::PlateHeatExchanger, m, p, x_in, x_out, q;
                          L = phx.L)
    ρl, ρv = ρ_bubble(r, p), ρ_dew(r, p)
    μl, λl, cpl = μ_bubble(r, p), λ_bubble(r, p), cp_bubble(r, p)
    Δh_lv = h_dew(r, p) - h_bubble(r, p)
    G = channel_mass_flux(phx, m)
    dh = 2phx.b
    hg = nusselt_film_htc(λl, ρl, μl, Δh_lv, L, q)
    hf = akers_htc(G, x_in, x_out, ρl, ρv, μl, λl, cpl, dh)
    xm = (x_in + x_out) / 2
    Re = reynolds_equivalent(G, xm, ρl, ρv, μl, dh)
    t = clamp((Re - model.Re_low) / (model.Re_high - model.Re_low), zero(Re), one(Re))
    w = t^2 * (3 - 2t)                       # C¹ weight of the forced-convection regime
    h = _sbg_correction(r, phx, m, p, xm, (1 - w) * hg + w * hf, ρv, Δh_lv)
    return (h = h, gravity = hg, forced = hf, Re_eq = Re)
end

"""
    condensation_htc(model::LongoJH, r, phx, m, p, x_in, x_out, q; L = phx.L) -> NamedTuple

As the [`Longo2010`](@ref) method (same arguments; `q` and `L` are not used: h does not depend on
them), from Longo's heat transfer factor. Returns `(; h, film, J_H, Re_eq)`: `film` before the
Silver–Bell–Ghaly correction of blends, `J_H` and `Re_eq` at the mean quality.
"""
function condensation_htc(model::LongoJH, r::AbstractRefrigerant, phx::PlateHeatExchanger, m, p, x_in, x_out, q;
                          L = phx.L)
    ρl, ρv = ρ_bubble(r, p), ρ_dew(r, p)
    μl, λl, cpl = μ_bubble(r, p), λ_bubble(r, p), cp_bubble(r, p)
    Δh_lv = h_dew(r, p) - h_bubble(r, p)
    G = channel_mass_flux(phx, m)
    dh = 2phx.b
    xm = (x_in + x_out) / 2
    Re = reynolds_equivalent(G, xm, ρl, ρv, μl, dh)
    J = longo_jh(model, Re)
    hf = J * λl * cbrt(cpl * μl / λl) / (dh * model.Φ)
    h = _sbg_correction(r, phx, m, p, xm, hf, ρv, Δh_lv)
    return (h = h, film = hf, J_H = J, Re_eq = Re)
end

# Silver–Bell–Ghaly for blends (pure fluids: no correction, no cost)
_sbg_correction(::AbstractRefrigerant, phx, m, p, xm, h, ρv, Δh_lv) = h
function _sbg_correction(r::ClapeyronRefrigerant, phx, m, p, xm, h, ρv, Δh_lv)
    ispure(r) && return h
    cpv = cp_dew(r, p)
    Zg = xm * cpv * (T_dew(r, p) - T_bubble(r, p)) / Δh_lv
    hv = h_platehex_wanniarachchi(xm * m / ρv, phx, (ρv, μ_dew(r, p), λ_dew(r, p), cpv))
    return inv(inv(h) + Zg / hv)
end
