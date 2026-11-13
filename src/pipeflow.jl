# Nusselt number for fully developed flow in a circular pipe.
# Laminar (Re < 2300): constant wall temperature value 3.657.
# Turbulent (Re ≥ 1e4): Dittus-Boelter, which requires knowing whether the fluid is heated or cooled
# (default: heated).
# Transition (2300 ≤ Re < 1e4): depends on the transition model, see `NusseltTransition`.
const RE_LAMINAR   = 2300.
const RE_TURBULENT = 1e4
const NU_LAMINAR   = 3.657

"""
    NusseltTransition

Model for the Nusselt number in the laminar-turbulent transition range `2300 ≤ Re < 1e4`.

- `SmoothTransition()` (default): linear interpolation in `Re` between the laminar value and
  Dittus-Boelter at `Re = 1e4` (Gnielinski 2013 / VDI Heat Atlas). `Nu` is continuous in `Re`,
  so model solutions vary continuously with flow rate and temperature.
- `SharpTransition()`: Gnielinski correlation from `Re = 2300`, i.e. `Nu` jumps at `Re = 2300`
  (×2–7 depending on `Pr`) and at `Re = 1e4`. This is the original model of the package.

The same two types select the friction-factor model of the pipe pressure drop, see `fDarcy`
(smooth: linear blend between `Re = 2300` and `4000`; sharp: original switch at `Re = 2300`).
"""
abstract type NusseltTransition end
struct SmoothTransition <: NusseltTransition end
struct SharpTransition  <: NusseltTransition end

nusselt_dittusboelter(Re, Pr; mode = "heating") = 0.023 * Re^0.8 * Pr^(mode == "heating" ? 0.3 : 0.4)

function nusselt_gnielinski(Re, Pr)
    f = (1.58*log(Re) - 3.28)^(-2)     # Fanning friction factor
    return (0.5*f*(Re - 1000.) * Pr)/(1 + 12.7*sqrt(0.5*f) * (Pr^(2/3) - 1))
end

nusselt_transition(::SmoothTransition, Re, Pr; mode = "heating") =
    (γ = (Re - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR);
     (1 - γ) * NU_LAMINAR + γ * nusselt_dittusboelter(RE_TURBULENT, Pr; mode = mode))
nusselt_transition(::SharpTransition, Re, Pr; mode = "heating") = nusselt_gnielinski(Re, Pr)

"""
`Nu = nusselt_circularpipe(Re, Pr; mode = "heating", transition = SmoothTransition())`

`Nu`,`Re` and `Pr` are respectively the Nusselt, Reynold and the Prandtl numbers. `mode` takes either
"heating" or "cooling" (Dittus-Boelter exponent). `transition` selects the model in the
laminar-turbulent transition range, see [`NusseltTransition`](@ref).

- `Re < 2300`: laminar, `Nu = 3.657`
- `2300 ≤ Re < 1e4`: transition model
- `Re ≥ 1e4`: turbulent, Dittus-Boelter
"""
function nusselt_circularpipe(Re, Pr; mode = "heating", transition::NusseltTransition = SmoothTransition())
    if Re < RE_LAMINAR
        return NU_LAMINAR + zero(Re * Pr)
    elseif Re < RE_TURBULENT
        return nusselt_transition(transition, Re, Pr; mode = mode)
    else
        return nusselt_dittusboelter(Re, Pr; mode = mode)
    end
end


# Evaluate the heat transfer coefficient  
"""
`h_circularpipe(mf, Tref, rp, fluid_properties; transition = SmoothTransition())`

Compute the heat transfer coefficient for circular pipe

`mf`   mass flow rate\\
`Tref` reference temperature for the evaluation of properties\\
`rp`   pipe radius\\
`fluid_properties` is an object containing a set of functions `f(Tref)` that returns various fluid properties\\
`transition` laminar-turbulent transition model, see [`NusseltTransition`](@ref)
"""
function h_circularpipe(mf, Tref, rp, fluid_properties; transition::NusseltTransition = SmoothTransition())
    # rp = get_rp(borehole)
    ρ= fluid_properties.ρ(Tref)
    μ = fluid_properties.μ(Tref)
    cp= fluid_properties.cp(Tref)
    k = fluid_properties.λ(Tref)
    w = mf/(ρ * π * rp^2)
    Re = 2 * ρ * w * rp / μ
    Pr = μ * cp / k
    Nu = nusselt_circularpipe(Re, Pr; transition = transition)
    return Nu * k /(2*rp)
end
