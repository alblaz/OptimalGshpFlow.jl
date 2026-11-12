# Nusselt number for fully developed flow in a circular pipe.
# Laminar (Re ≤ 2300): constant wall temperature value 3.657.
# Turbulent (Re ≥ 1e4): Dittus-Boelter, which requires knowing whether the fluid is heated or cooled
# (default: heated).
# Transition (2300 < Re < 1e4): linear interpolation in Re between the two values at the boundaries
# (Gnielinski 2013 / VDI Heat Atlas), so that Nu is continuous in Re. A discontinuous Nu makes the
# borehole equations lose their root when Re crosses a regime boundary.
const RE_LAMINAR   = 2300.
const RE_TURBULENT = 1e4
const NU_LAMINAR   = 3.657

nusselt_dittusboelter(Re, Pr; mode = "heating") = 0.023 * Re^0.8 * Pr^(mode == "heating" ? 0.3 : 0.4)

"""
`Nu = nusselt_circularpipe(Re, Pr; mode = "heating")`

`Nu`,`Re` and `Pr` are respectively the Nusselt, Reynold and the Prandtl numbers. The keyword argument takes either "heating" or "cooling".

- `Re ≤ 2300`: laminar, `Nu = 3.657`
- `Re ≥ 1e4`: turbulent, Dittus-Boelter
- in between: linear interpolation in `Re` between the laminar value and Dittus-Boelter at `Re = 1e4`
"""
function nusselt_circularpipe(Re, Pr; mode = "heating")
    if Re <= RE_LAMINAR
        return NU_LAMINAR + zero(Re * Pr)
    elseif Re >= RE_TURBULENT
        return nusselt_dittusboelter(Re, Pr; mode = mode)
    else
        γ = (Re - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR)
        return (1 - γ) * NU_LAMINAR + γ * nusselt_dittusboelter(RE_TURBULENT, Pr; mode = mode)
    end
end


# Evaluate the heat transfer coefficient  
"""
`h_circularpipe(mf, Tref, rp, fluid_properties)``

Compute the heat transfer coefficient for circular pipe

`mf`   mass flow rate\\
`Tref` reference temperature for the evaluation of properties\\
`rp`   pipe radius\\
`fluid_properties` is an object containing a set of functions `f(Tref)` that returns various fluid properties 
"""
function h_circularpipe(mf, Tref, rp, fluid_properties)
    # rp = get_rp(borehole)
    ρ= fluid_properties.ρ(Tref)
    μ = fluid_properties.μ(Tref)
    cp= fluid_properties.cp(Tref)
    k = fluid_properties.λ(Tref)
    w = mf/(ρ * π * rp^2)
    Re = 2 * ρ * w * rp / μ
    Pr = μ * cp / k
    Nu = nusselt_circularpipe(Re, Pr)
    return Nu * k /(2*rp)
end
