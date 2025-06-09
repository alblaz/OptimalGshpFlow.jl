# evaluate nusselt number for flow in pipe for given value of Reynolds and Prandtl number
# currently for Reynolds greater than 10000 we use Dittus-Boelter correlation which require
# knowledge on whether the fluid is heated or cooled. The choosen default is that the fluid is heated.
"""
`Nu = circularpipe_nusselt(Re, Pr; mode = "heating")``

`Nu`,`Re` and `Pr` are respectively the Nusselt, Reynold and the Prandtl numbers. The keyword argument takes either "heating" or "cooling"

"""
function nusselt_circularpipe(Re, Pr; mode = "heating")
    if Re < 2300
        Nu = 3.657 # consider the entry length (now underestimated Nu)
    elseif (Re >= 2300) && (Re < 10000)        
        f =  (1.58*log(Re) - 3.28)^(-2)
        Nu = (0.5*f*(Re - 1000.) * Pr)/(1 + 12.7*sqrt(0.5*f) * (Pr^(2/3) - 1))
    else 
        n = mode == "heating" ? 0.3 : 0.4  
        Nu = 0.023 * Re^0.8 * Pr^n
    end
    return Nu
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
