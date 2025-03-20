
"""
    `fDarcy(Re)`

compute the Darcy friction factor fD given the Reynolds number
"""
function fDarcy(Re)
    fD = Re<2300. ? 64/Re : (1/(lambertw(0.629*Re)))^2
    return fD
end


"""
    `pressure_drop(D,H,V,fluidproperties::SecondaryFluid,T)`

compute the the pressure drop for a circulat pipe

`D` is the pipe diameter,
`H` is the pipe Length,
`V` is the volume flow rate (m^3/s)
`fluidproperties::SecondaryFluid`,
`T` is the mean temperature along the pipe.
"""
function pressure_drop(H,V,D,fluidproperties,T)
    ρ = fluidproperties.ρ(T)
    μ = fluidproperties.μ(T)

    v = V/(π*D^2/4)
    Re = ρ*v*D/μ
    fD = fDarcy(Re) 
    
    Δp = fD * H/D * v^2 * ρ/2
    return Δp
end

"""
    `pumping_power(Δp, V, η)`

Compute the pumping power given the pressure loss `Δp`, 
the volume flow rate `V` and the pump efficiency `η`

    `pumping_power(η, H, V, D, fluidproperties, T)`

Compute the pumping power given 
`D` the pipe diameter,
`H` the pipe Length,
`V` the volume flow rate (m^3/s),
`fluidproperties::SecondaryFluid`,
`T` the mean temperature along the pipe.   
"""
function pumping_power(Δp, V, η)
    P = Δp * V / η
    return P
end

function pumping_power(η, H, V, D, fluidproperties, T)
    Δp = pressure_drop(H,V,D,fluidproperties,T)
    P = Δp * V / η
    return P
end