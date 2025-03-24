
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
    `f_Martin(Re)`

compute the Darcy friction factor fD given the Reynolds number and cheveron angle `β` based on Martin (1996)
"""
function f_Martin(Re,β)
    β=β*pi/180
    if Re < 2000
        f₀=64/Re
        f₁=597/Re+3.85
    else
        f₀=(1.8*log10(Re)-1.5)^(-2) # check log or log10
        f₁=39/Re^0.289
    end
    A=1/cos(β)*(0.18*sin(β)*(1+2*cos(β))+f₀)
    fD = 1/(cos(β)/(√(A))+(1-cos(β))/√(3.8*f₁))^2 # Martin 1996
    return fD
end

"""
    `pressure_drop_evap(V, platehex::PlateHeatExchanger, fluidproperties)`

Compute the pressure drop for a plate heat exchanger.

`V` is the volume flow rate (m^3/s),
`platehex::PlateHeatExchanger` is the plate heat exchanger object containing parameters such as number of plates, plate width, chevron angle, and plate length,
`fluidproperties` is an object of type `SecondaryFluid` that provides fluid properties like density and viscosity.
"""
function pressure_drop_evap(V, platehex::PlateHeatExchanger, fluidproperties,T;D_port=0.04,N_pass=1)    
    ρ = fluidproperties.ρ(T)
    μ = fluidproperties.μ(T)
    @unpack np,b,ϕ,β,L = platehex
    Nc=(np-1)/2 #number of channels
    Dₕ=2*b/ϕ # hydraulic diameter 4*A/P

    v = V/(L*b*Nc) # velocity per channel
    Re = ρ*v*Dₕ/μ #Reynolds number per channel
    fD = f_Martin(Re,β)

    # v_port = V/(π*D_port^2/4)
    # Δp_port=+ 1.5*ρ*v_port^2/2*N_pass
    # NB! port pressure drop accounted for in Martin
    
    Δp = fD * L/Dₕ * v^2 * ρ/2 
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