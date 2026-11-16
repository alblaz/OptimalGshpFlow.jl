
const RE_FRICTION_TURBULENT = 4000.   # start of the turbulent smooth-pipe fit in `fDarcy`

fdarcy_laminar(Re) = 64 / Re
fdarcy_turbulent(Re) = (1 / (0.838 * lambertw(0.629 * Re)))^2   # smooth pipe (Colebrook fit)

"""
    fDarcy(Re, transition = SmoothTransition())

Darcy friction factor of fully developed flow in a smooth circular pipe.

- `Re < 2300`: laminar, `64/Re`
- `Re ≥ 4000`: turbulent smooth-pipe fit `(1/(0.838 W(0.629 Re)))²` (`W` Lambert function)
- in between, `transition` (see [`NusseltTransition`](@ref)):
  `SmoothTransition()` interpolates linearly in `Re` between the laminar value at 2300 and the
  turbulent one at 4000, so `fD` (and the pressure drop) is continuous in the flow rate;
  `SharpTransition()` is the original model: the turbulent fit from `Re = 2300`, a ×1.65 jump.
"""
function fDarcy(Re, transition::NusseltTransition = SmoothTransition())
    Re < RE_LAMINAR && return fdarcy_laminar(Re)
    (transition isa SharpTransition || Re ≥ RE_FRICTION_TURBULENT) && return fdarcy_turbulent(Re)
    γ = (Re - RE_LAMINAR) / (RE_FRICTION_TURBULENT - RE_LAMINAR)
    return (1 - γ) * fdarcy_laminar(RE_LAMINAR) + γ * fdarcy_turbulent(RE_FRICTION_TURBULENT)
end


"""
    `pressure_drop(H, V, D, fluidproperties::SecondaryFluid, T; transition = SmoothTransition())`

compute the pressure drop for a circular pipe

`H` is the pipe length,
`V` is the volume flow rate (m^3/s),
`D` is the pipe diameter,
`fluidproperties::SecondaryFluid`,
`T` is the mean temperature along the pipe,
`transition` is the laminar-turbulent transition model of the friction factor, see [`fDarcy`](@ref).
"""
function pressure_drop(H,V,D,fluidproperties,T; transition::NusseltTransition = SmoothTransition())
    ρ = fluidproperties.ρ(T)
    μ = fluidproperties.μ(T)

    v = V/(π*D^2/4)
    Re = ρ*v*D/μ
    fD = fDarcy(Re, transition)

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
    # A=1/cos(β)*(0.18*sin(β)*(1+2*cos(β))+f₀)
    # fD = 1/(cos(β)/(√(A))+(1-cos(β))/√(3.8*f₁))^2 # Martin 1996
    # fD = 1/(cos(β)/(√(A))+(1-cos(β))/√(3.8*f₁))^2 # Martin 1996
    fD= (cos(β)/(0.18*tan(β) + 0.36 *sin(β) + f₀/cos(β))  + (1-cos(β))/√(3.8*f₁))^-2
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
`T` the mean temperature along the pipe,
`transition` the friction-factor transition model (keyword), see [`fDarcy`](@ref).
"""
function pumping_power(Δp, V, η)
    P = Δp * V / η
    return P
end

function pumping_power(η, H, V, D, fluidproperties, T; transition::NusseltTransition = SmoothTransition())
    Δp = pressure_drop(H,V,D,fluidproperties,T; transition = transition)
    P = Δp * V / η
    return P
end