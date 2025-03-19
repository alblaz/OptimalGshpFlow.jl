"""
`PlateHeatExchanger` is a type that contains all the properties that describe plate heat exchangers

`A` total heat exchange surface area of the plate heat exchanger\\
`L` plate length\\
`W` plate width\\
`H` heat exchanger height\\
`np` number of plates\\
`ϕ` constant ... \\
`b` pressing depth\\
`β` chevron angle\\
"""
struct PlateHeatExchanger{T <: Number}
   A::T
   L::T
   W::T
   H::T
   np::Int
   ϕ::T
   b::T
   β::T
end

"""
Constructor 
===========

PlateHeatExchanger(A,L,W; b=0.002, ϕ=1.189, β = 30.)
"""
function PlateHeatExchanger(A,L,W; b=0.002, ϕ=1.189, β = 30.)
    np =  round(Int, A/(L*W)/ϕ) + 2        # number of plates
    H = (np-1)/2 * b                  # total height of the secondary fluid channels 
    return PlateHeatExchanger(A,L,W,H,np,ϕ,b,β)
end

"""
`h_platehex_wanniarachchi(V, phx_geometrical_param, properties; β = 30.)`

compute the convective heat transfer coefficient for a liquid flowing in a plate heat exchanger according to the correlation suggested by Wanniarachchi.

`V`                              volume flow rate, 
`T`                              temperature, 
`platehex::PlateHeatExchanger`   plate heat exchanger geometrical parameters, 
`properties`                     fluid properties, 
`β`                              angle of the groove in the plate heat exchanger. 
"""
function h_platehex_wanniarachchi(V, platehex::PlateHeatExchanger, fluidproperties)    
    
    A, L, W, H, np, ϕ, b, β = getfield.(Ref(platehex), fieldnames(PlateHeatExchanger))    
    ρ, μ, λ, cp = fluidproperties

    us = V/(H*W)                  # speed of the secondary fluid 
    dh = 2*b                       # hydraulic diameter 

    Redh = ρ * us * dh / μ           # Reynolds

    Pr = cp * μ /  λ
    jNul = 3.65/(β^0.455 * Redh^0.339)
    jNut = 12.6/(β^1.142) * Redh^(0.642 + 0.00111*β)
    jNu  = (jNul^3 + jNut^3)^(1/3) 
    Nudh = jNu * Pr^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
    h = Nudh*λ/dh

    return h #, Redh, us, np
end



# function platehex_wanniarachchi(V, T, platehex::PlateHeatExchanger, fluidproperties)    
    
#     A, L, H, np, ϕ, b, β = getfield.(Ref(platehex), fieldnames(PlateHeatExchanger))    
#     ρ, μ, λ, cp = fluidproperties.ρ(T), fluidproperties.μ(T), fluidproperties.λ(T), fluidproperties.cp(T)

#     us = V/(H*W)                  # speed of the secondary fluid 
#     dh = 2*b                       # hydraulic diameter 

#     Redh = ρ * us * dh / μ           # Reynolds
#     # @show Redh
#     Pr = cp * μ /  λ
#     jNul = 3.65/(β^0.455 * Redh^0.339)
#     jNut = 12.6/(β^1.142) * Redh^(0.642 + 0.00111*β)
#     jNu  = (jNul^3 + jNut^3)^(1/3) 
#     Nudh = jNu * Pr^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
#     h = Nudh*λ/dh

#     return h, Redh, us, Pr, jNul, JNut, jNu, Nudh
# end



# getfield.(Ref(my_obj), fieldnames(MyStruct))
#platehex = PlateHeatExchanger(1.,0.2, 0.06)
# A, L, H, np, ϕ, b, β = getfield.(Ref(platehex), fieldnames(PlateHeatExchanger))


###############################################
# extract heat exchanger temperature profiles #
###############################################

ϵ_counterflow(x,NTU,Cr) = (1 - exp(-NTU*(x)*(1-Cr)))/(1-Cr*exp(-NTU*(x)*(1-Cr)))

# the fluid with lower capacity Cc is flowing according to x (from 0 to 1), while the hot fluid with higher capacity is flowing on the contrary direction (from 1 to 0). This is important to keep in mind if we want to plot the inverse representation where the fluid with lower capacity enter the heat exchanger at x=1. In this case we first would perform the calculation above and then we would remap the the values of the

# function counterflow_heatexchanger(Ch,Thi,Tho,Cc,Tci,Tco,UA; x =0:0.1:1, defaultdirection = true) 
function counterflow_heatexchanger(heatexchanger_solution::HeatExchangerSolution; n=11 , defaultdirection = true) 

    x = range(0,1,n)
    @unpack Cc, Ch, NTU, Tch, Tcl, Thh, Thl, UA = heatexchanger_solution

    Thi, Tho = Thh, Thl
    Tci, Tco = Tcl, Tch

    Cmin = min(Cc,Ch)
    Cmax = max(Cc,Ch)
    # NTU = UA/Cmin
    Cr = Cmin/Cmax

    if Cc > Ch
        Th = @.  Thi - ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)
        Tc = @.  Tco - ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)*Ch/Cc
    else
        ϵx = ϵ_counterflow.(x,NTU,Cr)
        Tc = @.  Tci +  ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)
        Th = @.  Tho +  ϵ_counterflow(x,NTU,Cr) * (Thi - Tci)*Cc/Ch
    end

    if defaultdirection == false
        reverse!(Tc) 
        reverse!(Th)
    end

    return x,Tc,Th
end
##.
