#=
In this file we define data types to store the solution of variables by component.
=#

"""
    HeatExchangerSolution

Contains all the parameters of the heat exchanger model. The model is based on the ϵ-NTU model.

# Arguments
- `Q`: total heat exchnage
- `Thh`: Highest temperature of the warmer fluid
- `Thl`: Lowest temperature of the warmer fluid
- `Tch`: Highest temperature of the colder fluid
- `Tcl`: Lowest temperature of the colder fluid
- `Ch`: thermal capacity of the fluid at higher temperature
- `Cc`: thermal capacity of the fluid at colder temperature
- `Vh`: Volume flow rate of the fluid at higher temperature
- `Vc`: Volume flow rate of the fluid at lower temperature
- `hh`: heat transfer coefficient between the warm fluid and the pipe
- `hc`: heat transfer coefficient between the cold fluid and the pipe
- `U`: U value of the heat exchanger
- `A`: surface area of the heat exchanger
- `ϵ`: efficiency
- `NTU`: number of transfer units
"""
struct HeatExchangerSolution
    Q
    Thh
    Thl
    Tch
    Tcl
    Ch
    Cc
    Vh
    Vc
    hh
    hc
    UA
    ϵ
    NTU
end

"""
    BoreholeSolution

Borehole solution based on Hellström model

# Arguments
- `borehole::Borehole`: borehole configuration
- `Tin`: inlet borehole temperature
- `Tout`: output borehole temperature
- `Tb`: borehole temperature
- `Vs`: Volume flow rate of the fluid circulating in the borehole
- `Cs`: Capacity of the fluid circulating in the borehole
- `Qb`: Heat exchnaged in the borehole
- `Rbstar`: Effective borehole resistance
"""
struct BoreholeSolution
    borehole::Borehole
    λground
    Tin
    Tout
    Tb
    Vs
    Cs
    Qb
    Rbstar
end

"""
Constructor
===========

    BoreholeSolution(borehole,Tin,Tout,Tb,Vs,properties)

# Arguments
- `borehole::Borehole`: borehole configuration
- `Tin`: inlet borehole temperature
- `Tout`: output borehole temperature
- `Tout`: borehole temperature
- `Vs`: Volume flow rate of the fluid circulating in the borehole
- `Cs`: Capacity of the fluid circulating in the borehole
- `Qb`: Heat exchnaged in the borehole
- `Rbstar`: Effective borehole resistance
"""
function BoreholeSolution(borehole,λground,Tin,Tout,Tb,Vs,properties)
    ρ = properties.ρ((Tin+Tout)/2)
    cp = properties.cp((Tin+Tout)/2)
    Cs = ρ*cp*Vs
    Qb = Cs*(Tout - Tin)
    Rbstar = (Tb - (Tin+Tout)/2 )/(Qb/borehole.H)
    return BoreholeSolution(borehole, λground, Tin, Tout, Tb, Vs, Cs, Qb, Rbstar)
end


"""
    EvaporatorSolution

# Arguments
- `plateheatexchanger::PlateHeatExchanger`: contains geometrical parameters of the plate heat exchanger
- `xs`: 
- `Tout`: output borehole temperature
- `Tb`: borehole temperature
- `Vs`: Volume flow rate of the fluid circulating in the borehole
- `Cs`: Capacity of the fluid circulating in the borehole
- `Qb`: Heat exchnaged in the borehole
- `Rbstar`: Effective borehole resistance
"""
struct EvaporatorSolution
    plateheatexchanger::PlateHeatExchanger
    xs::Real
    twophase_region::HeatExchangerSolution
    superheat_region::HeatExchangerSolution
    ΔT_l  
    ΔT_m
    ΔT_h
end

function EvaporatorSolution(plateheatexchanger::PlateHeatExchanger, xs, twophase_region::HeatExchangerSolution, superheat_region::HeatExchangerSolution)
    ΔT_l = twophase_region.Thl - twophase_region.Tcl
    ΔT_m = twophase_region.Thh - twophase_region.Tch   
    ΔT_h = superheat_region.Thh - superheat_region.Tch  
    return EvaporatorSolution(plateheatexchanger,xs, twophase_region,superheat_region, ΔT_l, ΔT_m, ΔT_h)
end


function temperatureprofiles(evaporator_solution::EvaporatorSolution; n = 11)
    @unpack xs, ΔT_m, ΔT_l, ΔT_h =  evaporator_solution
    x_twophase,Tc_twophase,Th_twophase  =  counterflow_heatexchanger(evaporator_solution.twophase_region; n=n, defaultdirection = ΔT_l > ΔT_m )
    x_superheat,Tc_superheat,Th_superheat = counterflow_heatexchanger(evaporator_solution.superheat_region; n=n,  defaultdirection = ΔT_m > ΔT_h )
    x_twophase  = collect(x_twophase * xs)
    x_superheat = @. x_superheat * (1-xs) + xs
    x  = vcat(x_twophase, x_superheat)
    Tc = vcat(Tc_twophase, Tc_superheat)
    Th = vcat(Th_twophase, Th_superheat) 
    return x, Tc, Th
end

struct EvaporatorBoreholeSolution
    evaporator::EvaporatorSolution
    borehole::BoreholeSolution
end




# ΔT1,ΔT2 = Thle - Tcle, Thhe - Tche
# default_direction_e = ΔT1>ΔT2 ? true : false 
# @show default_direction_e
# chx(Che,Thi,Tho,Cce,Tci,Tco,UAe; x = xe , defaultdirection = default_direction_e) 

# ΔT1,ΔT2 = Thlsup - Tclsup, Thhsup - Tchsup
# default_direction_sup = ΔT1>ΔT2 ? true : false 

# @show default_direction_sup
# Thi,Tho,Tci,Tco = Thhsup,Thlsup,Tclsup,Tclsup
# xsuper, Thsup,Tcsup  = chx(Chsup,Thi,Tho,Ccsup,Tci,Tco,UAsup; x = xsup , defaultdirection = default_direction_sup) 
