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
    U
    A
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
end


struct EvaporatorBoreholeSolution
    evaporator::EvaporatorSolution
    borehole::BoreholeSolution
end