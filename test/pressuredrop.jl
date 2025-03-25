
using OptimalGshpFlow, CoolProp, Parameters
# using NonlinearSolve, WGLMakie # This two packages are not dependencies so I load them from a separate environment.

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
evaporator = EvaporatorDesign(6000., 273.,273+40, refrigerant_name,secondaryfluid_name; U =2500., ΔTs = 3., ΔT_l = 1.5)

Trange = range(-15. + 273.15, 30 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.5, 0.2; b=0.002, ϕ=1.189, β = 30.)  # NB! chevron angle vs. reverse chev. ang. !

@unpack Vs,Thh,Thl = evaporator
T_mean=(Thh+Thl)/2
@show pressure_drop_evap(Vs, plateheatexchanger, secondaryfluid,T_mean)
@show pressure_drop(200*2,Vs,0.04*(1-2/17),secondaryfluid,T_mean) # SDR 17
@show pumping_power(0.5, 2*200, Vs, 0.04*(1-2/17), secondaryfluid, T_mean, plateheatexchanger)


plateheatexchanger2 = PlateHeatExchanger(evaporator.A, 0.5, 0.2; b=0.002, ϕ=1.189, β = 90. - 30.)  # NB! chevron angle vs. reverse chev. ang. !

@show pressure_drop_evap(Vs, plateheatexchanger2, secondaryfluid,T_mean)
