
using OptimalGshpFlow, CoolProp, Parameters
# using NonlinearSolve, WGLMakie # This two packages are not dependencies so I load them from a separate environment.

refrigerant_name, secondaryfluid_name = "R290", "water" #"INCOMP::MEA-28%"
evaporator = EvaporatorDesign(6000., 278.,273+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1)

Trange = range(0.1 + 273.15, 30 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

plateheatexchanger = PlateHeatExchanger(3.33, 0.52, 0.19; b=0.00206, ϕ=1.189, β = 30., Dport=0.06,N_pass=1)  # NB! chevron angle vs. reverse chev. ang. !
plateheatexchanger2 = PlateHeatExchanger(3.33, 0.52, 0.19; b=0.00206, ϕ=1.189, β = 30., Dport=0.06,N_pass=2)  # NB! chevron angle vs. reverse chev. ang. !

@unpack Vs,Thh,Thl = evaporator
T_mean=(Thh+Thl)/2
@show pressure_drop_evap(Vs, plateheatexchanger, secondaryfluid,T_mean)
@show pressure_drop_evap(Vs, plateheatexchanger2, secondaryfluid,T_mean)
@show pressure_drop(200*2,Vs,0.0354,secondaryfluid,T_mean) # SDR 17
@show pumping_power(0.4, 2*200, Vs, 0.0354, secondaryfluid, T_mean, plateheatexchanger)
@show Vs
# 0.04*(1-2/17)

plateheatexchanger2 = PlateHeatExchanger(evaporator.A,  0.52, 0.19; b=0.002, ϕ=1.189, β = 90. - 30.)  # NB! chevron angle vs. reverse chev. ang. !

@show pressure_drop_evap(Vs, plateheatexchanger2, secondaryfluid,T_mean)
