
using OptimalGshpFlow, CoolProp, Parameters, Revise
# using NonlinearSolve, WGLMakie # This two packages are not dependencies so I load them from a separate environment.

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
evaporator = EvaporatorDesign(2000., 273.,273+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1.5)

Trange = range(-15. + 273.15, 30 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.2, 0.06; b=0.002, ϕ=1.189, β = 30.)
# include("..\\pressuredrop.jl")
T_mean=(evaporator.Thh+evaporator.Thl)/2
pressure_drop_evap(evaporator.Vs, plateheatexchanger, secondaryfluid,T_mean)
pressure_drop(300*2,evaporator.Vs,0.02,secondaryfluid,T_mean)