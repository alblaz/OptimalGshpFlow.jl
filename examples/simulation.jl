using OptimalGshpFlow
using CoolProp, NonlinearSolve, Parameters

### TO-DO
# assume a condensing pressure and get the compressor power after optimizing for the evap-borehole model
# circ. pump model based on pressure drop and flow rate
# what to show as a simulation:
# - different borehole depths (?), 
# - different collector diameter
# - different 

##.

# DESIGN NOMINAL CONDITIONS 
# Tc = 40°; Te = -5°C; hei =hco; Q2 = 5000 W; defaults ΔTs = 3.; U = 2000 W/m^2;  we discard superheating in the design

K=273.15 #K to °C
########################################
# input parameters for the simulation #
########################################
        ################################
        # evaporator model #
        ################################

        ΔTsup = 5.;
        ΔTsub = 5.;
        ηᵥ=0.95;
        # Tsf_h = 278.;
        ################################
        # borehole model #
        ################################

        H = 100.
        D = 4.
        λg = 3.
        pipe_position = ((0.04, 0.0), (-0.04, 0.0))
        λground = 3.5
        Tin, Tb  = 0. + K , 5. + K


########################################
# simulation #
########################################

################################
## preparation phase ##
################################

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
evaporator = EvaporatorDesign(2000., 0. +K, 40. +K, refrigerant_name,secondaryfluid_name;
        U =2000., ΔTs = 3., ΔT_l = 1.5, ΔT_sub=ΔTsub, ΔT_sup=ΔTsup, ηᵥ=ηᵥ)

Trange = range(-15. + K, 30 + K, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

# we have the sizing parameter for a generic evaporator (A,Vs,Vref), and we now need to choose the size L and W of the plate heat exchanger
# PlateHeatExchanger(A,L,W; b=0.002, ϕ=1.189, β = 30.) --> hence giving th enumber of plates as output
@unpack Vref,Vs,hei,ΔT_sup,Thh,A=evaporator
plateheatexchanger = PlateHeatExchanger(A, 0.2, 0.06; b=0.002, ϕ=1.189, β = 30.)


####################
# evaporator model #
####################


p  = EvaporatorParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid,
                                        plateheatexchanger=plateheatexchanger, 
                                        hei = hei, Vs = Vs, Vref = Vref, 
                                        ΔTsup= ΔT_sup, Ts_h=Thh)
## could we just pass an evaporator type here?


# initial guesses for the u vector
pme = tanhmap_ab_to_R(PropsSI("P", "T",  evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
ϵme = ϵmsup =  tanhmap_01_to_R(0.5) 
xms = tanhmap_01_to_R(0.5)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Tsf_l, Tsf_m = evaporator.Thh - 4, evaporator.Thh - 2

u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m]

prob = NonlinearProblem(evaporator_model!,u0,p,abstol = 1e-7)  
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporator_solution = evaporator_model(sol.u,p)


##################
# borehole model #
##################
borehole = SingleUPipeBorehole(H=H, D=D, λg=λg, pipe_position = pipe_position)
u0 = [2. + K]
p = BoreholeParameters(secondaryfluid=secondaryfluid, Vs=Vs, Tin=Tin, borehole=borehole, λground=λground, Tb=Tb)
prob = NonlinearProblem(borehole_model!,u0,p,abstol = 1e-9)  
sol  = solve(prob, TrustRegion(); abstol = 1e-9)


#############################
# evaporator-borehole model #
#############################
# The model evaporatorborehole_model! combine the models above by imposing Tin = Ts_l and Tout = Ts_h

p = EvaporatorBoreholeParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
                                        hei = hei, Vs = Vs, Vref = Vref, 
                                        ΔTsup =  ΔT_sup, 
                                        borehole=borehole, λground=λground, Tb=Tb)
u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, evaporator.Thh]
prob = NonlinearProblem(evaporatorborehole_model!,u0,p,abstol = 1e-7) 
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporatorborehole_solution = evaporatorborehole_model(sol.u,p)

# evap_xs=evaporatorborehole_solution.evaporator.xs
# evap_2phase=evaporatorborehole_solution.evaporator.twophase_region
# evap_superheat=evaporatorborehole_solution.evaporator.superheat_region
# evap_T=evap_2phase.Tcl - K

#############################
# Plotting #
#############################
# x_2phase = range(0., evap_xs, length = 40)

# UA=evap_2phase.A*evap_2phase.U
# C_sf=evap_2phase.Ch

# T_sf_out_2phase=evap_2phase.Thl

# T_profile_2phase = evap_T .+ (T_sf_out_2phase-evap_T)*exp.(UA ./ (x_2phase*C_sf))


x, Tc, Th = temperatureprofiles(evaporator_solution; n = 15)

W  = compressor_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)
Q1 = condenser_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)

using Plots
plotly()
plot(x, Th .- K, label="secondary fluid")
plot!(x, Tc .- K, label="refrigerant")

# compute pumping power in the broehole
@unpack rpi, H = evaporatorborehole_solution.borehole.borehole
@unpack Tin, Tout = evaporatorborehole_solution.borehole
@unpack Vc = evaporatorborehole_solution.evaporator.superheat_region
Ppump = pumping_power(0.4, H, Vc ,2rpi, secondaryfluid, (Tin + Tout)/2,evaporatorborehole_solution.evaporator.plateheatexchanger)
