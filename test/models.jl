
#=
Sort of a unified interface/procedure to compute the solution for these three system types

EVAPORATOR MODEL
p = refrigerant, secondaryfluid, plateheatexchanger, hei, Vs, Vref, ΔTsup, Ts_h
u0 = pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m
prob = NonlinearProblem(evaporator_model!,u0,p,abstol = 1e-7)  #with evaporator_model!(du,u,p) 
@time sol = solve(prob, TrustRegion(); abstol = 1e-9) 
evaporator_solution = evaporator_model(sol.u,p)

BOREHOLE MODEL
p = secondaryfluid, Vs, Tin, borehole, λg, Tb
u0 = Tout
prob = NonlinearProblem(borehole_model!,u0,p,abstol = 1e-7)  
sol = solve(prob, TrustRegion(); abstol = 1e-9) 
borehole_solution = borehole_model(u,p)

EVAPORATORBOREHOLE MODEL
p = refrigerant, secondaryfluid, plateheatexchanger, borehole, hei, Vs, Vref, ΔTsup, borehole, λg, Tb
u0 = pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, Tsf_h
prob = NonlinearProblem(evaporatorborehole_model,u0,p,abstol = 1e-7)  #with evaporatorborehole_model(du,u,p) 
@time sol = solve(prob)
evaporatorborehole_solution = evaporatorborehole_model(sol.u,p)
=#


#####################
# preparation phase #
#####################

using OptimalGshpFlow, CoolProp, NonlinearSolve

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
evaporator = EvaporatorDesign(2000., 273.,273+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1.5)

Trange = range(-15. + 273.15, 30 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

# we have the sizing parameter for a generic evaporator (A,Vs,Vref), and we now need to choose the size L and W of the plate heat exchanger
# PlateHeatExchanger(A,L,W; b=0.002, ϕ=1.189, β = 30.)
plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.2, 0.06; b=0.002, ϕ=1.189, β = 30.)



####################
# evaporator model #
####################

# input parameters for the simulation
ΔTsup = 5.;  
Tsf_h = 278.;

p  = refrigerant, secondaryfluid, plateheatexchanger, evaporator.hei, evaporator.Vs, evaporator.Vref, ΔTsup, Tsf_h

pme = tanhmap_ab_to_R(PropsSI("P", "T",  evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
ϵme = ϵmsup =  tanhmap_01_to_R(0.5) 
xms = tanhmap_01_to_R(0.5)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Tsf_l, Tsf_m = Tsf_h - 4, Tsf_h - 2

u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m]

prob = NonlinearProblem(evaporator_model!,u0,p,abstol = 1e-7)  
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporator_solution = evaporator_model(sol.u,p)


##################
# borehole model #
##################
borehole = SingleUPipeBorehole(H = 100., D = 4., λg = 3., pipe_position = ((0.04, 0.0), (-0.04, 0.0)))
λground = 3.5
Tin, Tb  = 273.15, 278.15
Vs = evaporator.Vs
p = secondaryfluid, Vs, Tin, borehole, λground, Tb 
u0 = [275.]
prob = NonlinearProblem(borehole_model!,u0,p,abstol = 1e-9)  
sol  = solve(prob, TrustRegion(); abstol = 1e-9)


#############################
# evaporator-borehole model #
#############################
# The model evaporatorborehole_model! combine the models above by imposing Tin = Ts_l and Tout = Ts_h

p  = refrigerant, secondaryfluid, plateheatexchanger, evaporator.hei, evaporator.Vs, evaporator.Vref, ΔTsup, bh, λground , Tb
u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, Tsf_h]
prob = NonlinearProblem(evaporatorborehole_model!,u0,p,abstol = 1e-7) 
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporatorborehole_solution = evaporatorborehole_model(sol.u,p)





# uniform_Tb_coeffs(bh, 3.5, 0.5, 293.15, secondaryfluid)

# prob = NonlinearProblem(evaporator_model!,u0,p,abstol = 1e-7)  #with evaporator!(du,u,p) 
# nlcache = init(prob, TrustRegion())
# step!(nlcache)
# nlcache
# nlcache.J


# @time solve(prob, TrustRegion(); abstol = 1e-8)
# evaporator_solution = evaporatormodel(sol.u,p)

# u  = pme , evaporator.Te + 4. , Q2, 3000., 0.3, Te + 6, 500., 2000.

# p =  A, Vs, Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup, geometrical_param
# prob = NonlinearProblem(evaporator!,u0,p,abstol = 1e-7)  #with evaporator!(du,u,p) 
# @time sol = solve(prob)






# @testset "evaporator.jl" begin

# end