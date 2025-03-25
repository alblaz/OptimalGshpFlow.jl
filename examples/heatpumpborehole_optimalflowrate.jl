####################
# evaporator model #
####################
using OptimalGshpFlow, CoolProp, Parameters, LinearAlgebra

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
evaporator = EvaporatorDesign(6000., 273.,273+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1.5)

Trange = range(-15. + 273.15, 30 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

# we have the sizing parameter for a generic evaporator (A,Vs,Vref), and we now need to choose the size L and W of the plate heat exchanger
plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.5, 0.15; b=0.002, ϕ=1.189, β = 30.)  # NB! chevron angle vs. reverse chev. ang. !

##################
# borehole model #
##################
di = 0.04*(1-1/17) #internal pipe diameter
borehole = SingleUPipeBorehole(H = 200., D = 4., rpi = di/2,  rb = 0.115, λg = 2., pipe_position = ((0.035, 0.0), (-0.035, 0.0)))
@assert norm([borehole.pipe_position[1]...]) + borehole.rpo < borehole.rb  #check that at least the pipes are inside the borehole
@assert norm([borehole.pipe_position[2]...]) + borehole.rpo < borehole.rb  #check that at least the pipes are inside the borehole

λground = 3.

Tb = 273.15 - 5
ΔTsup = 0.5
ηpump = 0.3

# param_borehole_prototype = BoreholeParameters(secondaryfluid=secondaryfluid, Vs=Vs, Tin=Tin, borehole=borehole, λground=λground, Tb=Tb)
# u0 = [275.]
# prob = NonlinearProblem(borehole_model!,u0,param_borehole_prototype,abstol = 1e-9)  
# sol  = solve(prob, TrustRegion(); abstol = 1e-9)

#############################
# evaporator-borehole model #
#############################
# The model evaporatorborehole_model! combine the models above by imposing Tin = Ts_l and Tout = Ts_h

param_evaporatorborehole_prototype  = EvaporatorBoreholeParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
                                        hei = evaporator.hei, Vs = evaporator.Vs, Vref = evaporator.Vref, 
                                        ΔTsup =  ΔTsup, 
                                        borehole=borehole, λground=λground, Tb=Tb)

                                        
pme = sigmoid_ab_to_R(PropsSI("P", "T",  evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
ϵme = ϵmsup =  sigmoid_01_to_R(0.3) 
xms = sigmoid_01_to_R(0.3)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Ts_l, Ts_m, Ts_h = Tb - 6, Tb - 4, Tb - 2
u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m, Ts_h]


function solve_configuration(u0,p)
    prob = NonlinearProblem(evaporatorborehole_model!,u0,p,abstol = 1e-7) 
    sol  = solve(prob, TrustRegion(); abstol = 1e-9)
    return sol
end

params = [EvaporatorBoreholeParameters(param_evaporatorborehole_prototype,Vs = V) for V in  evaporator.Vs .* (2:-0.05:0.1)]
# sols = [solve_configuration(u0,p) for p in params]

function iteratively_solve(u0,params)
    sols = []
    for (k,p) in enumerate(params)
        current_sol = solve_configuration(u0,p)
        u0 =  current_sol.u   
        push!(sols, current_sol)
    end
    return sols
end

sols = iteratively_solve(u0,params)
retcodes = [sol.retcode for sol in sols]
@assert sum(x -> x == ReturnCode.Success, retcodes)  == length(retcodes)


function evaluate_solution(sol,p; Tcond = 273.15 + 40, ηis = 0.7)
    evaporatorborehole_solution = evaporatorborehole_model(sol.u,p)

    W  = compressor_power(Tcond, ηis, evaporatorborehole_solution.evaporator, refrigerant)
    Q1 = condenser_power(Tcond, ηis, evaporatorborehole_solution.evaporator, refrigerant)
    Q2 = evaporatorborehole_solution.evaporator.twophase_region.Q + evaporatorborehole_solution.evaporator.superheat_region.Q

    # compute pumping power in the borehole
    @unpack rpi, H = evaporatorborehole_solution.borehole.borehole
    @unpack Tin, Tout = evaporatorborehole_solution.borehole
    @unpack Vh = evaporatorborehole_solution.evaporator.superheat_region
    Ppump_borehole = pumping_power(ηpump, 2H, Vh ,2rpi, secondaryfluid, (Tin + Tout)/2 )
    #compute pumping power in the evaporator
    Tmean = (evaporatorborehole_solution.evaporator.superheat_region.Thh  + evaporatorborehole_solution.evaporator.twophase_region.Thl) / 2
    Δp_evap = pressure_drop_evap(evaporatorborehole_solution.evaporator.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean)
    Ppump_evap = pumping_power(Δp_evap, evaporatorborehole_solution.evaporator.superheat_region.Vh, ηpump)
    #compute total pumping Power
    Ppump = Ppump_borehole + Ppump_evap
    
    # COP1 = Q1/(W + Ppump)
    # COP2 = Q2/(W  + Ppump)
    COP1 = Q1/(W + Ppump)
    COP2 = Q2/(W  + Ppump)
    return COP1, COP2
end


COP1 = [evaluate_solution(sol,p)[1] for (sol,p) in zip(sols,params)]
COP2 = [evaluate_solution(sol,p)[2] for (sol,p) in zip(sols,params)]
Vs   = [p.Vs for p in params]

fig = Figure(; size = (600, 400))
ga = fig[1, 1] = GridLayout()
ax = Axis(fig[1, 1], xlabel = "Vs (l/s)", ylabel = "COP")
lines!(ax, Vs * 1000, COP1, color = :dodgerblue, linewidth = 2, label = "COP1")
lines!(ax, Vs * 1000, COP2, color = :indianred, linewidth = 2, label = "COP2")

fig