#=
effect of change of the secondary fluid flow rate
=#


#####################
# preparation phase #
#####################
using OptimalGshpFlow, NonlinearSolve, CoolProp, Parameters
# using NonlinearSolve, WGLMakie # This two packages are not dependencies so I load them from a separate environment.

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
# evaporator = EvaporatorDesign(5000., 273. -5 ,273+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1.5)
evaporator = EvaporatorDesign(6000., 273.,273+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1.5)


Trange = range(-15. + 273.15, 30 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

# we have the sizing parameter for a generic evaporator (A,Vs,Vref), and we now need to choose the size L and W of the plate heat exchanger
# PlateHeatExchanger(A,L,W; b=0.002, ϕ=1.189, β = 30.)
# plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.2, 0.06; b=0.002, ϕ=1.189, β = 30.)
# plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.25, 0.12; b=0.002, ϕ=1.189, β = 30.)
plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.5, 0.15; b=0.002, ϕ=1.189, β = 30.)  # NB! chevron angle vs. reverse chev. ang. !

####################
# evaporator model #
####################

# input parameters for the simulation
ΔTsup = 5.;  
Ts_h = 273.15 + 10;
ηpump = 0.3

# p  = refrigerant, secondaryfluid, plateheatexchanger, evaporator.hei, evaporator.Vs, evaporator.Vref, ΔTsup, Tsf_h

param_evaporator_prototype  = EvaporatorParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
                                        hei = evaporator.hei, Vs = evaporator.Vs, Vref = evaporator.Vref, 
                                        ΔTsup =  ΔTsup, Ts_h = Ts_h)

# you can use the template and just change one or two parameters instead of redifining everything all the time
# p2_ev = EvaporatorParameters(param_evaporator_prototype, Vs = 2*evaporator.Vs)
                                        
# @unpack plateheatexchanger = param_evaporator_prototype                           

pme = tanhmap_ab_to_R(PropsSI("P", "T",  evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
ϵme = ϵmsup =  tanhmap_01_to_R(0.95) 
xms = tanhmap_01_to_R(0.95)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Ts_l, Ts_m = Ts_h - 4, Ts_h - 2

u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m]

# prob = NonlinearProblem(evaporator_model!,u0, param_evaporator_prototype, abstol = 1e-7)  
prob = NonlinearProblem(evaporator_model!, sol.u, param_evaporator_prototype, abstol = 1e-7)  
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporator_solution = evaporator_model(sol.u, param_evaporator_prototype)

x, Tc, Th = temperatureprofiles(evaporator_solution; n = 11)

W  = compressor_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)
Q1 = condenser_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)
Q2 = evaporator_solution.twophase_region.Q + evaporator_solution.superheat_region.Q

Tmean = (evaporator_solution.superheat_region.Thh  + evaporator_solution.twophase_region.Thl) / 2
Δp_evap = pressure_drop_evap(evaporator_solution.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean)
Ppump = pumping_power(Δp_evap, evaporator_solution.superheat_region.Vh, ηpump)


p_ev2 = EvaporatorParameters(param_evaporator_prototype, Vs = 0.5*evaporator.Vs)
prob2 = NonlinearProblem(evaporator_model!,sol.u, p_ev2, abstol = 1e-7)  
sol2  = solve(prob2, TrustRegion(); abstol = 1e-9)
evaporator_solution2 = evaporator_model(sol2.u, p_ev2)
x2, Tc2, Th2 = temperatureprofiles(evaporator_solution2; n = 51)

W_2  = compressor_power(273.15 + 40, 0.7, evaporator_solution2, refrigerant)
Q1_2 = condenser_power(273.15 + 40, 0.7, evaporator_solution2, refrigerant)
Q2_2 = evaporator_solution2.twophase_region.Q + evaporator_solution2.superheat_region.Q

Tmean_2 = (evaporator_solution2.superheat_region.Thh  + evaporator_solution2.twophase_region.Thl) / 2
Δp_evap_2 = pressure_drop_evap(evaporator_solution2.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean_2)
Ppump_2 = pumping_power(Δp_evap_2, evaporator_solution2.superheat_region.Vh, ηpump)
# CopS = [round(x, digits=2) for x in [Q1/W, Q1_2/W_2] ]

p_ev3 = EvaporatorParameters(param_evaporator_prototype, Vs = 0.3*evaporator.Vs)
prob3 = NonlinearProblem(evaporator_model!,sol2.u, p_ev3, abstol = 1e-7)  
sol3  = solve(prob3, TrustRegion(); abstol = 1e-9)
evaporator_solution3 = evaporator_model(sol3.u, p_ev3)
x3, Tc3, Th3 = temperatureprofiles(evaporator_solution3; n = 51)

W_3  = compressor_power(273.15 + 40, 0.7, evaporator_solution3, refrigerant)
Q1_3 = condenser_power(273.15 + 40, 0.7, evaporator_solution3, refrigerant)
Q2_3 = evaporator_solution3.twophase_region.Q + evaporator_solution3.superheat_region.Q

Tmean_3 = (evaporator_solution3.superheat_region.Thh  + evaporator_solution3.twophase_region.Thl) / 2
Δp_evap_3 = pressure_drop_evap(evaporator_solution3.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean_3)
Ppump_3 = pumping_power(Δp_evap_3, evaporator_solution3.superheat_region.Vh, ηpump)

##.

CopS = [round(x, digits=2) for x in [Q1/(W + Ppump), Q1_2/(W_2 + Ppump_2), Q1_3/(W_3 + Ppump_3)] ]



##.

fig = Figure(; size = (600, 400))

ga = fig[1, 1] = GridLayout()
ax = Axis(fig[1, 1], xlabel = "x (%)", ylabel = "Temperature (°C)")

lref1 = lines!(ax, x, Tc .- 273.15, color = :grey, linewidth = 2, label = "refrigerant", visible = false)
ls1   = lines!(ax, x, Th .- 273.15, color = :grey, linestyle = :dash, linewidth =2, label = "secondary fluid", visible = false)

lines!(ax, x, Tc .- 273.15, color = :dodgerblue, linewidth = 2)
lines!(ax, x, Th .- 273.15, color = :dodgerblue, linestyle = :dash, linewidth =2)
vlines!(ax, [evaporator_solution.xs], color = :dodgerblue, linestyle = :dot, linewidth =1)
lines!(ax, x2, Tc2 .- 273.15, color = :chocolate1, linewidth =1)
lines!(ax, x2, Th2 .- 273.15, color = :chocolate1, linestyle = :dash, linewidth = 2)
vlines!(ax, [evaporator_solution2.xs], color = :chocolate1, linestyle = :dot, linewidth =1)
lines!(ax, x3, Tc3 .- 273.15, color = :green, linewidth = 2)
lines!(ax, x3, Th3 .- 273.15, color = :green, linestyle = :dash, linewidth =2 )
vlines!(ax, [evaporator_solution3.xs], color = :green, linestyle = :dot, linewidth = 1)
ylims!(ax,250 -273.15,290. -273.15)

axislegend(ax, [lref1, ls1], ["Tsf", "Tref"], color = [:black, :black],  position = :rb,
orientation = :horizontal)

fig

gb = fig[1, 2] = GridLayout()

# repeat([1,2,3],2)
# repeat([1,2,2],2)
# vcat([fill(i,3) for i=1:2]...)

nconf = 3

tbl = (power = [Q1, Q2, W + Ppump, Q1_2, Q2_2, W_2 + Ppump_2, Q1_3, Q2_3, W_3 + Ppump_3]./1000,
       grp_type          = repeat([1,2,3],nconf), #[1,2,3,1,2,3],
       grp_configuration = vcat([fill(i,3) for i=1:nconf]...), #[1,1,1,2,2,2],
       grp_test = repeat([1,2,2],nconf) #,       [1,2,2,1,2,2]
       )

ax2 = Axis(gb[1, 1], ylabel = "Power (kW)")
barplot!(ax2,
    tbl.grp_configuration, tbl.power,
    stack = tbl.grp_configuration,
    dodge = tbl.grp_test,
    color = tbl.grp_type,
    colormap = [:navyblue, :firebrick, :salmon],
)

ax2.xticks[] = ([1,2,3],["Vs","0.5Vs","0.3Vs"])
ylims!(ax2,0,9. )
hidexdecorations!(ax2, grid = false, label = false)


ax3 = Axis(gb[2, 1], ylabel = "COPh")
barplot!(ax3,
    1:nconf,CopS[1:nconf],
    bar_labels = string.(CopS[1:nconf]),
    flip_labels_at=0.85,
    color = [:dodgerblue, :chocolate1,:green],
)

ylims!(ax3,0,6.5 )
ax3.xticks[] = ([1,2,3],["Vs","0.5Vs","0.3Vs"])
linkxaxes!(ax2,ax3)

##
# Legend
labels = ["Q1", "Q2", "W"]
elements = [PolyElement(polycolor = c) for c in [:navyblue, :firebrick, :salmon]]
title = ""
Legend(gb[0,1], elements, labels, title,
        tellheight = true,
        tellwidth = false,
        # margin = (10, 10, 10, 10),
        halign = :center, 
        valign = :bottom, 
        orientation = :horizontal
        )

fig
##.

mypath = joinpath(@__DIR__,"figures")
save(joinpath(mypath,"variableVs_T$(round(Int,Ts_h)).png"), fig)
