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
Ts_h = 273.15 + 5 ;

# p  = refrigerant, secondaryfluid, plateheatexchanger, evaporator.hei, evaporator.Vs, evaporator.Vref, ΔTsup, Tsf_h

param_evaporator_prototype  = EvaporatorParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
                                        hei = evaporator.hei, Vs = evaporator.Vs, Vref = evaporator.Vref, 
                                        ΔTsup =  ΔTsup, Ts_h = Ts_h)

# you can use the template and just change one or two parameters instead of redifining everything all the time
# p2_ev = EvaporatorParameters(param_evaporator_prototype, Vs = 2*evaporator.Vs)
                                        
# @unpack plateheatexchanger = param_evaporator_prototype                           

pme = tanhmap_ab_to_R(PropsSI("P", "T",  evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
ϵme = ϵmsup =  tanhmap_01_to_R(0.7) 
xms = tanhmap_01_to_R(0.7)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Ts_l, Ts_m = Ts_h - 4, Ts_h - 2

u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m]

# prob = NonlinearProblem(evaporator_model!,u0, param_evaporator_prototype, abstol = 1e-7)  
prob = NonlinearProblem(evaporator_model!, u0, param_evaporator_prototype, abstol = 1e-7)  
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporator_solution = evaporator_model(sol.u, param_evaporator_prototype)

x, Tc, Th = temperatureprofiles(evaporator_solution; n = 11)

W  = compressor_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)
Q1 = condenser_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)
Q2 = evaporator_solution.twophase_region.Q + evaporator_solution.superheat_region.Q

Tmean = (evaporator_solution.superheat_region.Thh  + evaporator_solution.twophase_region.Thl) / 2
Δp_evap = pressure_drop_evap(evaporator_solution.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean)
Ppump = pumping_power(Δp_evap, evaporator_solution.superheat_region.Vh, 0.4)


##.

CopS = [round(x, digits=2) for x in [Q1/(W + Ppump)] ]



##.

fig = Figure(; size = (600, 400))

ga = fig[1, 1] = GridLayout()
ax = Axis(fig[1, 1], xlabel = "x (%)", ylabel = "Temperature (K)")
lines!(ax, x, Tc, color = :dodgerblue, linewidth = 2)
lines!(ax, x, Th, color = :dodgerblue, linestyle = :dash, linewidth =2)
vlines!(ax, [evaporator_solution.xs], color = :dodgerblue, linestyle = :dot, linewidth =1)
# lines!(ax, x2, Tc2, color = :chocolate1, linewidth =1)
# lines!(ax, x2, Th2, color = :chocolate1, linestyle = :dash, linewidth = 2)
# vlines!(ax, [evaporator_solution2.xs], color = :chocolate1, linestyle = :dot, linewidth =1)
# lines!(ax, x3, Tc3, color = :green, linewidth = 2)
# lines!(ax, x3, Th3, color = :green, linestyle = :dash, linewidth =2 )
# vlines!(ax, [evaporator_solution3.xs], color = :green, linestyle = :dot, linewidth = 1)
# ylims!(ax,272,280.)
ylims!(ax,250,290.)


fig

gb = fig[1, 2] = GridLayout()

# repeat([1,2,3],2)
# repeat([1,2,2],2)
# vcat([fill(i,3) for i=1:2]...)

nconf = 1

tbl = (power = [Q1, Q2, W + Ppump]./1000,
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
    color = [:dodgerblue],
)

ylims!(ax3,0,6.5 )
ax3.xticks[] = ([1],["Vs"])
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

##.

mypath = joinpath(@__DIR__,"figures")
save(joinpath(mypath,"heatpump_singleVs.png"), fig)
