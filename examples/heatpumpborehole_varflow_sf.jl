####################
# evaporator model #
####################

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

Tb = 273.15 + 10
ΔTsup = 5.
ηpump = 0.2

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
ϵme = ϵmsup =  sigmoid_01_to_R(0.5) 
xms = sigmoid_01_to_R(0.5)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Ts_l, Ts_m = Ts_h - 4, Ts_h - 2

u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m, Ts_h]

prob = NonlinearProblem(evaporatorborehole_model!,u0,param_evaporatorborehole_prototype,abstol = 1e-7) 
sol  = solve(prob, TrustRegion(); abstol = 1e-9)
evaporatorborehole_solution = evaporatorborehole_model(sol.u,param_evaporatorborehole_prototype)


x, Tc, Th = temperatureprofiles(evaporatorborehole_solution.evaporator; n = 11)

W  = compressor_power(273.15 + 40, 0.7, evaporatorborehole_solution.evaporator, refrigerant)
Q1 = condenser_power(273.15 + 40, 0.7, evaporatorborehole_solution.evaporator, refrigerant)
Q2 = evaporatorborehole_solution.evaporator.twophase_region.Q + evaporatorborehole_solution.evaporator.superheat_region.Q

# compute pumping power in the borehole
@unpack rpi, H = evaporatorborehole_solution.borehole.borehole
@unpack Tin, Tout = evaporatorborehole_solution.borehole
@unpack Vh = evaporatorborehole_solution.evaporator.superheat_region
Ppump_borehole = pumping_power(ηpump, 2H, Vh ,2rpi, secondaryfluid, (Tin + Tout)/2 )
#compute pumping power in the evaporator
Tmean = (evaporatorborehole_solution.evaporator.superheat_region.Thh  + evaporatorborehole_solution.evaporator.twophase_region.Thl) / 2
Δp_evap = pressure_drop_evap(evaporatorborehole_solution.evaporator.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean)
Ppump_evap = pumping_power(Δp_evap, evaporatorborehole_solution.evaporator.superheat_region.Vh, 0.4)
#compute total pumping Power
Ppump = Ppump_borehole + Ppump_evap

## example 2

param2 = EvaporatorBoreholeParameters(param_evaporatorborehole_prototype, Vs = 0.5evaporator.Vs)
prob2 = NonlinearProblem(evaporatorborehole_model!,u0,param2,abstol = 1e-7) 
sol2  = solve(prob2, TrustRegion(); abstol = 1e-9)
evaporatorborehole_solution2 = evaporatorborehole_model(sol2.u,param2)


x2, Tc2, Th2 = temperatureprofiles(evaporatorborehole_solution2.evaporator; n = 11)

W_2  = compressor_power(273.15 + 40, 0.7, evaporatorborehole_solution2.evaporator, refrigerant)
Q1_2 = condenser_power(273.15 + 40, 0.7, evaporatorborehole_solution2.evaporator, refrigerant)
Q2_2 = evaporatorborehole_solution2.evaporator.twophase_region.Q + evaporatorborehole_solution2.evaporator.superheat_region.Q

# compute pumping power in the borehole
@unpack rpi, H = evaporatorborehole_solution2.borehole.borehole
@unpack Tin, Tout = evaporatorborehole_solution2.borehole
@unpack Vh = evaporatorborehole_solution2.evaporator.superheat_region
Ppump_borehole_2 = pumping_power(ηpump, 2H, Vh ,2rpi, secondaryfluid, (Tin + Tout)/2 )
#compute pumping power in the evaporator
Tmean_2 = (evaporatorborehole_solution2.evaporator.superheat_region.Thh  + evaporatorborehole_solution2.evaporator.twophase_region.Thl) / 2
Δp_evap_2 = pressure_drop_evap(evaporatorborehole_solution2.evaporator.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean_2)
Ppump_evap_2 = pumping_power(Δp_evap_2, evaporatorborehole_solution2.evaporator.superheat_region.Vh, 0.4)
#compute total pumping Power
Ppump_2 = Ppump_borehole_2 + Ppump_evap_2

## example 3

param3 = EvaporatorBoreholeParameters(param_evaporatorborehole_prototype, Vs = 0.3evaporator.Vs)
prob3 = NonlinearProblem(evaporatorborehole_model!,sol2.u,param3,abstol = 1e-7) 
sol3  = solve(prob3, TrustRegion(); abstol = 1e-9)
evaporatorborehole_solution3 = evaporatorborehole_model(sol3.u,param3)


x3, Tc3, Th3 = temperatureprofiles(evaporatorborehole_solution3.evaporator; n = 11)

W_3  = compressor_power(273.15 + 40, 0.7, evaporatorborehole_solution3.evaporator, refrigerant)
Q1_3 = condenser_power(273.15 + 40, 0.7, evaporatorborehole_solution3.evaporator, refrigerant)
Q2_3 = evaporatorborehole_solution3.evaporator.twophase_region.Q + evaporatorborehole_solution3.evaporator.superheat_region.Q

# compute pumping power in the borehole
@unpack rpi, H = evaporatorborehole_solution3.borehole.borehole
@unpack Tin, Tout = evaporatorborehole_solution3.borehole
@unpack Vh = evaporatorborehole_solution3.evaporator.superheat_region
Ppump_borehole_3 = pumping_power(ηpump, 2H, Vh ,2rpi, secondaryfluid, (Tin + Tout)/2 )
#compute pumping power in the evaporator
Tmean_3 = (evaporatorborehole_solution3.evaporator.superheat_region.Thh  + evaporatorborehole_solution3.evaporator.twophase_region.Thl) / 2
Δp_evap_3 = pressure_drop_evap(evaporatorborehole_solution3.evaporator.superheat_region.Vh, plateheatexchanger, secondaryfluid,Tmean_3)
Ppump_evap_3 = pumping_power(Δp_evap_3, evaporatorborehole_solution3.evaporator.superheat_region.Vh, 0.4)
#compute total pumping Power
Ppump_3 = Ppump_borehole_3 + Ppump_evap_3

@assert sol.retcode == ReturnCode.Success
@assert sol2.retcode == ReturnCode.Success
@assert sol3.retcode == ReturnCode.Success

CopS = [round(x, digits=2) for x in [Q1/(W + Ppump), Q1_2/(W_2 + Ppump_2), Q1_3/(W_3 + Ppump_3)] ]


##.

fig = Figure(; size = (600, 400))

ga = fig[1, 1] = GridLayout()
ax = Axis(fig[1, 1], xlabel = "x (%)", ylabel = "Temperature (°C)")
lines!(ax, [0.,1.], (evaporatorborehole_solution.borehole.Tb -273.15) .* [1,1], color = :brown, linewidth = 2)

lref1 = lines!(ax, x, Tc .- 273.15, color = :grey, linewidth = 2, label = "refrigerant", visible = false)
ls1   = lines!(ax, x, Th .- 273.15, color = :grey, linestyle = :dash, linewidth =2, label = "secondary fluid", visible = false)

lines!(ax, x, Tc .- 273.15, color = :dodgerblue, linewidth = 2)
lines!(ax, x, Th .- 273.15, color = :dodgerblue, linestyle = :dash, linewidth =2)
vlines!(ax, [evaporatorborehole_solution.evaporator.xs], color = :dodgerblue, linestyle = :dot, linewidth =1)

lines!(ax, x2, Tc2 .- 273.15, color = :chocolate1, linewidth =2)
lines!(ax, x2, Th2 .- 273.15, color = :chocolate1, linestyle = :dash, linewidth = 2)
vlines!(ax, [evaporatorborehole_solution2.evaporator.xs], color = :chocolate1, linestyle = :dot, linewidth =1)

lines!(ax, x3, Tc3 .- 273.15, color = :green, linewidth = 2)
lines!(ax, x3, Th3 .- 273.15, color = :green, linestyle = :dash, linewidth =2)
vlines!(ax, [evaporatorborehole_solution3.evaporator.xs], color = :green, linestyle = :dot, linewidth = 1)
ylims!(ax,250 .- 273.15, 290 .- 273.15)

text!(ax, [0.1], [evaporatorborehole_solution.borehole.Tb .- 273.15 + 0.2 ],
    # text= rich(rich("Tb", font="Bold"), " Tb", font =:regular),
    text= rich("borehole temperature", font =:regular),
    align = (:left, :bottom), color = :brown)

axislegend(ax, [lref1, ls1], ["Tsf", "Tref"], color = [:black, :black],  position = :rb,
orientation = :horizontal)

# fig

gb = fig[1, 2] = GridLayout()

repeat([1,2,3],2)
repeat([1,2,2],2)
vcat([fill(i,3) for i=1:2]...)

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
# ax3.xticks[] = ([1,2,3],["Vs\nVref","0.5Vs\n0.5Vref","0.3Vs\n0.3Vref"])
ax3.xticks[] = ([1,2,3],["Vs","0.5Vs","0.3Vs"])
linkxaxes!(ax2,ax3)
fig

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
# save(joinpath(mypath,"boreholeheatpump_variableVs_T$(round(Int,Tb)).png"), fig)
save(joinpath(mypath,"boreholeheatpump_variableVs_T$(round(Int,Tb))_deltaTsup$(round(Int,div(ΔTsup,1)))_$(round(Int,mod(ΔTsup,1)*10)).png"), fig)
# fig