####################
# evaporator model #
####################
using OptimalGshpFlow, CoolProp, Parameters, LinearAlgebra

sec_fl_name="MEA-14%" # MEA-28% MEG-37%

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::"*sec_fl_name
Q_evap=6000.
evaporator = EvaporatorDesign(Q_evap, 273.15 - 8 ,273.15+40, refrigerant_name,secondaryfluid_name; U =2000., ΔTs = 3., ΔT_l = 1.5)

T_f=CoolProp.PropsSI("T_FREEZE", secondaryfluid_name) - 273.15

Trange = range(T_f + 1. + 273.15, 40 + 273.15, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

# we have the sizing parameter for a generic evaporator (A,Vs,Vref), and we now need to choose the size L and W of the plate heat exchanger
plateheatexchanger = PlateHeatExchanger(evaporator.A, 0.5, 0.15; b=0.002, ϕ=1.189, β = 30., Dport=0.05,N_pass=1)  # NB! chevron angle vs. reverse chev. ang. !

##################
# borehole model #
##################
dout = 0.04 #internal pipe diameter
H=Q_evap/24
borehole = SingleUPipeBorehole(H = H, D = 4., rpo = dout/2, SDR=17., rb = 0.115/2, λg = 2., pipe_position = ((0.035, 0.0), (-0.035, 0.0)))
@assert norm([borehole.pipe_position[1]...]) + borehole.rpo < borehole.rb  #check that at least the pipes are inside the borehole
@assert norm([borehole.pipe_position[2]...]) + borehole.rpo < borehole.rb  #check that at least the pipes are inside the borehole

λground = 3.

Tb = 273.15 + 0
ΔTsup = 5
ηpump = 0.4



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
borehole

                                        
pme = sigmoid_ab_to_R(PropsSI("P", "T",  evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
ϵme = ϵmsup =  sigmoid_01_to_R(0.4) 
xms = sigmoid_01_to_R(0.9)
Qe, Qsup = evaporator.Qnom* 0.9, evaporator.Qnom*0.1
Ts_l, Ts_m, Ts_h = Tb - 6, Tb - 4, Tb - 2
u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Ts_l, Ts_m, Ts_h]


function solve_configuration(u0,p)
    prob = NonlinearProblem(evaporatorborehole_model!,u0,p,abstol = 1e-7) 
    sol  = solve(prob, TrustRegion(); abstol = 1e-9)
    return sol
end


params = [EvaporatorBoreholeParameters(param_evaporatorborehole_prototype,Vs = V) for V in  evaporator.Vs .* (2:-0.05:0.3)]
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
    COP1 = Q1/W
    COP2 = Q1/(W  + Ppump)
    return COP1, COP2, Q1, evaporatorborehole_solution
end


COP1 = [evaluate_solution(sol,p)[1] for (sol,p) in zip(sols,params)]
COP2 = [evaluate_solution(sol,p)[2] for (sol,p) in zip(sols,params)]
Q1 = [evaluate_solution(sol,p)[3] for (sol,p) in zip(sols,params)]
ΔTs = [evaluate_solution(sol,p)[4].evaporator.superheat_region.Thh .- evaluate_solution(sol,p)[4].evaporator.twophase_region.Thl for (sol,p) in zip(sols,params)]
Tin_bh = [evaluate_solution(sol,p)[4].evaporator.twophase_region.Thl for (sol,p) in zip(sols,params)]
Tout_bh = [evaluate_solution(sol,p)[4].evaporator.superheat_region.Thh for (sol,p) in zip(sols,params)]
Tevap = [evaluate_solution(sol,p)[4].evaporator.twophase_region.Tcl for (sol,p) in zip(sols,params)]
Vs   = [p.Vs for p in params]
T_mean=(Tin_bh+Tout_bh)/2

# Q1.*(1 .- 1 ./COP1)
Re = @. evaporator.Vs .* (2:-0.05:0.3) / (π * (dout * (1 - 2 / 17))^2 / 4) * dout * (1 - 2 / 17) * secondaryfluid.ρ.(T_mean) / secondaryfluid.μ.(T_mean)
transflow_idx = findfirst(x->x<2300,Re)

ΔTsf = Tout_bh .- Tin_bh
ΔTsf_idx = findfirst(x-> x>3, ΔTsf)

##.
fig = Figure(; size = (600, 600))
ga = fig[1, 1] = GridLayout()
ax = Axis(fig[1, 1], xlabel = "Temperature difference [K]", ylabel = "COP")
lcop1 = lines!(ax, ΔTs, COP1, color = :dodgerblue, linewidth = 2, label = "COP1")
lcop2 = lines!(ax, ΔTs, COP2, color = :indianred, linewidth = 2, label = "COP2")
WGLMakie.ylims!(ax,3,5)

# Add a secondary axis for Q1
ax2 = Axis(fig[1, 1], xlabel = "Temperature difference [K]", ylabel = "Q1 [kW]", yaxisposition = :right)
lq1 = lines!(ax2, ΔTs, Q1/1000, color = :green, linewidth = 2, label = "Q1")
WGLMakie.ylims!(ax2,6,10)
fig

# Add dashed vertical line for the maximum of COP2
max_COP2_index = argmax(COP2)
lmax =vlines!(ax, ΔTs[max_COP2_index], linestyle = :dash, color = :black)

# ltransflow = vlines!(ax,ΔTs[transflow_idx-1:transflow_idx], linestyle = :dot, color = :gray)
lΔT3flow   = vlines!(ax,ΔTs[ΔTsf_idx], linestyle = :dash, color = :red)

# axislegend(ax, [lcop1, lcop2, lmax,ltransflow,lΔT3flow ], ["COP1", "COP2", "max(COP2)", "lam/turb", "ΔT=3°C control"], position = :rt)
axislegend(ax, [lcop1, lcop2, lq1], ["COP1", "COP2", "Q1"], position = :rt) # 

# axislegend(ax, [sc1, sc2], ["One", "Two"], "Selected Dots", position = :rb,
#     orientation = :horizontal)
# axislegend(ax, [lmax,ltransflow,lΔT3flow ], ["max(COP2)", "lam/turb", "ΔT=3°C control"], position = :rt)
# axislegend(ax2, [lq1], ["Q1"], position = :rb,)
hidespines!(ax2, :l, :b, :t)
hidexdecorations!(ax2,grid = true);
hideydecorations!(ax2,grid = true, label = false, ticklabels = false, ticks = false,);



# Save the first plot
mypath = joinpath(@__DIR__, "figures")
# save(joinpath(mypath, "optimalCOP_Tb_$(round(Int, Tb))_deltaTsup$(round(Int, div(ΔTsup, 1)))_$(round(Int, mod(ΔTsup, 1) * 10)).png"), fig)

# Add a subplot for T_mean vs flow rate under the existing plot
ax3 = Axis(fig[2, 1], xlabel = "Temperature difference [K]", ylabel = "Temperature [°C]")
ltin= lines!(ax3, ΔTs, Tin_bh .- 273.15, color = :blue, linewidth = 2, label = "T_in_bh")
ltout = lines!(ax3, ΔTs, Tout_bh .- 273.15, color = :red, linewidth = 2, label = "T_out_bh")
ltb = lines!(ax3, ΔTs, Tb*ones(length(Vs)) .- 273.15, color = :green, linewidth = 2, label = "T_out_bh")
lte = lines!(ax3, ΔTs, Tevap .- 273.15, color = :darkblue, linewidth = 2, label = "T_evap")
WGLMakie.ylims!(ax3,-15,11)

# Add a subplot for Reynolds number vs flow rate under the existing plots
# ax4 = Axis(fig[3, 1], xlabel = "Flow rate [l/s]", ylabel = "Reynolds Number")
# lre = lines!(ax4, Vs * 1000, Re, color = :purple, linewidth = 2, label = "Re")

# axislegend(ax4, [lre], ["Re"], position = :rb)

axislegend(ax3, [ltb, ltin, ltout, lte], ["Borehole wall", "Inlet borehole","Outlet borehole","Evaporation"], position = :rb)

# # Add a subplot for Reynolds number vs flow rate under the existing plots
# ax4 = Axis(fig[3, 1], xlabel = "Flow rate [l/s]", ylabel = "Reynolds Number")
# lre = lines!(ax4, Vs * 1000, Re, color = :purple, linewidth = 2, label = "Re")

# # add dotted lined to highlight laminar turbulent transition
# vlines!(ax4,Vs[transflow_idx-1:transflow_idx] * 1000, linestyle = :dot, color = :gray)
# hlines!(ax4,[2300.], linestyle = :dot, color = :gray)


##.
# Add dashed vertical line for the maximum of T_mean
# max_T_mean_index = argmax(T_mean)
# lmax_tmean = lines!(ax3, ones(2).*Vs[max_T_mean_index]*1000, [minimum(T_mean), maximum(T_mean)], linestyle = :dash, color = :black)


# Save the combined plot
WGLMakie.xlims!(ax,2,10)
WGLMakie.xlims!(ax2,2,10)
str_save="opt_combined_Tb$(round(Int, Tb - 273.15))°C_sup$(round(ΔTsup, digits=1))K_Q_evap$(round(Q_evap/1000;digits=1))kW_H$(round(H;digits=0))m_sf-$(sec_fl_name)"
save(joinpath(mypath, str_save*".png"), fig)
linkxaxes!(ax,ax2,ax3)
# WGLMakie.ylims!(ax,2,5.5)
# WGLMakie.ylims!(ax2,4,10)

WGLMakie.DataInspector(fig)
# fig


### Same plot in interactive form with Plotly
Plots.plotly()
# Create the plot with specified size
p = Plots.plot(ΔTs, COP1, label = "COP1", color = :dodgerblue, linewidth = 2, ylabel = "COP", legend = :topright, size = (600, 600))
Plots.plot!(ΔTs, COP2, label = "COP2", color = :indianred, linewidth = 2)
# Plots.plot!(twinx(p), Q1, label = "Q1", color = :green, linewidth = 2, ylabel = "Q1",axis=:right)

# Add dashed vertical line for the maximum of COP2
max_COP2_index = argmax(COP2)
Plots.plot!(ones(2).*ΔTs[max_COP2_index], [minimum(COP2), maximum(COP2)],label = "max", linestyle = :dash, color = :black)

# Add a subplot for T_mean vs temperature difference
p2 = Plots.plot(ΔTs, Tin_bh .- 273.15, label = "T_in_bh", color = :blue, linewidth = 2, xlabel = "Temperature difference (K)", ylabel = "Temperature [°C]", size = (600, 600))
Plots.plot!(ΔTs, Tout_bh .- 273.15, label = "T_out_bh", color = :red, linewidth = 2, xlabel = "Temperature difference (K)", ylabel = "Temperature [°C]")
Plots.plot!(ΔTs, Tb*ones(length(Vs)) .- 273.15, label = "T_b", color = :green, linewidth = 2, xlabel = "Temperature difference (K)", ylabel = "Temperature [°C]")
# Combine the plots into a single figure
p_combined = Plots.plot(p, p2, layout = (2, 1))

# Save the plot
mypath = joinpath(@__DIR__, "figures")
Plots.savefig(p, joinpath(mypath, str_save*".html"))

# # Display the plot
# p_combined
