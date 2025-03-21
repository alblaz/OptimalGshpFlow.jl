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

W  = compressor_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)
Q1 = condenser_power(273.15 + 40, 0.7, evaporator_solution, refrigerant)



#########################
Q2,Te,Tc = 8000., K - 5, 40 + K
# refrigerant = "R134a"
refrigerant = "R290"
# refrigerant = "R134a"
A,Vref,Vsf = evaporator_design(Q2, Te, Tc, refrigerant, "INCOMP::MEA-28%")

hco = PropsSI("H", "T",  Tc, "Q", 0, refrigerant) #enthalpy at the outlet of the condenser - nominal conditions
hei = hco         

# some more relevant properties 
pcrit = PropsSI.("Pcrit",  refrigerant) # critical pressure 
M = PropsSI.("M",  refrigerant) #molar mass
pr_min = 1e-4  #400. / PropsSI.("Pcrit",  refrigerant)

# build spline relevant splines for secondary fluid properties
Trange = range(-15. +  K, 30 + K, length = 40)
secondary_fluid =  "INCOMP::MEA-28%"
spl_sf = interpolate_sf_properties(Trange, secondary_fluid)

# build spline relevant splines for refrigerant properties
spl_saturation  = interpolate_ref_properties_twophaseregion(refrigerant, 1.; ns= 500)
spl_superheated = interpolate_ref_properties_superheatedregion(refrigerant)


# plate heat exchanger parameters
ϕ = 1.189
b,L,W = 0.002, 0.2, 0.06
np = round(A/(L*W)/ϕ) + 2        # number of plates
H = (np-1)/2 * b                 # total height of the secondary fluid channels 

b,L,W = 0.002, 0.35, 0.073
np = round(A/(L*W)) + 2        # number of plates
H = (np-1)/2 * b               # total height of the secondary fluid channels 

geometrical_param = b,L,W,H,np 
ΔTsup = 7.5


# some more relevant properties 
# pcrit = PropsSI.("Pcrit",  refrigerant) # critical pressure 
# M = PropsSI.("M",  refrigerant) #molar mass
# Tsf_h = Tsf_l + 3.
p = A, Vsf, Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup, geometrical_param
u  = inverse_mapping_p(PropsSI("P", "T",  Te, "Q", 1, refrigerant), pr_min, pcrit) , Te + 4, Q2, 3000., 0.3, Te + 6, 500., 2000.
du = zeros(length(u))
evaporatorborehole_model!(du,u,p)   
@show du
# hR = 55*(u[1]/pcrit)^(0.12)/((-log10((u[1]/pcrit)))^0.55 * M^0.5) * (u[3]/A)^0.67

ΔTsup =   5
p = A, Vsf, Vref, 273.15+0, hei, pcrit, pr_min, M, ΔTsup, geometrical_param

u0  = [inverse_mapping_p(PropsSI("P", "T",   273.15-5., "Q", 1, refrigerant), pr_min, pcrit),  # mapped evaporation pressure
        273.15+0, # Tsf_l
        4000,      # Q2_e
        1000.,     # Ue
        0.,       # y  mapping for x
        273.15+3 , # Tsf_m
        100. ,    # Q2_sup
        1000.      # U_sup 
        ]



prob = NonlinearProblem(f,u0,p,abstol = 1e-7)
@time sol = solve(prob, LevenbergMarquardt())
f(du,sol.u,p)   
@show du
@show u

# prob = NonlinearProblem(ff,u0,p,abstol = 1e-10)
# @time sol = solve(prob)
# ff(du,sol.u,p)   
# @show du
# @show u


# u0  = [ (-10.,10.),
#         (230.,300.),
#         (0.,20000.),      # Q2_e
#         (0,20000.),     # Ue
#         (-10.,10.),     # y  mapping for x
#         (230.,300.),    # Tsf_m
#         (0.,20000.),   # Q2_sup
#         (0,20000.),      # U_sup 
#         ]

# prob = IntervalNonlinearProblem(f,u0,p,abstol = 1e-10)
# @time sol = solve(prob)

# map_y2x(y) = (tanh(y) + 1)/2 
# ϵ

# prob = NonlinearProblem(ff,u0,p,abstol = 1e-10)
# @time sol = solve(prob)
# ff(du,sol.u,p)   
# @show du
# @show sol.u



# map_y2x(y) = (tanh(y) + 1)/2 
# ϵ
##.
xx,Th,Tc = postprocessing_counterflowHX(sol.u,p ;xe = 0:0.02:1, xsup = 0:0.02:1)


using WGLMakie

fig = Figure(size = (400, 400))

ax = Axis(fig[1,1],xlabel = "lenght (%)", ylabel ="T (K)", title = "superheat = $ΔTsup K")
ylims!(260,275)

lines!(ax, xx, Th, color = :dodgerblue, linestyle = :solid, linewidth = 1.)
lines!(ax, xx, Tc, color = :indianred, linestyle = :solid, linewidth = 1.)

Te = spl_saturation.T(log10(mapping_p(sol.u[1],pr_min,pcrit)))


# compute (Tch - Te) - ΔTsup This should be 0 be cause it is given as input.
# There some error here, which is probably introduced by the epsilon-NTU method.
Δϵ = Tc[end] - spl_saturation.T(log10(mapping_p(sol.u[1],pr_min,pcrit))) - ΔTsup 

@show Δϵ
@show Te

fig
##.

Ten = Te

u0 = sol.u
##.

ΔTsup_a =   2.5
p1 = A, Vsf*1.2, 0.5Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup_a, geometrical_param
prob1 = NonlinearProblem(f,u0,p1,abstol = 1e-9)
@time sol1 = solve(prob1)
xx1,Th1,Tc1 = postprocessing_counterflowHX(sol1.u,p1 ;xe = 0:0.02:1, xsup = 0:0.02:1)

p2 = A, 0.3Vsf, Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup_a, geometrical_param
prob2 = NonlinearProblem(f,u0,p2,abstol = 1e-10)
@time sol2 = solve(prob2)
xx2,Th2,Tc2 = postprocessing_counterflowHX(sol2.u,p2 ;xe = 0:0.02:1, xsup = 0:0.02:1)

ΔTsup_b =   4
p3 = A, Vsf, Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup_b, geometrical_param
prob3 = NonlinearProblem(f,u0,p3,abstol = 1e-10)
@time sol3 = solve(prob3)
xx3,Th3,Tc3 = postprocessing_counterflowHX(sol3.u,p3 ;xe = 0:0.02:1, xsup = 0:0.02:1)

p4 = A, 0.3Vsf, Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup_b, geometrical_param
prob4 = NonlinearProblem(f,u0,p4,abstol = 1e-10)
@time sol4 = solve(prob4)
xx4,Th4,Tc4 = postprocessing_counterflowHX(sol4.u,p4 ;xe = 0:0.02:1, xsup = 0:0.02:1)

##.

ΔTsup_c =   5
p5 = A, Vsf, Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup_c, geometrical_param
prob5 = NonlinearProblem(f,u0,p5,abstol = 1e-10)
@time sol5 = solve(prob5)
xx5,Th5,Tc5 = postprocessing_counterflowHX(sol5.u,p5 ;xe = 0:0.02:1, xsup = 0:0.02:1)


p6 = A, .5Vsf, 1Vref, 273.15+5, hei, pcrit, pr_min, M, ΔTsup_c, geometrical_param
prob6 = NonlinearProblem(f,u0,p6,abstol = 1e-10)
@time sol6 = solve(prob6)
xx6,Th6,Tc6 = postprocessing_counterflowHX(sol6.u,p6 ;xe = 0:0.02:1, xsup = 0:0.02:1)




##.
using WGLMakie

fig = Figure(size = (600, 400))

ax = Axis(fig[1,1],xlabel = "lenght (%)", ylabel ="T (K)", title = "superheat = $ΔTsup_c K")

lines!(ax, xx5, Th5, color = :dodgerblue, linestyle = :solid, linewidth = 1.)
lines!(ax, xx5, Tc5, color = :dodgerblue, linestyle = :dash, linewidth = 1.)

lines!(ax, xx6, Th6, color = :indianred, linestyle = :solid, linewidth = 1.)
lines!(ax, xx6, Tc6, color = :indianred, linestyle = :dash, linewidth = 1.)

fig
##.


##.

using WGLMakie

fig = Figure(size = (300, 500))

ax = Axis(fig[1,1],xlabel = "lenght (%)", ylabel ="T (K)", title = "superheat = $ΔTsup_a K")

lines!(ax, xx1, Th1, color = :dodgerblue, linestyle = :solid, linewidth = 1.)
lines!(ax, xx1, Tc1, color = :dodgerblue, linestyle = :dash, linewidth = 1.)

lines!(ax, xx2, Th2, color = :indianred, linestyle = :solid, linewidth = 1.)
lines!(ax, xx2, Tc2, color = :indianred, linestyle = :dash, linewidth = 1.)

ax2 = Axis(fig[2,1],xlabel = "lenght (%)", ylabel ="T (K)", title = "superheat = $ΔTsup_b K")

lines!(ax2, xx3, Th3, color = :dodgerblue, linestyle = :solid, linewidth = 1.)
lines!(ax2, xx3, Tc3, color = :dodgerblue, linestyle = :dash, linewidth = 1.)

lines!(ax2, xx4, Th4, color = :indianred, linestyle = :solid, linewidth = 1.)
lines!(ax2, xx4, Tc4, color = :indianred, linestyle = :dash, linewidth = 1.)

ax3 = Axis(fig[3,1],xlabel = "lenght (%)", ylabel ="T (K)", title = "superheat = $ΔTsup_c K")

lines!(ax3, xx5, Th5, color = :dodgerblue, linestyle = :solid, linewidth = 1.)
lines!(ax3, xx5, Tc5, color = :dodgerblue, linestyle = :dash, linewidth = 1.)

lines!(ax3, xx6, Th6, color = :indianred, linestyle = :solid, linewidth = 1.)
lines!(ax3, xx6, Tc6, color = :indianred, linestyle = :dash, linewidth = 1.)


linkaxes!(ax,ax2,ax3)

fig
##.




# # evaporator section
# x  = map_y2x(sol.u[5])

# xe = 0:0.1:1
# UAe = sol.u[4]*p[1]*x   #U*A*x
# Tcle = Tche =  PropsSI("T", "P",  mapping_p(sol.u[1], pr_min,pcrit), "Q", 1, refrigerant)
# Thle = sol.u[2]
# Thhe = sol.u[6]
# Tsme = (sol.u[2] + sol.u[6])/2    #mean temperature of the secondary fluid 
# ρsfe = spl_sf.ρ(Tsm)              #evaluate density at mean temperature 
# cpsfe = spl_sf.cp(Tsm)            #evaluate cp at mean temperature 
# Che = Vsf*ρsfe*cpsfe
# Cce = Inf
# The,Tce = temperatures_counterflowHX(Che,Thle,Thhe,Cce,Tcle,Tche,UAe; x =xe)

# xxe = x.*xe

# # superheater section
# xsup = 0:0.1:1
# UAsup = sol.u[8]*p[1]*(1-x)   #U*A*x
# Tclsup = Tcle
# Tchsup = Tcle + p[9]
# Thlsup = sol.u[6]
# Thhsup = p[4]
# Tsmsup = (Thlsup + Thhsup)/2    #mean temperature of the secondary fluid 
# ρsfsup = spl_sf.ρ(Tsm)              #evaluate density at mean temperature 
# cpsfsup = spl_sf.cp(Tsm)            #evaluate cp at mean temperature 
# Chsup = Vsf*ρsf*cpsf
# Tcsup = (Tclsup + Tchsup)/2    #mean temperature of the secondary fluid 
# ρrefsup = spl_superheated.ρ(log10(mapping_p(sol.u[1], pr_min,pcrit)), p[9])              #evaluate density at mean temperature 
# cprefsup = (spl_superheated.h(log10(mapping_p(sol.u[1], pr_min,pcrit)), p[9]) -  spl_superheated.h(log10(mapping_p(sol.u[1], pr_min,pcrit)), 0))/p[9]            #evaluate cp at mean temperature 
# Chsup = Vsf*ρsf*cpsf
# Ccsup = Vref*ρrefsup*cprefsup
# Thsup,Tcsup = temperatures_counterflowHX(Chsup,Thlsup,Thhsup,Ccsup,Tclsup,Tchsup,UAsup; x =xsup)

# xxsup = x .+  (1-x) .* xsup

# xx = vcat(xxe, xxsup[2:end]) 
# Th = vcat(The, Thsup[2:end]) 
# Tc = vcat(Tce, Tcsup[2:end]) 

#


# ##.
# # superheater section
#     pme, Tsf_l, Q2e, Ue, y, Tsf_m,Q2sup, Usup= sol.u

#     # p is a tuple of parameters. 
#     A, Vsf, Vref, Tsf_h, hei, pcrit, pr_min, M, ΔTsup, geometrical_param = p
    
#     pe  = mapping_p(pme, pr_min, pcrit)
#     x  = map_y2x(y)
    
#     UAsup = Usup*A*(1-x)   #U*A*x
#     Tclsup = Tcle
#     Tchsup = Tcle + ΔTsup
#     Thlsup = Tsf_m
#     Thhsup = Tsf_h
#     Tsmsup = (Thlsup + Thhsup)/2    #mean temperature of the secondary fluid 
#     ρsfsup = spl_sf.ρ(Tsmsup)              #evaluate density at mean temperature 
#     cpsfsup = spl_sf.cp(Tsmsup)            #evaluate cp at mean temperature 
#     Chsup = Vsf*ρsf*cpsf
#     # Tcsup = (Tclsup + Tchsup)/2                        #mean temperature of the secondary fluid 
#     ρrefsup = spl_superheated.ρ(log10(pe), ΔTsup)              #evaluate density at mean temperature 
#     cprefsup = (spl_superheated.h(log10(pe), ΔTsup) -  spl_superheated.h(log10(pe), 0))/ΔTsup           
#     Ccsup = Vref*ρrefsup*cprefsup
#     Thsup,Tcsup = temperatures_counterflowHX(Chsup,Thlsup,Thhsup,Ccsup,Tclsup,Tchsup,UAsup; x =xsup)

#     Thsup
#     Tcsup


# ##..


