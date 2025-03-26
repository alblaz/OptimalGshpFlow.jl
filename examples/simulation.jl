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
        Q_evap=6000.;
        L=0.5;
        W=0.15;
        T_evap_design= 0. + K
        # Tsf_h = 278.;
        ################################
        # borehole model #
        ################################

        H = Q_evap/30. #30 W/m
        Dpo=H > 250 ? 0.05 : 0.04
        SDR=17.
        Dpi=Dpo*(1-2/SDR)
        D = 4.
        λg = 3. #grout
        pipe_position = ((0.04, 0.0), (-0.04, 0.0))
        λground = 3.5 #ground
        Tin, Tb  = 1. + K , 5. + K

        ################################
        # compressor #
        ################################


        ################################
        # condenser #
        ################################
        T_cond=40. + K

        ################################
        # simulation parameters #
        ################################
        reltol=1e-3

########################################
# simulation #
########################################

################################
## preparation phase ##
################################

refrigerant_name, secondaryfluid_name = "R290", "INCOMP::MEA-28%"
evaporator = EvaporatorDesign(Q_evap, T_evap_design, T_cond, refrigerant_name,secondaryfluid_name;
        U =2000., ΔTs = 3., ΔT_l = 1.5, ΔT_sub=ΔTsub, ΔT_sup=ΔTsup, ηᵥ=ηᵥ)

Trange = range(-15. + K, 30 + K, length = 40) #temperature range used to build secondary fluid interpolants
secondaryfluid = SecondaryFluid(Trange, secondaryfluid_name)
refrigerant = Refrigerant(refrigerant_name)

# we have the sizing parameter for a generic evaporator (A,Vs,Vref), and we now need to choose the size L and W of the plate heat exchanger
# PlateHeatExchanger(A,L,W; b=0.002, ϕ=1.189, β = 30.) --> hence giving th enumber of plates as output
@unpack Vref,Vs,hei,ΔT_sup,Thl,Thh,A=evaporator
plateheatexchanger = PlateHeatExchanger(A, L, W; b=0.002, ϕ=1.189, β = 30.)


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
Tsf_l, Tsf_m = evaporator.Thl,  evaporator.Thl+(evaporator.Thh -  evaporator.Thl)*3/4

u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m]

prob = NonlinearProblem(evaporator_model!,u0,p,reltol = reltol)  
sol_evap  = solve(prob, TrustRegion(); reltol = reltol)
evaporator_solution = evaporator_model(sol_evap.u,p)


##################
# borehole model #
##################
borehole = SingleUPipeBorehole(H=H, D=D, λg=λg, pipe_position = pipe_position, rpo=Dpo/2, SDR=SDR)
u0 = [1. + K]
p = BoreholeParameters(secondaryfluid=secondaryfluid, Vs=Vs, Tin=Tin, borehole=borehole, λground=λground, Tb=Tb)
prob = NonlinearProblem(borehole_model!,u0,p,reltol = reltol)  
sol_bh  = solve(prob, TrustRegion(); reltol = reltol)


#############################
# evaporator-borehole model #
#############################
# The model evaporatorborehole_model! combine the models above by imposing Tin = Ts_l and Tout = Ts_h

p = EvaporatorBoreholeParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
                                        hei = hei, Vs = Vs, Vref = Vref, 
                                        ΔTsup =  ΔT_sup, 
                                        borehole=borehole, λground=λground, Tb=Tb)
u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, evaporator.Thh]
prob = NonlinearProblem(evaporatorborehole_model!,u0,p,reltol = reltol)
try 
        sol_evapbh  = solve(prob, TrustRegion(); reltol = reltol)
catch e
end

a=true
i=0
max_iter=20
call_i(;i=i,max_iter=max_iter)=i,max_iter
while a || (sol_evapbh.retcode!=ReturnCode.Success && i<=max_iter)
        i,max_iter=call_i()
        @show i
        i += 1
        a=false
        try
            prob = NonlinearProblem(evaporatorborehole_model!, u0, p, reltol=reltol)
            sol_evapbh = solve(prob, TrustRegion(); reltol=reltol)
        catch e
            # @show tanhmap_R_to_01(u0[4]),tanhmap_R_to_01(u0[7])
        end
        u0[4] = u0[5] = tanhmap_01_to_R(min(max(rand(1)[1],0.03),0.97))
        u0[7]=tanhmap_01_to_R(min(max(rand(1)[1],0.03),0.97))
        if i == max_iter
            @show sol_evapbh.retcode
            println("/!\\ not converged")
        elseif sol_evapbh.retcode == ReturnCode.Success
            println(":) converged") 
        end
end
evaporatorborehole_solution = evaporatorborehole_model(sol_evapbh.u,p)

#############################
# Plotting #
#############################



x, Tc, Th = temperatureprofiles(evaporator_solution; n = 15)

W  = compressor_power(T_cond, 0.7, evaporator_solution, refrigerant)
Q1 = condenser_power(T_cond, 0.7, evaporator_solution, refrigerant)

COP₁=Q1/W

using Plots
plotly()
plot(x, Th .- K, label="secondary fluid")
plot!(x, Tc .- K, label="refrigerant")

# compute pumping power in the broehole
@unpack rpi, H = evaporatorborehole_solution.borehole.borehole
@unpack Tin, Tout = evaporatorborehole_solution.borehole
@unpack Vc = evaporatorborehole_solution.evaporator.superheat_region
Ppump = pumping_power(0.4, H, Vc ,2rpi, secondaryfluid, (Tin + Tout)/2,evaporatorborehole_solution.evaporator.plateheatexchanger)

COP₂=Q1/(W+Ppump)

#############################
# Range of flow rates #
#############################

# function to pass initial value in local scope
function pass_u0(Vs,Q_evap,ρ_cₚ,to_pass="sol_evapbh";sol_evap=sol_evap, sol_bh=sol_bh, sol_evapbh=sol_evapbh,u0_def=[pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, evaporator.Thh])
        if to_pass == "sol_evapbh"
                u0=sol_evapbh.u
        elseif to_pass == "sol_bh"
                u0=sol_bh.u
        elseif to_pass == "sol_evap"
                u0=sol_evap.u
        else
                u0=u0_def
                u0[9]=u0_def[7] + Q_evap/(Vs * ρ_cₚ)
                u0[8]=u0_def[7] + (u0[9]-u0_def[7]) * 3/4
        end
        return u0
end

T_mean=(Thl+Thh)/2
ΔT_range=range(0.5, 7, length=61)
V_range=Q_evap./(ΔT_range * secondaryfluid.ρ(T_mean) * secondaryfluid.cp(T_mean))


####################
# evaporator model #
####################

evaporator_solutions = []
borehole_solutions=[]
evapbh_solutions=[]
x,Tc,Th,W,Q₁,COP₁,Ppump,COP₂=[],[],[],[],[],[],[],[]


function simulation(evaporator_solutions, borehole_solutions, evapbh_solutions, x,Tc,Th,W,Q₁,COP₁,Ppump,COP₂,Q_evap
        ;sol_evap=sol_evap,sol_bh=sol_bh,sol_evapbh=sol_evapbh,borehole = borehole, Tin=Tin, λground=λground, Tb=Tb
        )

    for Vs in V_range
        # p = EvaporatorParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid,
        
        @show Vs
        

        p = EvaporatorBoreholeParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger,
                                         hei=hei, Vs=Vs, Vref=Vref,
                                         ΔTsup=ΔT_sup,
                                         borehole=borehole, λground=λground, Tb=Tb)

        u0 = pass_u0(Vs,Q_evap,p.secondaryfluid.ρ(Tin)*p.secondaryfluid.cp(Tin); sol_evapbh=sol_evapbh)
        a=true
        i=0
        max_iter=20
        while a || (sol_evapbh.retcode!=ReturnCode.Success && i<=max_iter)
                i += 1
                a=false
                try
                        prob = NonlinearProblem(evaporatorborehole_model!, u0, p, reltol=reltol)
                        sol_evapbh = solve(prob, TrustRegion(); reltol=reltol)
                catch e
                        # @show tanhmap_R_to_01(u0[4]),tanhmap_R_to_01(u0[7])
                end
                u0[4] = u0[5] = tanhmap_01_to_R(min(max(rand(1)[1],0.03),0.97))
                u0[7]=tanhmap_01_to_R(min(max(rand(1)[1],0.03),0.97))
                if i == max_iter
                        @show sol_evapbh.retcode
                        println("/!\\ not converged")
                elseif sol_evapbh.retcode==ReturnCode.Success
                        println(":) converged") 
                end
        end
        push!(evapbh_solutions,evaporatorborehole_model(sol_evapbh.u, p))
        # this could be done outside the for loop
        push!(x, temperatureprofiles(evapbh_solutions[end].evaporator; n = 15)[1])
        push!(Tc, temperatureprofiles(evapbh_solutions[end].evaporator; n = 15)[2])
        push!(Th, temperatureprofiles(evapbh_solutions[end].evaporator; n = 15)[3])
        push!(W, compressor_power(T_cond, 0.7, evapbh_solutions[end].evaporator, refrigerant))
        push!(Q₁, condenser_power(T_cond, 0.7, evapbh_solutions[end].evaporator, refrigerant))
        push!(COP₁, Q₁[end]/ W[end])


        # compute pumping power in the borehole
        @unpack rpi, H = evapbh_solutions[end].borehole.borehole
        @unpack Tin, Tout = evapbh_solutions[end].borehole
        push!(Ppump, pumping_power(0.5, H, Vs, 2rpi, secondaryfluid, (Tin + Tout) / 2, evapbh_solutions[end].evaporator.plateheatexchanger))

        push!(COP₂, Q₁[end] / (W[end] + Ppump[end]))
        # @show Q1 / (W + Ppump)
    end
end
simulation(evaporator_solutions, borehole_solutions, evapbh_solutions, x,Tc,Th,W,Q₁,COP₁,Ppump,COP₂,Q_evap)


plot(ΔT_range, COP₁, label="COP₁", xlabel="Temperature difference [K]", ylabel="COP")
plot!(ΔT_range, COP₂, label="COP₂")



##### Trying to solve with broadcasting


# p = EvaporatorParameters.(refrigerant=refrigerant, secondaryfluid=secondaryfluid,
#                               plateheatexchanger=plateheatexchanger, 
#                               hei=hei, Vs=Vs, Vref=Vref, 
#                               ΔTsup=ΔT_sup, Ts_h=Thh)
#     ## could we just pass an evaporator type here?

#     # initial guesses for the u vector
#     pme = tanhmap_ab_to_R(PropsSI("P", "T", evaporator.Te, "Q", 1., refrigerant.name), refrigerant.pr_min, refrigerant.pcrit)
#     ϵme = ϵmsup = tanhmap_01_to_R(0.5)
#     xms = tanhmap_01_to_R(0.5)
#     Qe, Qsup = evaporator.Qnom * 0.9, evaporator.Qnom * 0.1
#     Tsf_l, Tsf_m = evaporator.Thh - 4, evaporator.Thh - 2

#     u0 = sol_evap.u
#     # [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m]

#     prob = NonlinearProblem.(evaporator_model!, u0, p, reltol=1e-7)
#     sol_evap = solve.(prob, TrustRegion(); reltol=1e-9)
#     evaporator_solution = evaporator_model(sol_evap.u,p)

#     ##################
#     # borehole model #
#     ##################
#     borehole = SingleUPipeBorehole(H=H, D=D, λg=λg, pipe_position=pipe_position, rpo=Dpo/2, SDR=SDR)
#     u0 = [2. + K]
#     p = BoreholeParameters(secondaryfluid=secondaryfluid, Vs=Vs, Tin=Tin, borehole=borehole, λground=λground, Tb=Tb)
#     prob = NonlinearProblem(borehole_model!, u0, p, reltol=1e-9)
#     sol = solve(prob, TrustRegion(); reltol=1e-9)

#     #############################
#     # evaporator-borehole model #
#     #############################
#     # The model evaporatorborehole_model! combine the models above by imposing Tin = Ts_l and Tout = Ts_h

#     p = EvaporatorBoreholeParameters(refrigerant=refrigerant, secondaryfluid=secondaryfluid, plateheatexchanger=plateheatexchanger, 
#                                      hei=hei, Vs=Vs, Vref=Vref, 
#                                      ΔTsup=ΔT_sup, 
#                                      borehole=borehole, λground=λground, Tb=Tb)
#     u0 = [pme, Qe, Qsup, ϵme, ϵmsup, xms, Tsf_l, Tsf_m, evaporator.Thh]
#     prob = NonlinearProblem(evaporatorborehole_model!, u0, p, reltol=1e-7)
#     sol = solve(prob, TrustRegion(); reltol=1e-9)
#     evaporatorborehole_solution = evaporatorborehole_model(sol.u, p)
