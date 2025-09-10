### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 6f32ff1e-5faa-4551-a898-01167021e4a4
using AlgebraOfGraphics: VLines, color

# ╔═╡ f5ba6b5c-9f00-43a2-ae87-4d86444c4b3f

using Colors: RGBA, red, green, blue

# ╔═╡ 9539a581-42e3-43dc-996d-4d9316aa9c92
using AlgebraOfGraphics: filled_contours, contours

# ╔═╡ 8527ed5d-3ff6-435d-a93d-548c824dee78
using Interpolations

# ╔═╡ eadd49a1-7bca-40a2-abdb-064a0fb710e6
using Roots

# ╔═╡ d47bf527-6960-4e14-b9e5-a422f8888e66
using PlutoDevMacros

# ╔═╡ 910c36c9-4ec1-4b0f-b044-823a229fa536
using PlutoTest

# ╔═╡ 00bcb75f-8da5-46a8-9c20-6b94861342c4
using SpecialFunctions

# ╔═╡ 1499b12a-8a87-4354-b4b5-d181ffd63f1b
# ╠═╡ skip_as_script = true
#=╠═╡
using DataFrames
  ╠═╡ =#

# ╔═╡ fa83b49d-d4ad-470a-88f8-524ae0e3191f
using ModelingToolkit

# ╔═╡ 3ef17b0f-7ba5-4cef-98fe-2abf3711c907
using NonlinearSolve

# ╔═╡ 00996e13-d7a2-4f9e-be84-4eec83428486
using Printf

# ╔═╡ 2aa4e5cf-f6d9-4608-8000-ec078c7ae7b9
using CairoMakie

# ╔═╡ bbaed5c6-445a-43ed-82f2-c7b4b17fbf19
using AlgebraOfGraphics

# ╔═╡ fe20b40a-ce2a-408d-89c4-3bd17255ad62
using LaTeXStrings

# ╔═╡ f0fa9e95-4369-480d-a6d4-d6b545504686
using ProgressMeter

# ╔═╡ 68df8cdc-5652-4fd7-9391-b9de65722f10
using SharedArrays

# ╔═╡ 6e948a61-d8bb-4d93-b6e6-da19256cb24f
using Base.Threads

# ╔═╡ f3d29aa1-ee85-46a9-8258-f6b9f03d87a2
using MLJ

# ╔═╡ 164a20af-4421-4cb6-a4fe-0f8b00ad7059
using SymbolicRegression  # provides SRRegressor via MLJ interface

# ╔═╡ 734d74db-7b19-42fd-89ad-78816a939565
@fromparent begin
	using ..Properties
	using ..Thermal
	using ..Hydraulics
end

# ╔═╡ 4ef8beb3-621c-41f4-97cf-5a634ec87769
# ╠═╡ show_logs = false
using QuasiMonteCarlo  # Sobol sampler

# ╔═╡ 011f212e-63b1-4df3-a646-37bcf23717f8
using Random

# ╔═╡ 12ca0c64-931b-4989-9594-6bd09ed6fea9
# ╠═╡ show_logs = false
using LoopVectorization

# ╔═╡ 0aa504e7-2687-4fae-a9a0-3614455e6db4
using DynamicExpressions


# ╔═╡ 0ef59446-cc78-47f8-8196-32fdfe199980
using DynamicQuantities


# ╔═╡ f6b71da4-db72-4896-9e34-b7afc2e5677c
using Statistics


# ╔═╡ e9a4105f-abec-47ca-9170-ca6aef7a9c49
using MonteCarloMeasurements

# ╔═╡ 4bbe463f-7e79-44f8-8abe-5f104b346e45
using Distributions

# ╔═╡ 5870e251-b3da-42b7-b72f-211e74c669fc
using LinearAlgebra

# ╔═╡ 16436692-04a1-4915-8e89-b8c0805b6bd5
using PlutoUI

# ╔═╡ cedbf94b-e85d-4662-a64e-16fba9ef2ed4
import Base: lock, unlock



# ╔═╡ a2cd3b4a-713b-4d5d-8739-139d0a64214c
println("Threads available: ", nthreads())

# ╔═╡ f3070a9a-d98b-49b8-b300-631a4065211e
set_theme!(
    Theme(
        fontsize = 18,
        Axis = (
            xlabelsize = 22,
            ylabelsize = 22,
            titlesize = 24,
            xgridstyle = :dash,
            ygridstyle = :dash
        ),
        Lines = (linewidth = 3,),
        Scatter = (markersize = 10, strokewidth = 1),
        Legend = (labelsize = 20, framevisible = true)
    )
)


# ╔═╡ 2ab2c592-3e29-4031-8d19-8856f3cdbc6e
# ╠═╡ disabled = true
#=╠═╡
sanity_check_real()
  ╠═╡ =#

# ╔═╡ 004dd778-a212-4147-a9c0-ad60f1eb539b
"""
    cooper_pool_boiling_htc(ρₗ, ρᵥ, h_fg, μₗ, σ, ΔT; g=9.81)

Compute the nucleate-boiling heat-transfer coefficient h_nb [W/m²·K] using the Cooper correlation:

h_nb = 55 * [ g*(ρₗ - ρᵥ)*h_fg*ρₗ^2 / (μₗ*σ) ]^0.33 * ΔT^0.67

# Arguments
- ρₗ    Liquid density [kg/m³]  
- ρᵥ    Vapor density [kg/m³]  
- h_fg  Latent heat of vaporization [J/kg]  
- μₗ    Liquid dynamic viscosity [Pa·s]  
- σ     Surface tension [N/m]  
- ΔT    Wall superheat ΔT_sat = T_wall - T_sat [K]  
- g     Gravity [m/s²] (default 9.81)

# Returns
- h_nb  Nucleate-boiling HTC [W/m²·K]

# References
- Cooper, M. G. (1982). Saturation and Nucleate Pool Boiling – A Simple Correlation.
"""
function cooper_pool_boiling_htc(ρₗ,
                                 ρᵥ,
                                 h_fg,
                                 μₗ,
                                 σ,
                                 ΔT;
                                 g=9.81)
    return 55 * ((g*(ρₗ - ρᵥ)*h_fg*ρₗ^2) / (μₗ*σ))^0.33 * ΔT^0.67
end

# ╔═╡ 63cf7ccc-e17f-48d7-bbf9-2419a6ebe336
# 1) Nucleate‐boiling (Cooper) wrapper for plate channels
function h_nb_plate(ρ_b, ρ_d,    # bubble/liquid and dew/vapor densities
                    h_fg,       # latent heat
                    μ_b,        # bubble‐viscosity (liquid)
                    σ,          # surface tension
                    ΔT;         # wall superheat
                    g=9.81)     # gravity
    return cooper_pool_boiling_htc(
        ρ_b, ρ_d, h_fg, μ_b, σ, ΔT;
        g = g
    )
end

# 2) Convective‐boiling (Amalfi et al.) wrapper for plate channels

# ╔═╡ 6d536fb0-e23f-4d51-8839-8d9dc2422cbd
function h_cb_plate(ṁ,
                    b, np, W, ϕ, β,      # geometry
                    μ_b, k_b, Pr_b,      # liquid properties at bubble
                    h_fg, q_flux,        # latent heat & heat flux
                    ρ_b, ρ_d, x,         # densities + quality
                    λ_pitch)             # corrugation pitch
    D_h  = 2*b
    A_ch = W * 2*b * np * ϕ
    return amalfi_convective_boiling_htc(
        D_h, λ_pitch, β,
        ṁ, A_ch,
        μ_b, k_b, Pr_b,
        h_fg, q_flux,
        ρ_b, ρ_d, x
    )
end

# 3) Smooth‐transition blender between h_nb and h_cb

# ╔═╡ 96465cb7-82db-4287-bdb7-9260d39e12ef
# ╠═╡ disabled = true
#=╠═╡
function h_evap_plate(ṁ,
                      b, np, W, ϕ, β,
                      ρ_b, ρ_d, h_fg, μ_b, σ, ΔT,
                      k_b, Pr_b,
                      q_flux, x, λ_pitch;
                      Re_low=3000, Re_high=8000)

    # hydraulic diameter & channel area
    D_h  = 2*b
    A_ch = W * 2*b * np * ϕ

    # equivalent two‐phase Reynolds
    G_eq  = (ṁ / A_ch)*((1 - x) + x*sqrt(ρ_d/ρ_b))
    Re_eq = G_eq * D_h / μ_b

    # limits
    h_lam  = h_nb_plate(ρ_b, ρ_d, h_fg, μ_b, σ, ΔT)
    h_turb = h_cb_plate(ṁ, b, np, W, ϕ, β,
                        μ_b, k_b, Pr_b,
                        h_fg, q_flux,
                        ρ_b, ρ_d, x,
                        λ_pitch)

    # smooth‐step transition
    ΔRe   = Re_high - Re_low
    center = (Re_low + Re_high)/2
    α     = 2 * erfinv(0.98) / ΔRe
    tf    = 0.5*(1 + erf(α*(Re_eq - center)))

    return (1 - tf)*h_lam + tf*h_turb
end
  ╠═╡ =#

# ╔═╡ 917d72da-5985-42b2-a079-9dd14b97a791

function h_boiling_Huang_Sheer(
   ρₗ, ρᵥ,
   μₗ, λₗ,
   h_vap,
   σ,  Cₚₗ,
   q,  T_sat;
   angle = 35.0
   )

   # Physical constants

   g = 9.81

   # Characteristic bubble diameter

   d_o = 0.0146 * angle * sqrt((g*(ρₗ - ρᵥ)^2) / σ)  # \[m] ([ht.readthedocs.io][1])

   # Liquid thermal diffusivity and Prandtl number

   αₗ  = λₗ / (ρₗ * Cₚₗ)
   Prₗ = μₗ / (ρₗ * αₗ)

   # Break the formula into parts for clarity

   term1 = d_o * λₗ
   term2 = (λₗ * T_sat / (q * d_o))^0.56
   term3 = ((αₗ^2 * h_vap) / (d_o^2))^0.31

   # Final h

   return 1.87e-3 * term1 * term2 * term3 * Prₗ^0.33
end


# ╔═╡ b2fdf520-379d-426a-8271-7170ee7c7ea5
function ExpansionValve(; name,
	fluid= "R134A",
)
	
	ps = @parameters begin
	end
	sts = @variables begin
		logp_i = log(3e5)
		logp_o = log(1e5)
	end
	
	@named i = Properties.RefrigerantProps(fluid = fluid)
	@named o = Properties.RefrigerantProps(fluid = fluid)
	
	eqs = [
		logp_i ~ i.logp,
		logp_o ~ o.logp,
		i.h ~ o.h, # isenthalpic
	    ];
	
	
    # Create the parent system with connected components
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [i, o], name=name)
end

# ╔═╡ 1e443538-f5b0-4de7-8ec1-df8737e5ba70
function test_expansion()
    # Create the borehole model
	@named exp = ExpansionValve()

	eqs=[
		exp.i.h ~ 400e3,
		exp.i.logp ~ log(5e5),
		exp.o.logp ~ log(1e5),
    ]

    
    # Create the parent system with connected components
    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [exp])

    # Simplify and solve
    simpl_sys = structural_simplify(sys)
    problem = NonlinearProblem(simpl_sys, [])
    sol = NonlinearSolve.solve(problem,  NewtonRaphson(autodiff=AutoFiniteDiff()), abstol = 1e-6)
	return SciMLBase.successful_retcode(sol)#.retcode
	
end

# ╔═╡ e4c92fbe-af30-48e9-b501-1ab71bc0d182
@test test_expansion()

# ╔═╡ 39697eef-1b94-46ac-9d75-943997cbb913
function Compressor(; name,
	fluid= "R134A",
	η_is = 0.7,
	V = 0.003/50,
	# pars
	#m = 1.0,
	# vars
	#xy = [0.0, 0.0], u = [0.0, 0.0]
)
	
	ps = @parameters begin
		V = V
		η_is = η_is
	end
	sts = @variables begin
		T_i = 300
		T_o = 300
		T_o_is = 300
		f = 50.0
		W = 0.0
		V̇ = 1.0
		ṁ = 1.0
		logp_i = log(1e5)
		logp_o = log(10e5)
	end
	
	@named i = Properties.RefrigerantProps(fluid = fluid)
	@named o_is = Properties.RefrigerantProps(fluid = fluid)
	@named o = Properties.RefrigerantProps(fluid = fluid)
	
	eqs = [
		# pulling up variables
		T_i ~ i.T,
		T_o ~ o.T,
		T_o_is ~ o_is.T,
		
		# Efficiency definition: ratio of the isentropic enthalpy change to the actual enthalpy change.
		η_is*(o.h - i.h) ~ (o_is.h - i.h),

		logp_i ~ i.logp,
		
		ṁ ~ V̇ * i.ρ,
		
		V̇ ~ f*V,
		
		W ~ ṁ*(o.h-i.h),

		# Enforce the isentropic condition:
		logp_o ~ o.logp,
		
		o_is.logp ~ o.logp,		
		o_is.s ~ i.s,
	    ];
	
	
    # Create the parent system with connected components
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [i, o_is, o], name=name)
end

# ╔═╡ 43489277-7150-4ad1-aabc-b76bd5472b8d
function test_compressor()
    # Create the borehole model
	@named compressor = Compressor()

	eqs=[
		compressor.ṁ ~ 0.03,
		#compressor.f ~ 50,
		compressor.i.h ~ 400e3,
		compressor.logp_i ~ log(2.8e5),
		compressor.logp_o ~ log(10e5),
    ]

    
    # Create the parent system with connected components
    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [compressor])

    # Simplify and solve
    simpl_sys = structural_simplify(sys)
    problem = NonlinearProblem(simpl_sys, [])
    sol = NonlinearSolve.solve(problem,  NewtonRaphson(autodiff=AutoFiniteDiff()), abstol = 1e-6)
	println(sol[compressor.ṁ])
	println(sol[compressor.f])
	println(sol[compressor.W])
	println(sol[compressor.T_i]-273.15)
	println(sol[compressor.T_o]-273.15)
	println(sol[compressor.T_o_is]-273.15)
	println(sol[compressor.o.hd])
	println(sol[compressor.o.hb])
	return SciMLBase.successful_retcode(sol)#.retcode
end

# ╔═╡ 18475dd7-225b-48d5-b88a-81f4223b8dd8
# ╠═╡ disabled = true
#=╠═╡
test_compressor()
  ╠═╡ =#

# ╔═╡ e9829f7f-cfd0-4a4a-ae4c-00ff405cca42
function htc(ρ, ṁ, dh, μ, β, b, np)
	S = (b*np)/2 #assuming half volume/section in each channel
	V̇ = abs(ṁ)/ρ
	
	v = V̇/S                  # speed of the secondary fluid 
	
	Re = ρ * v * dh / μ           # Reynolds
	#very wrong correlation ⚠
	Pr = cp * μ / λ
	Nul = 3.65/(β^0.455 * Re^0.339)
	Nut = 12.6/(β^1.142) * Re^(0.642 + 0.00111*β)
	Nu  = (Nul^3 + Nut^3)^(1/3) * Pr^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
	htc = Nu*λ/dh
end

# ╔═╡ 0c1298e5-440c-467b-b334-a4e0929dd83d
function EvaporatorA(; name,
	refrigerant= "R134A",
	fluid= "MEA",
	χ = 0.28,
	L = 0.2,
	Per = 2*pi*0.05,
)
	
	ps = @parameters begin
		χ = χ
		L = L
		Per = Per
		U_SH = 100000
		U_EV = 100000
		cp_ref_SH = 1000
		cp_sf_SH = 1000
		cp_sf_EV = 1000
		ṁ_ref = 1.0
		ṁ_sf = 1.0
		αL = pi/4
	end
	sts = @variables begin
		T_i_ref = 300
		T_mid_ref = 300
		T_o_ref = 300
		
		h_i_ref = 300
		h_mid_ref = 300
		h_o_ref = 300
		
		T_i_sf = 300
		T_mid_sf = 300
		T_o_sf = 300

		
		p_ref = 1e5
	end

	
	@variables x
	L_EV = sin(αL)^2 * L
	L_SH = cos(αL)^2 * L


	UA_SH = U_SH * L_SH * Per
	P_r = UA_SH / ṁ_ref / (cp_ref_SH)
	α_SH = UA_SH / ṁ_sf / cp_sf_SH
	
	UA_EV = U_EV * L_EV * Per
	α_EV = UA_EV / ṁ_sf / cp_sf_EV
	
	# ---
	
	T_evap = T_mid_ref
	T_b_SH = T_evap + (T_i_sf - T_evap)*exp(-α_SH*(L_SH-x))
	T_b_interface = substitute(T_b_SH, x => 0.0)
	T_b_EV = T_evap + (T_b_interface - T_evap)*exp(-α_EV*(L_EV-x))
	T_r_SH = T_evap + (P_r*(T_i_sf - T_evap)*exp(-α_SH*L_SH))/(P_r + α_SH) * (exp((P_r + α_SH)*x) - 1)/exp(P_r*x)
	T_r_out = substitute(T_r_SH, x => L_SH)
	T_b_out = substitute(T_b_EV, x => 0.0)
	T_r_EV = T_evap
	h_dew = h_mid_ref
	h_r_EV = h_dew - ṁ_sf*cp_sf_EV/ṁ_ref*(T_b_interface-T_b_EV)
	h_r_in = substitute(h_r_EV, x=> 0.0) |> simplify;
	
	@named i_ref = Properties.RefrigerantProps(fluid = refrigerant)
	@named mid_ref = Properties.RefrigerantProps(fluid = refrigerant)
	@named o_ref = Properties.RefrigerantProps(fluid = refrigerant)
	
	@named i_sf = Properties.IncompressibleProps(fluid = fluid)
	@named mid_sf = Properties.IncompressibleProps(fluid = fluid)
	@named o_sf = Properties.IncompressibleProps(fluid = fluid)

	
	eqs = [
		# all pressurees are the same: simplification
		p_ref ~ i_ref.p,
		p_ref ~ o_ref.p,
		p_ref ~ mid_ref.p,

		# Mid point condition on the refrigrant: dew point
		mid_ref.h ~ mid_ref.hd, # this point is fully characterized in h and p
			
		#expressions
		T_b_interface ~ T_mid_sf,
		T_b_out ~ T_o_sf,

		h_r_in ~ h_i_ref,
		h_o_ref ~ 
		
		# pulling up variables for convenience
		T_i_ref ~ i_ref.T,
		T_o_ref ~ o_ref.T,
		T_mid_ref ~ mid_ref.T,
		h_i_ref ~ i_ref.h,
		h_o_ref ~ o_ref.h,
		h_mid_ref ~ mid_ref.h,
		T_i_sf ~ i_sf.T,
		T_o_sf ~ o_sf.T,
		T_mid_sf ~ mid_sf.T,
	]
	


    # Create the parent system with connected components
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [i_ref, mid_ref, o_ref, i_sf, mid_sf, o_sf], name=name)
end

# ╔═╡ 86df6265-c57f-4ef8-bc0c-522957653fba
squash(x) = 1/pi*atan(x) + 1/2

# ╔═╡ 53fb2d9b-13ea-44c5-9cff-03d9562dbfb1
zap(x) = 1/pi*atan(x) + 1/2

# ╔═╡ 2f2fa303-81a3-4bc4-b8a5-caa21851bf23
unzap(x) = tan((x-1/2)*pi)

# ╔═╡ 6a83b5bc-4e25-4781-afa8-8b378d6c838e
"""
SinglePhaseHX: NTU-based model for counter-current single-phase to single-phase heat exchange.

Arguments:
  name: Unique system name
  A: Heat transfer area (default: 1.0)
"""
function SinglePhaseHX(; name)
    # ——— Parameters ———
    ps = @parameters begin
    end

    # ——— State variables ———
    sts = @variables begin
		A = 1.0
		
        # Heat transfer coefficients
        #h_1 = 1000.0      # Heat transfer coefficient side 1 (W/m²·K)
        #h_2 = 1000.0      # Heat transfer coefficient side 2 (W/m²·K)
        U = 500.0         # Overall heat transfer coefficient (W/m²·K)
        
        # Heat capacity rates
        C_1 = 2000.0      # Heat capacity rate side 1 (W/K)
        C_2 = 1500.0      # Heat capacity rate side 2 (W/K)
        C_min = 1500.0    # Minimum heat capacity rate (W/K)
        C_max = 2000.0    # Maximum heat capacity rate (W/K)
        C_r = 0.75        # Heat capacity ratio (C_min/C_max)
        
        # Heat transfer parameters
        NTU = 1.0         # Number of transfer units
        ε = 0.5           # Effectiveness
        
        # Temperatures (K)
        T_in_1 = 320.0    # Inlet temperature side 1
        T_out_1 = 310.0   # Outlet temperature side 1
        T_in_2 = 290.0    # Inlet temperature side 2
        T_out_2 = 300.0   # Outlet temperature side 2
        
        # Heat transfer
        Q = 5000.0        # Heat transfer rate (W)
		
    end
    
    # ——— Governing equations ———
    eqs = [
        # Determine min and max heat capacity rates using ifelse
        C_min ~ ifelse(C_1 <= C_2, C_1, C_2),
        C_max ~ ifelse(C_1 >= C_2, C_1, C_2),
        C_r ~ C_min / C_max,
        
        # NTU definition
        NTU ~ U * A / C_min,
        
        # Effectiveness for counterflow
        ε ~ (1 - exp(-NTU * (1 - C_r))) / (1 - C_r * exp(-NTU * (1 - C_r))),
        
        # Heat transfer calculation using effectiveness
        Q ~ ε * C_min * (T_in_1 - T_in_2),
        
        # Energy balances for both streams
        Q ~ C_1 * (T_in_1 - T_out_1),
        Q ~ C_2 * (T_out_2 - T_in_2)
    ]
    
    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name = name)
    return sys
end

# ╔═╡ 87178d51-5926-4ae4-aea4-cbfc211d7200
#=╠═╡
unknowns(simple)
  ╠═╡ =#

# ╔═╡ 010a0e27-2e68-427b-9fef-569a32590d87
#=╠═╡
full_equations(simple)
  ╠═╡ =#

# ╔═╡ 11df4133-c9f6-436a-9b2b-ada02f872bf0
#=╠═╡
equations(simple)
  ╠═╡ =#

# ╔═╡ 5818a465-6685-40da-b0eb-a7e035c0170c
# ╠═╡ disabled = true
#=╠═╡
test_sys,simple = test_single_phase_hx_temps();
  ╠═╡ =#

# ╔═╡ d70ebc3d-7082-4043-925d-220f0ac753f7
"""
Test function for the SinglePhaseHX model with temperatures
"""
function test_single_phase_hx_temps()
    @named hx = SinglePhaseHX()
    
    # Set basic conditions
    eqs = [
		hx.A ~ 1.0,

		hx.U ~ 800.0,      # Heat transfer coefficient side 2
        
        # Fixed heat capacity rates
        #hx.C_1 ~ 1250,     # Heat capacity rate side 1 (W/K)
        hx.C_2 ~ 1500.0,     # Heat capacity rate side 2 (W/K)
        
        # Inlet temperatures
        hx.T_in_1 ~ 280.0,   # Inlet temperature side 1 (K) (hot side)
        hx.T_out_1 ~ 285.0,   # Inlet temperature side 1 (K) (hot side)
        hx.T_out_2 ~ 290.0    # Inlet temperature side 2 (K) (cold side)
    ]
    
    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [hx])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])
    sol = NonlinearSolve.solve(prob, NewtonRaphson())
    
    println("=== Heat Exchanger Results with Temperatures ===")
    
    println("Heat Capacity Rates:")
    println("  C_1 = ", sol[hx.C_1], " W/K")
    println("  C_2 = ", sol[hx.C_2], " W/K")
    println("  C_min = ", sol[hx.C_min], " W/K")
    println("  C_max = ", sol[hx.C_max], " W/K")
    println("  C_r = ", sol[hx.C_r])
    
    println("Heat Transfer Parameters:")
    println("  A = ", sol.ps[hx.A], " m²")
    println("  NTU = ", sol[hx.NTU])
    println("  ε = ", sol[hx.ε])
    
    println("Temperatures:")
    println("  T_in_1 = ", sol[hx.T_in_1], " K (", sol[hx.T_in_1] - 273.15, " °C)")
    println("  T_out_1 = ", sol[hx.T_out_1], " K (", sol[hx.T_out_1] - 273.15, " °C)")
    println("  ΔT_1 = ", sol[hx.T_in_1] - sol[hx.T_out_1], " K")
    println("  T_in_2 = ", sol[hx.T_in_2], " K (", sol[hx.T_in_2] - 273.15, " °C)")
    println("  T_out_2 = ", sol[hx.T_out_2], " K (", sol[hx.T_out_2] - 273.15, " °C)")
    println("  ΔT_2 = ", sol[hx.T_out_2] - sol[hx.T_in_2], " K")
    
    println("Heat Transfer:")
    println("  Q = ", sol[hx.Q], " W")
    
    # Energy balance verification
    Q1 = sol[hx.C_1] * (sol[hx.T_in_1] - sol[hx.T_out_1])
    Q2 = sol[hx.C_2] * (sol[hx.T_out_2] - sol[hx.T_in_2])
    println("\nEnergy Balance Verification:")
    println("  Q from side 1 = ", Q1, " W")
    println("  Q from side 2 = ", Q2, " W")
    println("  Difference = ", abs(Q1 - Q2), " W")
    
    return sol, simpl
end

# ╔═╡ 65e682c7-af7f-43ec-aa70-e12553026888
md"""
## Same with properties
"""

# ╔═╡ 62cea410-f0af-4219-a7df-92a0e574d74d
# ╠═╡ disabled = true
#=╠═╡
"""
    nusselt_plate(Re, Pr, β)

Compute a combined laminar+turbulent Nusselt number
for a chevron plate channel.
"""
function nusselt_plate(Re::Real, Pr::Real, β::Real)
    # laminar contribution
    Nul = 3.65 / (β^0.455 * Re^0.339)
    # turbulent contribution
    Nut = 12.6 / (β^1.142) * Re^(0.642 + 0.00111*β)
    # blending
    return (Nul^3 + Nut^3)^(1/3)
end
  ╠═╡ =#

# ╔═╡ 2ed52c9e-a7ec-4793-a4c4-4a1a1611636f
# ╠═╡ disabled = true
#=╠═╡
"""
    h_plate(ṁ, ρ, cp, μ, λ, b, np, W, ϕ, β)

Compute heat‐transfer coefficient in a plate channel.
- ṁ: mass‐flow rate (kg/s)
- ρ, cp, μ, λ: fluid properties
- b: half‐channel depth (m)
- np, W, ϕ: plate geometry
- β: chevron angle (rad)
"""
function h_plate(ṁ, ρ, cp, μ, λ, b, np, W, ϕ, β, Re_low=2300, Re_high=4000)
    # 1) Calculate channels for ONE fluid circuit
    n_channels = (np - 1) / 2  # Channels per fluid
    
    # 2) Flow area for ONE fluid (not entire PHE)
    A_ch_fluid = n_channels * W * (2b)  # ϕ NOT included! [m²]
    
    # 3) Velocity and dimensionless numbers
    v = abs(ṁ) / (ρ * A_ch_fluid)
    dh = 4b / ϕ  # Correct hydraulic diameter [m]
    Re = (ρ * v * dh) / μ
    Pr = (cp * μ) / λ
    
    # 4) Calculate blend weights for transition region
    blend = blend_weight(Re, Re_low, Re_high)
    
    # 5) Nusselt numbers for different regimes
    Nu_lam = 3.65 / (β^-0.455 * Re^0.339)  # Laminar
    Nu_turb = 12.6 / (β^-1.142 * Re^(0.642 + 0.00111*β)) # Turbulent
    
    # 6) Blended Nusselt number
    Nu = (1 - blend) * Nu_lam + blend * Nu_turb
    
    return Nu * λ / dh  # [W/m²·K]
end
  ╠═╡ =#

# ╔═╡ 70ca56b4-dcde-44a2-aa6c-cd5b84c28c16
# Define ranges for ṁ and β, fix everything else
ranges = Dict(
  :ṁ => range(0.01, .9, length=20),
  :β  => deg2rad.(range(25, 60, length=5))
)

# ╔═╡ ce7c4e4d-e0ef-4b49-8ae8-0f63acfe4dc4
fixed = Dict(
  :ρ  => 1000.0,
  :cp => 4180.0,
  :μ  => 1e-3,
  :λ  => 0.6,
  :b  => 0.002,
  :np => 12,
  :W  => 0.07,
  :ϕ  => 1.18
)

# ╔═╡ 6be93b5a-be7d-47ad-81ef-72e3d7350fd6
# ╠═╡ disabled = true
#=╠═╡
df = sweep_h_plate(ranges; fixed)
  ╠═╡ =#

# ╔═╡ b41466ba-fcff-4def-957b-2b96032643aa
# ╠═╡ disabled = true
#=╠═╡
simp = test_single_phase_plate_hx()
  ╠═╡ =#

# ╔═╡ 143ccb02-56b0-4b50-a1e5-0d06b53f1b46
# ╠═╡ disabled = true
#=╠═╡
test_PHX_Superheater()
  ╠═╡ =#

# ╔═╡ f7eeabcd-23e6-4733-81aa-a6bbde82837e
md"""
## Base problem for the Phase change NTU
"""

# ╔═╡ 9d502c24-7e38-4c3f-96ed-0e18425d82d9
"""
PhaseChangeHX: NTU-based model for counter-current heat exchange with one stream undergoing phase change (infinite heat capacity).

Arguments:
  name: Unique system name
  A: Heat transfer area (m²)
  U: Overall heat transfer coefficient (W/m²·K)
  C_other: Heat capacity rate of finite side (W/K)
  tol: Tolerance to identify phase-change condition
"""
function PhaseChangeHX(; name, A=1.0, U=500.0, C_other=1500.0, tol=1e-6)
    # ——— Parameters ———
    ps = @parameters begin
    end

    # ——— State variables ———
    sts = @variables begin
        A = 1.0
        U = 1.0
        C_1 = 100.0
		T_in_1 = 300
        T_out_1 = 300
        T_2 = 310
        Q = 10000
        NTU = 1.0
        ε = 0.5
    end

    # ——— Governing equations ———
    eqs = [
        NTU      ~ U * A / C_1,
        ε        ~ 1 - exp(-NTU),
        Q        ~ ε * C_1 * (T_in_1 - T_2),
        Q        ~ C_1 * (T_in_1 - T_out_1),
    ]

    return NonlinearSystem(eqs, sts, ps; name=name)
end

# ╔═╡ 0907ee2b-5ce8-4789-936a-ec14a4c3c059
"""
Test function for the PhaseChangeHX model: checks infinite-capacity behavior.
"""
function test_phase_change_hx()
    @named hx = PhaseChangeHX()

    # Input and parameters
    eqs = [
        hx.A        ~ 1.0,
        hx.U        ~ 800.0,
        hx.C_1  ~ 1500.0,
        hx.T_2  ~ 280.0,
        hx.T_in_1   ~ 300.0
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [];
                                name=:test_pc), [hx])
    simpl = structural_simplify(sys)
    prob  = NonlinearProblem(simpl, [])
    sol   = NonlinearSolve.solve(prob, NewtonRaphson())

    println("=== PhaseChangeHX Test ===")
    println("NTU      = ", sol[hx.NTU])
    println("ε        = ", sol[hx.ε])
    println("Q        = ", sol[hx.Q], " W")
    println("T_sat = ", sol[hx.T_2], " K (== T_in_pc)")
    println("T_in_1  = ", sol[hx.T_in_1], " K")
    println("T_out_1  = ", sol[hx.T_out_1], " K")
    return sol
end

# run test

# ╔═╡ 039f80ce-620b-4863-aa98-1f0f2ed49a7e
# ╠═╡ disabled = true
#=╠═╡
test_phase_change_hx()
  ╠═╡ =#

# ╔═╡ eb386e80-6a11-4b38-b55a-50bbad670474
"""
    blend_weight(Re, Re_low, Re_high)

Smooth weight in [0,1], transitioning over [Re_low, Re_high].
"""
function blend_weight(Re, Re_low, Re_high)
    Re_mid = 0.5*(Re_low + Re_high)
    σ = (Re_high - Re_low) / (2 * erfinv(0.98))
    return 0.5*(1 + erf((Re - Re_mid)/σ))
end

# ╔═╡ e31e4221-54e9-4049-b80e-eca118e9bb00
"""
    h_boiling_Amalfi(m, x, D_h, ρ_l, ρ_g, μ_l, μ_g,
                    κ_l, H_vap, σ, q, A_ch; θ_ch=45.0, θ_max=45.0)

Two-phase boiling heat-transfer coefficient [W/(m²·K)]  
for a plate-and-frame exchanger (Amalfi correlation).

# Arguments
- `m` : mass flow rate [kg/s]
- `x` : vapor quality [–]
- `D_h` : hydraulic diameter [m]
- `ρ_l`, `ρ_g` : liquid and gas densities [kg/m³]
- `μ_l`, `μ_g` : liquid and gas viscosities [Pa·s]
- `κ_l` : liquid thermal conductivity [W/(m·K)]
- `H_vap` : latent heat of vaporization [J/kg]
- `σ` : surface tension [N/m]
- `q` : heat flux [W/m²]
- `S_ch` : channel flow area [m²]
- `θ_ch` : chevron angle [°] (default 45)
- `θ_max` : reference max angle [°] (default 45)

# Returns
- `h` : boiling HT coefficient [W/(m²·K)]
"""
function h_boiling_Amalfi(m, x, D_h,
                          ρ_l, ρ_g,
                          μ_l, μ_g,
                          λ_l, H_vap,
                          σ, q̇, S_ch, β;
						  β_max=pi/4)

    # 1. Mass flux, homogeneous density
    G   = m / S_ch
	ρ_m = 1 / ( x/ρ_g + (1 - x)/ρ_l )
    #ρ_m = x*ρ_g + (1-x)*ρ_g

    # 2. Dimensionless groups
    Bo   = q̇ / (G * H_vap)               # Boiling number
    Bd   = (ρ_l - ρ_g) * 9.81 * D_h^2 / σ # Bond number


    # 3. Two asymptotic Nu correlations (per Amalfi et al.)
    h_tiny = 982 * λ_l / D_h * (β/β_max)^1.101 * (G^2*D_h/ρ_m/σ)^0.315*(ρ_l / ρ_g)^(-0.224)*Bo^0.320
	
    Re_lo    = G * D_h / μ_l
    Re_g     = G * x * D_h / μ_g
    h_macro = 18.495 * λ_l / D_h * (β/β_max)^0.248 * Re_g^0.135 * Re_lo^0.351 * (ρ_g / ρ_l)^(-0.223) * Bd^0.235 * Bo^0.198
    # 4. Smooth blend around Bd ≈ 4 (e.g. from 3.5 to 4.5)
    w  = blend_weight(Bd, 3.9, 4.1)
    h = (1 - w)*h_tiny + w*h_macro

    # 5. Return h = k_l/D_h · Nu
    return h
end

# ╔═╡ eacfb867-b8fe-4166-84df-038848a9d837
function sanity_check_real()
    ### ——— Geometry & Flow from your MTK test ———
    b        = 0.002        # plate spacing [m] (1.5 mm)
    np       = 15            # number of plates
    W        = 0.05          # plate width [m] (10.0 cm)
    L        = 0.2          # plate length [m] (20.0 cm)
    ϕ        = 1.2           # enlargement factor
    β        = π/6           # chevron angle [rad] (30°)
    D_h      = 2*b/ϕ         # hydraulic diameter [m] = 2.5 mm
    CS       = np * b * W    # cross-sectional area [m²] = 3750 mm²
    ṁ        = 0.0174        # refrigerant mass flow [kg/s]
    S_ch     = np*W*ϕ*L      # heat transfer area [m²] = 0.6 m²
    q_flux   = 5810.0        # heat flux [W/m²]
    
    ### ——— R-134a properties @ Tsat = 0°C (from your test) ———
    ρ_l   = 1294.8      # kg/m³, liquid density
    ρ_g   = 14.43       # kg/m³, vapor density
    μ_l   = 281.7e-6    # Pa·s, liquid viscosity (281.7 μPa·s)
    μ_g   = 10.54e-6    # Pa·s, vapor viscosity (10.54 μPa·s)
    κ_l   = 93.4e-3     # W/(m·K), liquid thermal conductivity (93.4 mW/m·K)
    H_vap = 398.6e3     # J/kg, vapor enthalpy (398.6 kJ/kg)
    σ     = 11.428e-3   # N/m, surface tension (11.428 mN/m)

    println("Real‐Case Sanity‐Check: h_boiling_Amalfi")
    println("Geometry:")
    println("  D_h = $(round(D_h*1e3,digits=2)) mm, CS = $(round(CS*1e6,digits=1)) mm²")
    println("  S_ch = $(round(S_ch,digits=2)) m², q″ = $(round(q_flux,digits=0)) W/m²")
    println("  ṁ = $(round(ṁ*1e3,digits=1)) g/s")
    println()
    println("R134A Properties @ 0°C:")
    println("  ρₗ = $(ρ_l) kg/m³, ρᵥ = $(ρ_g) kg/m³")
    println("  μₗ = $(round(μ_l*1e6,digits=1)) μPa·s, μᵥ = $(round(μ_g*1e6,digits=1)) μPa·s")
    println("  κₗ = $(round(κ_l*1e3,digits=1)) mW/m·K, σ = $(round(σ*1e3,digits=3)) mN/m")
    println()

    @printf("%6s | %8s | %12s\n", " x   ", "h [W/m²·K]", "h [kW/m²·K]")
    println("------+----------+--------------")
    for x in (0.2, 0.4, 0.6, 0.8, 1.0)
        h = h_boiling_Amalfi(
            ṁ, x, D_h,
            ρ_l, ρ_g,
            μ_l, μ_g,
            κ_l, H_vap,
            σ, q_flux,
            CS, β
        )
        @printf("  %.1f  |   %6.0f   |     %6.2f\n", x, h, h/1e3)
    end
    
end

# ╔═╡ 07f01adf-8cc7-4a40-aa14-a90cef491b95
"""
    soft_clamp(x, x_low, x_high; tol=0.01)

Smoothly constrain `x` into `[x_low, x_high]` in a C^∞ way with no branching.
- `tol` is the *fraction* of `(x_high - x_low)` over which each edge transition occurs.
  e.g. `tol=0.01` means each “ramp” is 1% of the full width.
"""
function soft_clamp(x::Real, x_low::Real, x_high::Real; tol::Real=0.01)
    # transition half-width
    δ = tol * (x_high - x_low)
    # logistic‐style bumps at the low & high edges:
    s_low  = 1 / (1 + exp(-(x_low + δ - x)  / δ))
    s_high = 1 / (1 + exp(-(x       - (x_high - δ)) / δ))
    # blend:    identity + (x_low - x)·s_low + (x_high - x)·s_high
    return x + (x_low  - x)*s_low  + (x_high - x)*s_high
end

# ╔═╡ a1eec254-6b86-440d-a2a4-2151d588dc23
"""
    phx_friction_factor(Re, β; Re_low = 1800, Re_high = 2200)

Smoothed **Fanning** friction factor for a chevron‑plate channel following
Holger Martin’s VDI‑Heat‑Atlas (2010) correlations.  Martin provides separate
Darcy (
    fᴰₗₐₘ = 64/Re⋅(1 + 0.22/ sinβ)
    fᴰₜᵤʳ = 0.316 Re^–0.25⋅(1 + 2.8/ sinβ)
) with a hard switch at *Re≈2000*.  We blend them over a user‑set window so the
function is smooth.

**Arguments**
* `Re`  – Reynolds number based on hydraulic diameter.
* `β`   – chevron angle [rad].
* `Re_low`, `Re_high` – start/end of the blend window (defaults 1800→2200).

Returns the **Fanning** friction factor 𝑓.
"""
function phx_friction_factor(Re, β, Re_low = 1800, Re_high = 2200)
    β = soft_clamp(β, deg2rad(30), deg2rad(60); tol = 0.02)
    G = sin(β)                      # Martin’s geometry factor
    # Darcy branches (VDI equations N6‑28 & N6‑30)
    fD_lam  = 64 / Re * (1 + 0.22 / G)
    fD_turb = 0.316 * Re^(-0.25) * (1 + 2.8 / G)
    # Smooth blend
    w = blend_weight(Re, Re_low, Re_high)
    fD = (1 - w) * fD_lam + w * fD_turb
    return fD / 4                   # convert to Fanning
end

# ╔═╡ e51c75dc-07f3-4682-9ce2-ed2e72977ea5
"""
    phx_nusselt(Re, Pr, β; Re_low = 1800, Re_high = 2200)

Nusselt number for a chevron plate using Martin’s analogy:

    Nu = 0.205 Pr^{1/3} [ fᴰ Re² sin(2β) ]^{0.374}

The same Darcy friction factor produced by `phx_friction_factor` is used,
ensuring thermo‑hydraulic consistency.
"""
function phx_nusselt(Re, Pr, β, Re_low = 1800, Re_high = 2200)
    fF = phx_friction_factor(Re, β, Re_low, Re_high)
    fD = 4fF
    return 0.205 * Pr^(1/3) * (fD * Re^2 * sin(2β))^0.374
end


# ╔═╡ d31ad2dc-528d-46f6-a19f-90658a1615d0

"""
    h_plate(ṁ, ρ, cp, μ, λ, b, np, W, β; φ = 1.18,
            Re_low = 1800, Re_high = 2200)

Overall single‑phase heat‑transfer coefficient (one side) for a brazed / welded
plate heat exchanger using the Martin correlations.

**Geometry**
* `b`   – half‑gap [m] (one corrugation valley ↔ plate)
* `W`   – plate width [m]
* `np`  – number of plates (=> (np−1)/2 channels per side)
* `φ`   – surface enlargement factor (defaults 1.18)

**Flow / fluid**
* `ṁ`   – mass flow rate [kg s⁻¹]
* `ρ`, `μ`, `cp`, `λ` – liquid properties.
* `β`   – chevron angle [rad]

Returns h [W m⁻² K].
"""
function h_plate(ṁ, ρ, cp, μ, λ, b, np, W, φ, β, Re_low = 1800, Re_high = 2200)
    n_channels = (np - 1) / 2
    A_flow     = n_channels * W * (2b)
    v          = abs(ṁ) / (ρ * A_flow)
    dh         = 2b / φ
    Re         = ρ * v * dh / μ
    Pr         = cp * μ / λ

    Nu = phx_nusselt(Re, Pr, β, Re_low, Re_high)
    return Nu * λ / dh
end


# ╔═╡ 19a9664b-cfa8-4a1a-b2bb-ed3f8c357100
#=╠═╡
"""
    sweep_h_plate(ranges; fixed=Dict(), kwargs...)

Perform a parameter sweep of `h_plate` over the named `ranges`, 
while keeping any other parameters fixed via `fixed`.

# Arguments
- `ranges::Dict{Symbol,AbstractVector}`  
  Each key is a positional argument name of `h_plate`, and its
  value is a vector of values you want to try.
- `fixed::Dict{Symbol,Real}`  (optional)  
  Any other positional arguments you wish to hold constant.
- `kwargs...`  
  Passed through to `h_plate` (e.g. `Re_low`, `Re_high`).

# Returns
A `DataFrame` with one column per input and one column `h` for
the computed heat‑transfer coefficient.
"""
function sweep_h_plate(ranges;
                       fixed=Dict{Symbol,Real}(), kwargs...)
    # positional argument order
    argnames = (:ṁ, :ρ, :cp, :μ, :λ, :b, :np, :W, :ϕ, :β)
    
    # build vectors for each arg in order
    vals = Vector{AbstractVector}(undef, length(argnames))
    for (i,name) in enumerate(argnames)
        if haskey(ranges, name)
            vals[i] = ranges[name]
        elseif haskey(fixed, name)
            vals[i] = [fixed[name]]
        else
            error("You must supply values for $name either in ranges or fixed.")
        end
    end
    
    # prepare output DataFrame
    df = DataFrame()
    for name in argnames
        df[!, Symbol(name)] = Any[]  # placeholder
    end
    df.h = Float64[]
    
    # loop over all combinations
    for combo in Iterators.product(vals...)
        args = Tuple(combo)
        # call h_plate with positional args and any kwargs
        hval = h_plate(args...; kwargs...)
        push!(df, (; 
            ṁ=args[1], ρ=args[2], cp=args[3], μ=args[4], λ=args[5],
            b=args[6], np=args[7], W=args[8], ϕ=args[9], β=args[10],
            h=hval
        ))
    end
    
    return df
end

  ╠═╡ =#

# ╔═╡ e438b966-15eb-42e4-a610-0fc426385bc7

"""
SinglePhasePlateHX_WithProperties: NTU-based counterflow HX with fluid properties
and plate-channel correlations for film coefficients.

Arguments:
  name: system name
  Per: heat-transfer area per length (m²/m)
  b: half-channel depth (m)
  np: number of plates
  W: plate width (m)
  ϕ: surface enlargement factor
  β: chevron angle (rad)
  th: wall thickness (m)
  λ_Cu: wall conductivity (W/m·K)
"""
function SinglePhasePlateHX_WithProperties(; name, fluid1 = "MEA", fluid2 = "MEA")
    # ——— Parameters ———
    pars = @parameters begin
	end

    # ——— State variables ———
    vars = @variables begin
		# parameters but variables in order to solve sizing problems
		b = 2e-3
		np = 30
		W = 10e-2
		ϕ = 1.138
		β = pi/6
		th = 5e-4
		λ_pl = 16

		χ1 = 0.28
		χ2 = 0.28
		
    	L = 1.0                    # plate length (m)
	    γṁ_1 = 0.0
		γṁ_2 = 0.0             # log mass flows (water, MEA)
	    T_in_1 = 300
		T_out_1 = 300        # hot-side temps
	    T_in_2 = 300
		T_out_2 = 300        # cold-side temps
	    Q = 1e4                     # heat transfer rate (W)
	    h_1 = 600
		h_2 = 600
		U =300
		A = 1.0
    end

    # ——— Subsystems ———
    @named base    = SinglePhaseHX()
    @named prop1i  = Properties.IncompressibleProps(fluid=fluid1, χ=χ1)
    @named prop1o  = Properties.IncompressibleProps(fluid=fluid1, χ=χ1)
    @named prop2i  = Properties.IncompressibleProps(fluid=fluid2, χ=χ2)
    @named prop2o  = Properties.IncompressibleProps(fluid=fluid2, χ=χ2)

    # ——— Derived expressions ———
    ṁ_1     = exp(γṁ_1)
    ṁ_2     = exp(γṁ_2)
    dh       = 2*b
	
    # ——— Equations ———
    eqs = [
        # geometry mapping
        A           ~ L*np*W*ϕ,
        A           ~ base.A,
        # property temperatures
        prop1i.T    ~ T_in_1,  prop1i.T    ~ base.T_in_1,
        prop1o.T    ~ T_out_1, prop1o.T    ~ base.T_out_1,
        prop2i.T    ~ T_in_2,  prop2i.T    ~ base.T_in_2,
        prop2o.T    ~ T_out_2, prop2o.T    ~ base.T_out_2,
        # capacity rates
        base.C_1         ~ ṁ_1 * 0.5*(prop1i.cp + prop1o.cp),
        base.C_2         ~ ṁ_2 * 0.5*(prop2i.cp + prop2o.cp),
        # film coefficients
        h_1         ~ h_plate(ṁ_1,
							  0.5*(prop1i.ρ+prop1o.ρ),
							  0.5*(prop1i.cp+prop1o.cp),
							  0.5*(prop1i.μ+prop1o.μ),
							  0.5*(prop1i.λ+prop1o.λ),
							  b, np, W, ϕ, β, 1000, 2000),
        h_2         ~ h_plate(ṁ_2,
							  0.5*(prop2i.ρ+prop2o.ρ),
							  0.5*(prop2i.cp+prop2o.cp),
							  0.5*(prop2i.μ+prop2o.μ),
							  0.5*(prop2i.λ+prop2o.λ),
							  b, np, W, ϕ, β, 1000, 2000),
        # overall coefficient
        U           ~ 1/(1/h_1 + th/λ_pl + 1/h_2),
        U           ~ base.U,
        # heat transfer
        Q           ~ base.Q,
    ]

    sys = NonlinearSystem(eqs, vars, pars, name = name)

    return ModelingToolkit.compose(sys, [base, prop1i, prop1o, prop2i, prop2o], name = name)
end

# ╔═╡ 96966964-d666-4dc2-a554-8c359e3db528

"""
Test function for SinglePhasePlateHX_WithProperties
"""
function test_single_phase_plate_hx()
    @named hx = SinglePhasePlateHX_WithProperties(fluid1 = "MEA", fluid2 = "MEA")

    eqs = [
		hx.b ~ 0.001,
		hx.np ~ 12,
		hx.W ~ 0.05,
		hx.ϕ ~1.0,
		hx.β ~ π/6,
		hx.th ~ 0.0005,
		hx.λ_pl ~ 26,
        # geometry
        hx.L ~ 0.2,
        # mass flows
        hx.γṁ_1 ~ log(0.5),    # 0.1 kg/s water
        hx.γṁ_2 ~ log(0.5),   # 0.15 kg/s MEA
        # inlet temps
        hx.T_in_1 ~ 320.0,
        #hx.T_out_1 ~ 300.0,
        hx.T_in_2 ~ 290.0
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name=:test), [hx])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])
    sol = NonlinearSolve.solve(prob, NewtonRaphson())

    println("=== Plate HX with Properties Test ===")
    println("Q = ", sol[hx.Q], " W")
	
    Q1 = sol[hx.base.C_1] * (sol[hx.base.T_in_1] - sol[hx.base.T_out_1])
    Q2 = sol[hx.base.C_2] * (sol[hx.base.T_out_2] - sol[hx.base.T_in_2])
    println("\nEnergy Balance Verification:")
    println("  Q from side 1 = ", Q1, " W")
    println("  Q from side 2 = ", Q2, " W")
    println("  Difference = ", abs(Q1 - Q2), " W")
	
    println("m1 = ", exp(sol[hx.γṁ_1]), " kg/s")
    println("m2 = ", exp(sol[hx.γṁ_2]), " kg/s")
    println("C1 = ", sol[hx.base.C_1], " ")
    println("C2 = ", sol[hx.base.C_2], " ")
    println("h1 = ", sol[hx.h_1], " W/m2/K")
    println("h2 = ", sol[hx.h_2], " W/m2/K")
    println("U = ", sol[hx.U], " W")
    println("A = ", sol[hx.A], " m2")
    println("eff = ", sol[hx.base.NTU], " %")
    println("eff = ", sol[hx.base.ε], " %")
    println("T_in_1 = ", sol[hx.T_in_1]-273.15, " C")
    println("T_out_1 = ", sol[hx.T_out_1]-273.15, " C")
    println("T_in_2 = ", sol[hx.T_in_2]-273.15, " C")
    println("T_out_2 = ", sol[hx.T_out_2]-273.15, " C")
    return sol
end


# ╔═╡ 9b117e8b-7a38-49ba-9586-f3aacaafa0f3
function PHX_Superheater(; name, sec_fluid = "MEA", refrigerant = "R134A",
                         Re_low=3000,        # blender limits
                         Re_high=8000)
    # ——— Parameters ———
    pars = @parameters begin
        Re_low = Re_low
        Re_high = Re_high
    end

    # ——— State variables ———
    vars = @variables begin
		
		b = 2e-3
		np = 25
		W = 10e-2
		ϕ = 1.138
		β = pi/6
		th = 5e-4
		λ_pl = 16
        λ_pitch = 4e-3

		χ = 0.28
		
        L = 1.0                # plate length (m)
        γṁ_1 = 0.0
        γṁ_2 = 0.0             # log mass flows
        T_in_1 = 300
        T_out_1 = 300          # hot side temps
        T_in_2 = 300
        T_out_2 = 300          # cold/refrig side temps
        Q = 1e4                # heat transfer rate (W)
        h_1 = 600
        h_2 = 600
        U = 300
        A = 1.0
        # refrigerant inlet/outlet enthalpies
        h2i = 300e3
        h2o = 300e3
        logpref = log(1e5)
    end

    # ——— Subsystems ———
    @named base    = SinglePhaseHX()
    @named prop1i  = Properties.IncompressibleProps(fluid=sec_fluid, χ=χ)
    @named prop1o  = Properties.IncompressibleProps(fluid=sec_fluid, χ=χ)
    @named propref = Properties.RefrigerantProps(fluid = refrigerant)

    # ——— Derived expressions ———
    ṁ_1 = exp(γṁ_1)
    ṁ_2 = exp(γṁ_2)
    dh   = 2*b


    # ——— Equationss ———
    eqs = [
        # geometry mapping
        A           ~ L*(np-1)*W*ϕ,
        A           ~ base.A,

		propref.logp ~ logpref,
		propref.h ~ h2i,
		
        # property temperatures mapping
        prop1i.T    ~ T_in_1,
        prop1o.T    ~ T_out_1,
        T_in_2      ~ base.T_in_2,
        T_in_2      ~ propref.T,
        T_out_2     ~ base.T_out_2,
        T_in_1      ~ base.T_in_1,
        T_out_1     ~ base.T_out_1,

        # capacity rates
        base.C_1    ~ ṁ_1 * 0.5*(prop1i.cp + prop1o.cp),
        base.C_2    ~ ṁ_2 * propref.cpd,

        # film coefficients
        h_1 ~ h_plate(
            ṁ_1,
            0.5*(prop1i.ρ+prop1o.ρ),
            0.5*(prop1i.cp+prop1o.cp),
            0.5*(prop1i.μ+prop1o.μ),
            0.5*(prop1i.λ+prop1o.λ),
            b, np, W, ϕ, β, 1000, 2000
        ),

        h_2 ~ h_plate(
            ṁ_2,
            propref.ρ,
            propref.cp,
            propref.η,
            propref.λ,
            b, np, W, ϕ, β, 1000, 2000
        ),


		
        # overall U‐value
        U ~ 1/(1/h_1 + th/λ_pl + 1/h_2),
        U ~ base.U,

        # heat‐transfer balances
        Q       ~ base.Q,
        Q       ~ ṁ_2*(h2o - h2i),
    ]

    sys = NonlinearSystem(eqs, vars, pars, name = name)
    return ModelingToolkit.compose(sys, [base, prop1i, prop1o, propref], name = name)
end


# ╔═╡ 98c20212-3447-481c-b0d4-326a1703997c
"""
Test function for PHX_Superheater
"""
function test_PHX_Superheater()
    # — instantiate with 2 plates, small width for a quick solve
    @named hx = PHX_Superheater(
        Re_low = 3000,
        Re_high= 8000
    )

    # — build a tiny set of boundary conditions to fully specify the superheater
    eqs = [
		
		hx.b ~ 0.001,
		hx.np ~ 25,
		hx.W ~ 0.05,
		hx.ϕ ~1.0,
		hx.β ~ π/6,
		hx.th ~ 0.0005,
		hx.λ_pl ~ 16,
		hx.λ_pitch ~ 26,
		hx.χ ~ 0.28,
        # geometry
        hx.L         ~ 0.2,

        # mass flows (0.5 kg/s each side)
        hx.γṁ_1      ~ log(0.5),
        hx.γṁ_2      ~ log(0.5),

        # hot‐side temps (MEA side)
        hx.T_in_1    ~ 293.15,

        # refrigerant pressure and enthalpies
        hx.logpref      ~ log(2e5),       # 2 bar
        hx.h2i       ~ 400e3,     # inlet enthalpy [J/kg]
        #hx.h2o       ~ 250e3,     # outlet enthalpy [J/kg]
    ]

    # compose and solve
    @named sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name=:test), [hx])
    simpl        = structural_simplify(sys)
    prob         = NonlinearProblem(simpl, [])
    solver       = NewtonRaphson(
                      linesearch = BackTracking(maxstep=20),
                      autodiff    = AutoFiniteDiff()
                   )
    sol = NonlinearSolve.solve(prob, solver; abstol=1e-6)

    # print key results
    println("=== PHX Superheater Test ===")
    println("Q     = ", sol[hx.Q],        " W")
    println("A     = ", sol[hx.A],        " m²")
    println("U     = ", sol[hx.U],        " W/m²·K")
    println("h1     = ", sol[hx.h_1],        " W/m²·K")
    println("h2     = ", sol[hx.h_2],        " W/m²·K")
    println("NTU   = ", sol[hx.base.NTU])
    println("ε     = ", sol[hx.base.ε])
    println("T₁,in = ", sol[hx.T_in_1]-273.15, " °C")
    println("T₁,out= ", sol[hx.T_out_1]-273.15, " °C")
    println("T₂,in = ", sol[hx.T_in_2]-273.15, " °C")
    println("T₂,out= ", sol[hx.T_out_2]-273.15, " °C")
    println("h2i   = ", sol[hx.h2i],       " J/kg")
    println("h2o   = ", sol[hx.h2o],       " J/kg")
    println("p_ref = ", exp(sol[hx.logpref])/1e5,  " bar")
    println("T_sat = ", sol[hx.propref.Tsat]-273.15, " °C")
    println("ρ_b   = ", sol[hx.propref.ρb], " kg/m³")
    println("ρ_d   = ", sol[hx.propref.ρd], " kg/m³")
    println("h_fg  = ", sol[hx.propref.hd - hx.propref.hb], " J/kg")
    return sol
end


# ╔═╡ 11823168-8288-49db-88e8-319df4535273
function PhaseChangePlateHX_WithProperties(; name, sec_fluid = "MEA", refrigerant = "R134A",
                         Re_low=1000,        # blender limits
                         Re_high=2000)


		
    # ——— Parameters ———
    pars = @parameters begin
        Re_low = Re_low
        Re_high = Re_high
    end

    # ——— State variables ———
    vars = @variables begin
		
		b = 2e-3
		np = 30
		W = 7e-2
		ϕ = 1.118
		β = pi/6
		th = 5e-4
		λ_pl = 16
        λ_pitch = 4e-3

		χ = 0.28
		
        L = 1.0
        γṁ_2 = 0.0          # log mass flow for refrigerant
        γṁ_1 = 0.0          # log mass flow for incompressible fluid
        h_2 = 7000.0         # 
        h_2_raw = 7000.0         # 
        h_1 = 2000.0
        T_in_1 = 280.0
        T_out_1 = 280.0
        T_2    = 300.0      # refrigerant bulk T
        logp_2    = log(2.5e5)        # refrigerant pressure
        #X_in    = 0.5
        h2i    = 300e3
        h2o    = 310e3
        Q      = 1000.0
    end

    # ——— Subsystems ———
    @named base   = PhaseChangeHX()
    @named propref= Properties.RefrigerantProps(fluid = refrigerant)
    @named prop1i = Properties.IncompressibleProps(fluid = sec_fluid, χ = χ)
    @named prop1o = Properties.IncompressibleProps(fluid = sec_fluid, χ = χ)

    # ——— Derived expressions ———
    ṁ_1 = exp(γṁ_1)
    ṁ_2 = exp(γṁ_2)
	D_h = 2b/ϕ
	CS = (np-1)/2*b*W

    # ——— Equations ———
    eqs = [
        # geometry
        base.A       ~ L * (np-1) * W * ϕ,

        # incompressible side mapping
        T_in_1       ~ base.T_in_1,
        prop1i.T     ~ T_in_1,
        T_out_1      ~ base.T_out_1,
        prop1o.T     ~ T_out_1,

        # refrigerant side mapping
        base.T_2     ~ propref.Tsat,     # refrigerant bulk = saturation T
        T_2          ~ base.T_2,
        propref.h       ~ h2i,           # set to inlet enthalpy for sat props
        propref.logp       ~ logp_2,           # 

        # capacity rates (incompressible only)
        base.C_1     ~ ṁ_1 * 0.5*(prop1i.cp + prop1o.cp),

		#X_in ~ soft_clamp(propref.X, 0.0, 1.0),#ifelse(propref.X<1.0, ifelse(propref.X>0.0, propref.X, -1.0), 2.0),
		h_2_raw ~ mean(map(x->h_boiling_Amalfi(
           ṁ_2,                                       # mass‐flow [kg/s]
           x,                                       # fixed quality
           D_h,                                     # D_h [m]
           propref.ρb, propref.ρd,                              # ρₗ, ρᵥ [kg/m³]
           propref.ηb, propref.ηd,                         # μₗ, μᵥ [Pa·s]
           propref.λb,                                      # κₗ [W/(m·K)]
           propref.hd,                                  # H_vap [J/kg]
           propref.σ,                                    # σ [N/m]
           abs(Q)/base.A,                             # q″ [W/m²]
           CS,
			β,
		), range(0.2, 1.0, length = 10))),####FIXME
		
		h_2 ~  soft_clamp(h_2_raw, 1000, 20000),
		
        h_1 ~ h_plate(
            ṁ_1,
            0.5*(prop1i.ρ + prop1o.ρ),
            0.5*(prop1i.cp + prop1o.cp),
            0.5*(prop1i.μ + prop1o.μ),
            0.5*(prop1i.λ + prop1o.λ),
            b, np, W, ϕ, β, Re_low, Re_high
        ),

        # overall U linking
        base.U       ~ 1/(1/h_1 + th/λ_pl + 1/h_2),

        # heat-transfer balances
        Q            ~ base.Q,
        Q            ~ ṁ_2*(h2o - h2i),
    ]

    sys = NonlinearSystem(eqs, vars, pars; name = name)
    return ModelingToolkit.compose(sys, [base, propref, prop1i, prop1o]; name = name)
end


# ╔═╡ 409aa2af-c6b7-451c-8557-ced24eccc9c1
"""
Test function for PhaseChangePlateHX_WithProperties:
prints geometry, flows, HTCs –and now all the temperatures in °C – for quick debugging.
"""
function size_pcpwp_simple()
    @named hx = PhaseChangePlateHX_WithProperties(
        Re_low  = 1000.0,
        Re_high = 2000.0
    )

	vars = @variables begin
		
	end
	
    eqs = [
		
		hx.b ~ 0.002,
		hx.np ~ 20,#*factor,
		hx.W ~ 0.07,
		hx.ϕ ~1.18,
		hx.β ~ π/6,
		hx.th ~ 0.5e-3,
		hx.λ_pl ~ 16,
		hx.λ_pitch ~ 4e-3,
		hx.χ ~ 0.28,
        #hx.L          ~ 0.2,              # plate length [m]
		
		hx.γṁ_1       ~ log(0.4),         # brine flow [kg/s]
        hx.γṁ_2       ~ log(0.025),         # refrigerant flow [kg/s]
        hx.T_in_1     ~ 273.15 + 10.0,            # brine inlet [K]
        hx.T_out_1     ~ 273.15 + 7.0,            # brine outlet [K]
        #hx.h2i        ~ 250e3,    #
        hx.h2i        ~ hx.propref.hb,    # 
        hx.h2o        ~  hx.propref.hd,    # start evaporator at bubble‐point enthalpy
    ]

    @named sys = ModelingToolkit.compose(
        NonlinearSystem(eqs, vars, []; name=:test_plate),
        [hx]
    )
    simplified = structural_simplify(sys)

	
    prob = NonlinearProblem(simplified, [])
	
	solver = RobustMultiNewton(autodiff=AutoFiniteDiff());
	
	sol = NonlinearSolve.solve(prob, solver, abstol = 1e-6)
	println(sol[sys.hx.L])
	return sol
end

# ╔═╡ 51de15da-535e-4908-b346-d69d34e2ce01
# ╠═╡ disabled = true
#=╠═╡
sol_size = size_pcpwp_simple()
  ╠═╡ =#

# ╔═╡ 9b9c3b3b-6e9e-49c0-8498-99217be15b2e
#=╠═╡
sol_size.retcode
  ╠═╡ =#

# ╔═╡ fcfeffb2-35aa-4209-a495-ce79949d08e0
#=╠═╡
sol_pcpwp.retcode
  ╠═╡ =#

# ╔═╡ 93114c1c-393c-4569-bcad-936ad8938e44
#=╠═╡
sol_pcpwp.resid
  ╠═╡ =#

# ╔═╡ 86eeb20e-7da2-489d-84af-317f31f63ef4
# ╠═╡ disabled = true
#=╠═╡
sol_pcpwp, simp_pcpwd = test_pcpwp_simple();
  ╠═╡ =#

# ╔═╡ 7d1b4554-571a-4113-afdd-243a24260b7c


# ╔═╡ d01db8eb-17dc-405e-b566-63fa0e86682b
"""
Test function for PhaseChangePlateHX_WithProperties:
prints geometry, flows, HTCs –and now all the temperatures in °C – for quick debugging.
"""
function test_pcpwp_simple()
    @named hx = PhaseChangePlateHX_WithProperties(
        Re_low  = 2000.0,
        Re_high = 2100.0
    )

	vars = @variables begin
		#γfactor = log(10.0)
	end
	#factor = exp(γfactor)
    # drive the model
    eqs = [
		
		hx.b ~ 0.002,
		hx.np ~ 20,#*factor,
		hx.W ~ 0.10,
		hx.ϕ ~1.2,
		hx.β ~ π/6,
		hx.th ~ 0.5e-3,
		hx.λ_pl ~ 16,
		hx.λ_pitch ~ 4e-3,
		hx.χ ~ 0.28,
        hx.L          ~ 0.2,              # plate length [m]
        
		
		#hx.γṁ_1       ~ log(0.9),         # brine flow [kg/s]
        hx.γṁ_2       ~ log(0.025),         # refrigerant flow [kg/s]
        hx.T_in_1     ~ 273.15 + 10.0,            # brine inlet [K]
        hx.T_out_1     ~ 273.15 + 7.0,            # brine outlet [K]
        #hx.h2i        ~ 250e3,    #
        hx.h2i        ~ hx.propref.hb,    # 
        hx.h2o        ~  hx.propref.hd,    # start evaporator at bubble‐point enthalpy
        #hx.p_2  ~ 2.8e5               # refrigerant pressure [Pa]
        #hx.propref.Tsat  ~ hx.T_in_1 - 10,               # refrigerant evaporation temperature
    ]

    @named sys = ModelingToolkit.compose(
        NonlinearSystem(eqs, vars, []; name=:test_plate),
        [hx]
    )
    simplified = structural_simplify(sys)

	
    prob = NonlinearProblem(simplified, [])
	
	solver = RobustMultiNewton(autodiff=AutoFiniteDiff());
	sol = NonlinearSolve.solve(prob, solver, abstol = 1e-6)
	#sol = NonlinearSolve.solve(prob, abstol = 1e-6)
	
    # Extract geometry
    b = sol[hx.b]
    np = sol[hx.np]
    W = sol[hx.W]
    ϕ = sol[hx.ϕ]
    β_deg = sol[hx.β] * 180/π
    th = sol[hx.th]
    
    # Extract flows and heat transfer
    ṁ1 = exp(sol[hx.γṁ_1])
    ṁ2 = exp(sol[hx.γṁ_2])
    A = sol[hx.base.A]
    Q = sol[hx.base.Q]
    ho = sol[hx.h2o] 
    hi = sol[hx.h2i] 
    Δh = ho - hi 
    
    # Extract temperatures
    T_in_1 = sol[hx.T_in_1] - 273.15      # brine inlet °C
    T_out_1 = sol[hx.T_out_1] - 273.15    # brine outlet °C
    T_sat = sol[hx.propref.Tsat] - 273.15 # evaporation temperature °C
    
    # Extract pressure
    pev = sol[hx.propref.p]/1e5
    
    # Extract HTCs
    h1 = sol[hx.h_1]
    h2 = sol[hx.h_2]
    h2_raw = sol[hx.h_2_raw]

    # ===== EXTRACT ALL REFRIGERANT PROPERTIES USED IN h2_raw =====
    ρb = sol[hx.propref.ρb]        # bubble point density [kg/m³]
    ρd = sol[hx.propref.ρd]        # dew point density [kg/m³] 
    ηb = sol[hx.propref.ηb]        # bubble point viscosity [Pa·s]
    ηd = sol[hx.propref.ηd]        # dew point viscosity [Pa·s]
    λb = sol[hx.propref.λb]        # bubble point thermal conductivity [W/(m·K)]
    hd = sol[hx.propref.hd]        # dew point enthalpy [J/kg]
    σ = sol[hx.propref.σ]          # surface tension [N/m]
    
    # Additional useful properties
    Tsat_K = sol[hx.propref.Tsat]  # saturation temperature [K]
    p_ref = sol[hx.propref.p]      # refrigerant pressure [Pa]
    
    # Derived quantities used in h2_raw calculation
    D_h = 2 * b / ϕ               # hydraulic diameter [m]
    CS = np * b * W               # cross-sectional area [m²]
    q_flux = abs(Q) / A           # heat flux [W/m²]

    # ===== COMPREHENSIVE OUTPUT =====
    println("=== Enhanced PhaseChangePlateHX Test ===")
    println()
    println("GEOMETRY:")
    println("  b       = ", round(b*1000, digits=2), " mm (plate spacing)")
    println("  np      = ", round(Int, np), " (number of plates)")
    println("  W       = ", round(W*100, digits=1), " cm (plate width)")
    println("  L       = ", round(sol[hx.L]*100, digits=1), " cm (plate length)")
    println("  ϕ       = ", round(ϕ, digits=3), " (enlargement factor)")
    println("  β       = ", round(β_deg, digits=1), "° (chevron angle)")
    println("  th      = ", round(th*1000, digits=2), " mm (plate thickness)")
    println("  D_h     = ", round(D_h*1000, digits=2), " mm (hydraulic diameter)")
    println("  CS      = ", round(CS*1e6, digits=2), " mm² (cross-sectional area)")
    println("  A       = ", round(A, digits=3), " m² (heat transfer area)")
    println()
    
    println("FLOWS & OPERATING CONDITIONS:")
    println("  ṁ₁ (brine)     = ", round(ṁ1, digits=3), " kg/s")
    println("  ṁ₂ (refrig)    = ", round(ṁ2, digits=4), " kg/s")
    println("  Q              = ", round(Q, digits=0), " W")
    println("  q″             = ", round(q_flux, digits=0), " W/m² (heat flux)")
    println("  Δh (refrig)    = ", round(Δh/1000, digits=1), " kJ/kg")
    println("  p_evap         = ", round(pev, digits=4), " bar")
    println()
    
    println("TEMPERATURES:")
    println("  T_brine_in     = ", round(T_in_1, digits=2), " °C")
    println("  T_brine_out    = ", round(T_out_1, digits=2), " °C")
    println("  T_evaporation  = ", round(T_sat, digits=4), " °C")
    println("  p_evaporation  = ", round(p_ref/1e5, digits=4), " bar")
    println("  X  = ", round(sol[hx.propref.X]/1e2, digits=4) , " %")
	
    println("  ΔT_approach    = ", round(T_in_1 - T_sat, digits=2), " K")
    println()
    
    println("=== REFRIGERANT PROPERTIES (inputs to h2_raw calculation) ===")
    println()
    println("BUBBLE POINT (liquid) PROPERTIES:")
    println("  ρₗ = ", round(ρb, digits=1), " kg/m³ (liquid density)")
    println("  ηₗ = ", round(ηb*1e6, digits=1), " μPa·s (liquid viscosity)")
    println("  λₗ = ", round(λb*1000, digits=1), " mW/(m·K) (liquid thermal conductivity)")
    println()
    
    println("DEW POINT (vapor) PROPERTIES:")
    println("  ρᵥ = ", round(ρd, digits=2), " kg/m³ (vapor density)")
    println("  ηᵥ = ", round(ηd*1e6, digits=2), " μPa·s (vapor viscosity)")
    println("  hᵥ = ", round(hd/1000, digits=1), " kJ/kg (vapor enthalpy)")
    println()
    
    println("INTERFACE PROPERTIES:")
    println("  σ = ", round(σ*1000, digits=3), " mN/m (surface tension)")
    println()
    
    println("HEAT TRANSFER COEFFICIENTS:")
    println("  h₁ (brine)      = ", round(h1, digits=0), " W/(m²·K)")
    println("  h₂_raw (calc)   = ", round(h2_raw, digits=0), " W/(m²·K) [unclamped]")
    println("  h₂ (final)      = ", round(h2, digits=0), " W/(m²·K) [clamped]")
    println("  Clamp effect    = ", round((h2_raw-h2)/h2_raw*100, digits=1), "% reduction")
    println()
    
    println("QUALITY RANGE FOR h2_raw (0.2 to 1.0):")
    for x in [0.2, 0.4, 0.6, 0.8, 1.0]
        h_x = h_boiling_Amalfi(ṁ2, x, D_h, ρb, ρd, ηb, ηd, λb, hd, σ, q_flux, CS, sol[hx.β])
        println("  x = ", x, " → h = ", round(h_x, digits=0), " W/(m²·K)")
    end

    return sol, simplified
end

# ╔═╡ 9418798c-4603-498a-ab7a-a45d509aa67f
#=╠═╡
sol_pcpwp.retcode
  ╠═╡ =#

# ╔═╡ 42bda7b4-ea1c-4717-a514-a6d97aa68f93
#=╠═╡
sol_pcpwp.resid
  ╠═╡ =#

# ╔═╡ 1bd62dd3-e84e-4e0b-8731-890cc33286b8
function Evaporator(; name,
        sec_fluid="MEA", χ=0.28, refrigerant = "R134A")
    # ——— Parameters ———
    pars = @parameters begin
    end

    # ——— State variables ———
    vars = @variables begin
		
		b = 2e-3
		np = 30
		W = 10e-2
		ϕ = 1.138
		β = pi/6
		th = 5e-4
		λ_pl = 16
        λ_pitch = 4e-3

		χ = 0.28
		
		L = 0.2
		
        γṁ_2 = 0.0          # log mass flow for ref
        γṁ_1 = 0.0          # log mass flow for incompressible fluid

		uz_fr_sh = 0.0 # fraction of superheater, unbounded variable R -> 0 1 interval, 0.0 -> 0.5

        #Q = 1e4
        h2i = 600.0
        h2mid = 600.0
        h2o = 600.0
        # PhaseChangeHX variables
        T_in_1 = 280.0
        T_mid_1 = 280.0
        T_out_1 = 280.0
        T_out_2 = 280.0
        T_evap = 280.0
        
		logp_2 = log(1e5)

		Q = 3000
    end

    # ——— Derived expressions ———
	fr_sh = zap(uz_fr_sh)
	fr_ev = 1-zap(uz_fr_sh)
    ṁ_1 = exp(γṁ_1)
    ṁ_2 = exp(γṁ_2)
	#p_2 = exp(logp_2)
    dh     = 2*b

    @named ev = PhaseChangePlateHX_WithProperties(sec_fluid = sec_fluid, refrigerant = refrigerant)

    @named sh = PHX_Superheater(sec_fluid = sec_fluid, refrigerant = refrigerant)

    eqs = [
		b ~ ev.b,
		b ~ sh.b,
		
		np ~ ev.np,
		np ~ sh.np,
		
		W ~ ev.W,
		W ~ sh.W,
		
		ϕ ~ ev.ϕ,
		ϕ ~ sh.ϕ,
		
		β ~ ev.β,
		β ~ sh.β,
		
		th ~ ev.th,
		th ~ sh.th,
		
		λ_pl ~ ev.λ_pl,
		λ_pl ~ sh.λ_pl,
        
		λ_pitch ~ ev.λ_pitch,
        λ_pitch ~ sh.λ_pitch, 

		χ ~ ev.χ,
		χ ~ sh.χ,
    
		# geometry
        sh.L ~ L * fr_sh,
        ev.L ~ L * fr_ev,
		
        # mass flows
        sh.γṁ_1 ~ γṁ_1,
        sh.γṁ_2 ~ γṁ_2,
    	ev.γṁ_1 ~ γṁ_1,
        ev.γṁ_2 ~ γṁ_2,
		
        # inlet temps
        sh.T_in_1 ~ T_in_1,
        sh.T_out_1 ~ T_mid_1,
		sh.T_out_2 ~ T_out_2,
        ev.T_in_1 ~ T_mid_1,
        ev.T_out_1 ~ T_out_1,

		
		sh.logpref ~ logp_2,
		ev.logp_2 ~ logp_2,
		
		sh.h2o ~ h2o,
		sh.h2i ~ h2mid,
		ev.h2o ~ h2mid,
		ev.h2i ~ h2i,
		T_evap ~ ev.propref.Tsat,
		
		h2mid ~ sh.propref.hd, # dew point
		
		Q ~ ev.Q + sh.Q,
    ]
	

    sys = NonlinearSystem(eqs, vars, pars; name=name)
    return ModelingToolkit.compose(sys, [ev, sh]; name=name)
end


# ╔═╡ 2d05d7cf-fa80-4552-a2ca-a5ceac3939f1
"""
Test function for SinglePhasePlateHX_WithProperties
"""
function test_ev()
    @named hx = Evaporator()

    eqs = [
		
		hx.b ~ 0.002,
		hx.np ~ 25,#*factor,
		hx.W ~ 0.10,
		hx.ϕ ~1.2,
		hx.β ~ π/6,
		hx.th ~ 0.5e-3,
		hx.λ_pl ~ 16,
		hx.λ_pitch ~ 4e-3,
		hx.χ ~ 0.28,
		hx.L ~ 0.2,
		
        #hx.γṁ_1   ~ log(.365),
        hx.γṁ_2   ~ log(0.025),
        hx.T_in_1 ~ 273.15 + 10,
        hx.T_out_1 ~ hx.T_in_1 - 3,
        hx.h2i ~ hx.ev.propref.hb,
        #hx.h2i ~ 200.0e3,
        hx.T_out_2 ~ hx.T_evap + 7,
		#hx.logp_2 ~ log(2.8e5),
    ]

    @named sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], []; name=:test_plate), [hx])
    simplified = structural_simplify(sys)
	#return simplified
    prob = NonlinearProblem(simplified, [])

	solver = NewtonRaphson(; 
		linesearch = BackTracking(; maxstep=20), 
	#	autodiff = AutoFiniteDiff()
	)

	sol = NonlinearSolve.solve(prob, solver, abstol = 1e-5)
	#sol = NonlinearSolve.solve(prob, abstol = 1e-6)

    println("=== Plate HX with Properties Test ===")
	
    println("Return code: ", sol.retcode)
    println("Q = ", sol[hx.Q], " W")
    println("Q  ev = ", sol[hx.ev.Q], " W")
    println("Q sh = ", sol[hx.sh.Q], " W")
    println("Tin brine= ", sol[hx.T_in_1]-273.15, " C")
    println("T mid brine= ", sol[hx.T_mid_1]-273.15, " C")
    println("T out brine= ", sol[hx.T_out_1]-273.15, " C")
    println("T evaporation = ", sol[hx.sh.propref.Tsat]-273.15, " C")
    println("T out ev = ", sol[hx.sh.T_out_2]-273.15, " C")
    println("h in = ", sol[hx.h2i], " J/kg")
    println("hmid = ", sol[hx.h2mid], " J/kg")
    println("h out = ", sol[hx.h2o], " J/kg")
    println("=== Length fraction ===")
    println("L = ", sol[hx.L], " m")
    println("EV = ", (1-zap(sol[hx.uz_fr_sh]))*100, " %")
    println("SH = ", zap(sol[hx.uz_fr_sh])*100, " %")
    println("=== HTC ===")
    println("EV h1 = ", sol[hx.ev.h_1], " W/m2/K")
    println("EV h2 = ", sol[hx.ev.h_2], " W/m2/K")
    println("SH h1 = ", sol[hx.sh.h_1], " W/m2/K")
    println("SH h2 = ", sol[hx.sh.h_2], " W/m2/K")
    println("=== Flow ===")
    println("p2 = ", exp(sol[hx.logp_2]), " Pa")
    println("p2 = ", sol[hx.sh.propref.p], " Pa")
	
    println("ṁ 1 = ", exp(sol[hx.γṁ_1]), " kg/s")
    println("ṁ 2 = ", exp(sol[hx.γṁ_2]), " kg/s")
    return sol
end

# ╔═╡ 6e831997-cf76-4ea4-815f-8f5ffd48d3e1
# ╠═╡ disabled = true
#=╠═╡
sol_evap = test_ev()
  ╠═╡ =#

# ╔═╡ e8abe6ea-9eda-4f83-abdc-da2bd91dd26a
#=╠═╡
sol_evap.retcode
  ╠═╡ =#

# ╔═╡ 6e3b84d1-afe1-4c56-98ba-50b2f13cfc46
#=╠═╡
sol_evap.resid
  ╠═╡ =#

# ╔═╡ c81e49af-a294-45d2-b2fd-c605dad70ef3
#=╠═╡
sol_evap.stats
  ╠═╡ =#

# ╔═╡ 1adba409-b8cb-430d-9ab6-c29aad55035f

# 1) Factor out your problem‐construction into its own function:
function make_evap_comp_prob()
    @named hx = Evaporator()
    @named compressor = Compressor()
    
    eqs = [
        # ===== EVAPORATOR SETUP =====
        hx.b ~ 0.002,
        hx.np ~ 25,
        hx.W ~ 0.10,
        hx.ϕ ~ 1.2,
        hx.β ~ π/6,
        hx.th ~ 0.5e-3,
        hx.λ_pl ~ 16,
        hx.λ_pitch ~ 4e-3,
        hx.χ ~ 0.28,
        hx.L ~ 0.2,
        
        # Fix evaporator mass flow (like the working test_ev)
        #hx.γṁ_2 ~ log(0.025),
        hx.T_in_1 ~ 273.15 + 10,
        hx.T_out_1 ~ hx.T_in_1 - 3,
        hx.h2i ~ hx.ev.propref.hb,
        hx.T_out_2 ~ hx.T_evap + 7,
        
        # ===== COUPLING =====
        exp(hx.γṁ_2) ~ compressor.ṁ,           # Mass flow coupling
        compressor.logp_i ~ hx.logp_2,         # Pressure connection
        compressor.i.h ~ hx.h2o,               # Enthalpy connection
        compressor.logp_o ~ log(10e5),         # Outlet pressure
        compressor.f ~ 28.81,                   # REMOVED - let frequency be solved
    ]
	
    sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name=:evap_comp), [hx, compressor])
    simp  = structural_simplify(sys)
    NonlinearProblem(simp, [])            # <-- no u0/p here yet
end

# ╔═╡ fddf235e-b540-4bc9-99d6-b75a7c49bdf6
"""
    solve_with_cache(prob; solver=TrustRegion(autodiff=AutoFiniteDiff()),
                          sol=nothing,
                          kwargs...)

Solve `prob` with `solver`, passing through any `kwargs` to `solve`.  
If you pass a previous `sol` (a `NonlinearSolution`), its `.u` is used as the new `u0` for a warm start.
"""
function solve_with_cache(prob;
                          solver = TrustRegion(autodiff = AutoFiniteDiff()),
                          sol    = nothing,
                          kwargs...)
    # pick up last solution vector if given
    u0 = isnothing(sol) ? nothing : sol.u

    # delegate to NonlinearSolve.solve
    return solve(prob, solver; u0 = u0, kwargs...)
end


# ╔═╡ 6466aac2-8c05-4bea-a0bc-d2c4723529e0
prob_1 = make_evap_comp_prob();

# ╔═╡ 192d94f7-1dba-4ccf-8a14-7f1e090e84f5
solver_1 = TrustRegion(autodiff = AutoFiniteDiff());


# ╔═╡ e4f83718-32dd-4170-bfae-36c5edd28997
# ╠═╡ disabled = true
#=╠═╡
sol_5 = NonlinearSolve.solve(prob_1, solver_1)
  ╠═╡ =#

# ╔═╡ e5614c61-6260-414f-b100-3a3683de0524
#=╠═╡
sol_5.u
  ╠═╡ =#

# ╔═╡ 832139e7-83e3-496e-8e84-3499c4b9c04b
function ThermodynamicCycle(; name, refrigerant="R134A", fluid="MEA")
    # —— Variables with initial guesses ——
    vars = @variables begin
        # Geometry variables
        b = 0.002
        np = 25
        W = 0.10
        ϕ = 1.2
        β = π/6
        th = 0.5e-3
        λ_pl = 16
        λ_pitch = 4e-3
        L = 0.2
		χ = 0.28
        
        # Boundary/state variables
        T_i_b = 283.15   # 10°C
		T_o_b = T_i_b - 3
        logp_l = log(2e5)
        logp_h = log(10e5)
		γṁ_sf = log(0.4)
		γṁ_ref = (0.04)
        f = 28.81         # Compressor frequency
        h_i_r = 200e3       # Refrigerant inlet enthalpy
        h_o_r = 300e3       # Refrigerant inlet enthalpy

		
        #T_i_b ~ cycle.T_i_b,
        #T_m_b ~ cycle.T_o_b,
        #T_o_b ~ cycle.T_o_b,
        #h_i_r ~ cycle.h_i_r,
        #h_o_r ~ cycle.h_o_r,
			
		log_ΔTSH = log(7)
		log_ΔTSC = log(3)
		log_ΔTSF = log(3)
		
		Q_evap = 4000
		W_comp = 1000
		Q_cond = -5000

		COP1 = 5
    end

    # —— Subcomponents ——
    @named evaporator = Evaporator(sec_fluid=fluid, refrigerant=refrigerant)
    @named compressor = Compressor(fluid=refrigerant)
	@named inlet_evap = Properties.RefrigerantProps(fluid = refrigerant)
	@named inlet_comp = Properties.RefrigerantProps(fluid = refrigerant)
	@named outlet_cond = Properties.RefrigerantProps(fluid = refrigerant)
	@named outlet_comp = Properties.RefrigerantProps(fluid = refrigerant)

    # —— System equations (remove redundant χ) ——
    eqs = [
        # Geometry setup
        evaporator.b ~ b,
        evaporator.np ~ np,
        evaporator.W ~ W,
        evaporator.ϕ ~ ϕ,
        evaporator.β ~ β,
        evaporator.th ~ th,
        evaporator.λ_pl ~ λ_pl,
        evaporator.λ_pitch ~ λ_pitch,
        evaporator.L ~ L,
        evaporator.χ ~ χ,
        
        # Evaporator boundary conditions
        evaporator.T_in_1 ~ T_i_b,
        evaporator.T_out_1 ~ T_o_b,
		
        evaporator.T_out_1 ~ T_i_b - exp(log_ΔTSF),   # Fixed temperature drop
        evaporator.h2i ~ h_i_r,
        evaporator.h2o ~ h_o_r,
        evaporator.T_out_2 ~ evaporator.T_evap + exp(log_ΔTSH),  # Superheat
        
        # Component coupling
        exp(γṁ_ref) ~ compressor.ṁ,           # Mass flow connection
		logp_l ~ compressor.logp_i,
        compressor.logp_i ~ evaporator.logp_2, # Pressure connection
        compressor.i.h ~ evaporator.h2o,       # Enthalpy connection
        compressor.logp_o ~ logp_h,           # Outlet pressure
        
        h_i_r ~ outlet_cond.h,
		
        # Compressor control
        compressor.f ~ f,
        
        # Connect to evaporator's mass flow
        evaporator.γṁ_2 ~ γṁ_ref,
        evaporator.γṁ_1 ~ γṁ_sf,

		inlet_evap.h ~ h_i_r,
		inlet_evap.logp ~ logp_h,

		inlet_comp.h ~ compressor.i.h,
		inlet_comp.logp ~ logp_h,

		outlet_cond.T ~ outlet_cond.Tsat - exp(log_ΔTSC),
		outlet_cond.logp ~ logp_h,

		outlet_comp.h ~ compressor.o.h,
		outlet_comp.logp ~ logp_h,

		
		
		Q_evap ~ exp(γṁ_ref) * (inlet_comp.h - inlet_evap.h),
		W_comp ~ exp(γṁ_ref) * (outlet_comp.h - inlet_comp.h),
		Q_cond ~ exp(γṁ_ref) * (outlet_cond.h - outlet_comp.h),
			
		COP1 ~ -Q_cond/W_comp,


    ]

    # —— Assemble system ——
    sys = NonlinearSystem(eqs, vars, [], systems=[evaporator, compressor, inlet_evap, outlet_cond, outlet_comp, inlet_comp], name=name)
    return sys
end

# ╔═╡ e09caba7-c987-4c7b-af32-beec786c2504
#=╠═╡
SciMLBase.successful_retcode(sol_cycle.retcode)
  ╠═╡ =#

# ╔═╡ fbfcb718-1d8d-462f-8c9d-c7e340152d50
#=╠═╡
sol_cycle.retcode
  ╠═╡ =#

# ╔═╡ 3d841b73-85ed-40dd-b78a-4c85dc471164
#=╠═╡
println.(sol_cycle.resid)
  ╠═╡ =#

# ╔═╡ 93a62e8e-200c-401f-bcbf-31d4fa7c9275
#=╠═╡
sol_cycle.stats
  ╠═╡ =#

# ╔═╡ 7ac5bbd6-d4a9-4b17-bd0e-b473d97c8a57
# ╠═╡ disabled = true
#=╠═╡
sol_cycle = test_thermodynamic_cycle()
  ╠═╡ =#

# ╔═╡ b12af1c4-53f5-49ed-94f9-dd9f8434bcf0
#=╠═╡
sol_cycle.retcode
  ╠═╡ =#

# ╔═╡ 5eae6b7f-ee19-4f7e-8943-e8828ef027fb
function test_thermodynamic_cycle()
    refrigerant = "R134A"
    fluid = "MEA"

    @named cycle = ThermodynamicCycle(refrigerant=refrigerant, fluid=fluid)

    # Second solve: Fixed frequency
    eqs2 = [
        cycle.b ~ 0.002,
        cycle.np ~ 25,
        cycle.W ~ 0.10,
        cycle.ϕ ~ 1.2,
        cycle.β ~ π/6,
        cycle.th ~ 0.5e-3,
        cycle.λ_pl ~ 16,
        cycle.λ_pitch ~ 4e-3,
        cycle.χ ~ 0.28,
        cycle.L ~ 0.2,
        cycle.T_i_b ~ 283.15,
        cycle.logp_h ~ log(20e5),
		
		cycle.log_ΔTSH ~ log(7),
		cycle.log_ΔTSC ~ log(3),
		cycle.log_ΔTSF ~ log(3),
        cycle.f ~ 50,  # Fix frequency to found value
    ]

    sys2 = ModelingToolkit.compose(NonlinearSystem(eqs2, [], [], name=:sys2), [cycle])
    @time simplified2 = structural_simplify(sys2)
    
    # Create proper initial guess for second system
    u0_second = [
        cycle.γṁ_ref => log(0.043)  # Previous mass flow as initial guess
        # Other variables will use default initial conditions
    ]

	
    solver       = NewtonRaphson(
                      linesearch = BackTracking(maxstep=20),
                      autodiff    = AutoFiniteDiff()
                   )
    prob2 = NonlinearProblem(simplified2, u0_second)
    @time sol = NonlinearSolve.solve(prob2, solver, abstol=1e-5)
    println("KPI")
    println(sol[cycle.Q_evap])
    println(sol[cycle.Q_cond])
    println(sol[cycle.W_comp])
    println(sol[cycle.COP1])
    println("Refrigerant")
	
    println(exp(sol[cycle.γṁ_ref]))
    println(sol[cycle.inlet_comp.h])
    println(sol[cycle.inlet_evap.h])
    println(sol[cycle.outlet_comp.h])
    println(sol[cycle.outlet_cond.h])
    println("brine")
    println(exp(sol[cycle.γṁ_sf]), " kg/s")
    println(sol[cycle.evaporator.T_in_1]-273.15, " C")
    println(sol[cycle.evaporator.T_mid_1]-273.15, " C")
    println(sol[cycle.evaporator.T_out_1]-273.15, " C")
    return sol
end

# ╔═╡ f66a54da-5889-4dfb-a093-5a8eb8964b42
# 1) Factor out your problem‐construction into its own function:
function make_thermodynamic_cycle_prob()
    refrigerant = "R134A"
    fluid = "MEA"

    @named cycle = ThermodynamicCycle(refrigerant=refrigerant, fluid=fluid)

    # Second solve: Fixed frequency
    eqs2 = [
        cycle.b ~ 0.002,
        cycle.np ~ 25,
        cycle.W ~ 0.10,
        cycle.ϕ ~ 1.2,
        cycle.β ~ π/6,
        cycle.th ~ 0.5e-3,
        cycle.λ_pl ~ 16,
        cycle.λ_pitch ~ 4e-3,
        cycle.χ ~ 0.28,
        cycle.L ~ 0.2,
        cycle.T_i_b ~ 283.15,
        cycle.logp_h ~ log(20e5),
		
		cycle.log_ΔTSH ~ log(7),
		cycle.log_ΔTSC ~ log(3),
		cycle.log_ΔTSF ~ log(3),
        cycle.f ~ 50,  # Fix frequency to found value
    ]

    sys2 = ModelingToolkit.compose(NonlinearSystem(eqs2, [], [], name=:sys2), [cycle])
    @time simplified2 = structural_simplify(sys2)
    
    # Create proper initial guess for second system
    u0_second = [
        cycle.γṁ_ref => log(0.043)  # Previous mass flow as initial guess
        # Other variables will use default initial conditions
    ]

	
    solver       = NewtonRaphson(
                      linesearch = BackTracking(maxstep=20),
                      autodiff    = AutoFiniteDiff()
                   )
    prob2 = NonlinearProblem(simplified2, u0_second)
end

# ╔═╡ d2b14ef0-9b53-4846-a169-27de59cb9f0f
prob_cycle = make_thermodynamic_cycle_prob()

# ╔═╡ d3cc6b23-4a6d-451c-8a6e-81bf1b71abfb
# ╠═╡ disabled = true
#=╠═╡
solver3       = NewtonRaphson(
                      linesearch = BackTracking(maxstep=20),
                      autodiff    = AutoFiniteDiff()
                   )
  ╠═╡ =#

# ╔═╡ d930fddf-92b8-468c-a556-f0126a93eb43
#=╠═╡
sol_4 = NonlinearSolve.solve(prob_cycle, solver3)
  ╠═╡ =#

# ╔═╡ b5d46480-3583-4eff-aeec-861e2b1bd1cb
#=╠═╡
sol_4.retcode
  ╠═╡ =#

# ╔═╡ 757853b5-e844-432a-9538-b7afab10bb72
#=╠═╡
sol_5 = solve_with_cache(prob_cycle, solver = solver3, sol=sol_cycle)
  ╠═╡ =#

# ╔═╡ b9a9a9a0-dca2-45ee-8250-24bb680d839e
# ╠═╡ disabled = true
#=╠═╡
"""
    solve_with_cache(prob; solver=TrustRegion(autodiff=AutoFiniteDiff()),
                          sol=nothing,
                          kwargs...)

Solve `prob` with `solver`, passing through any `kwargs` to `solve`.  
If you pass a previous `sol` (a `NonlinearSolution`), its `.u` is used as the new `u0` for a warm start.
"""
function solve_with_cache(prob;
                          solver = TrustRegion(autodiff = AutoFiniteDiff()),
                          sol    = nothing,
                          kwargs...)
    # pick up last solution vector if given
    u0 = isnothing(sol) ? nothing : sol.u

    # delegate to NonlinearSolve.solve
    return solve(prob, solver; u0 = u0, kwargs...)
end
  ╠═╡ =#

# ╔═╡ e1c53c82-3ede-40a7-ae3f-e848fc763372


# ╔═╡ f60e1452-e13c-472d-85bb-fbb64c3e1311


# ╔═╡ 364eaab7-4812-4467-8319-777f0b12d1c7


# ╔═╡ 819264a6-e4b7-406c-9d97-8ac0f36c3ea8
"""

"""
function GSHP(; name,
                     refrigerant = "R134A",
                     fluid       = "MEA",
                     χ           = 0.28,)
    # ——— Parameters ———
    ps = @parameters begin
    end

    # ——— State variables ———
    sts = @variables begin
		T_i_b = 9+273.15
		T_o_b = 9+273.15
		h_i_r = 200e3
		h_o_r = 200e3
		log_ΔTSF = 3
		log_ΔTSH = 7
		log_ΔTSC = 5
		γṁ_sf = 0.1
		γṁ_ref = 0.1
		#T_cond = 45
		#p_cond = 2e5
		logp_l = 1e5
		logp_h = 10e5

		W_comp = 1000
		P_elec_pump = 1000
		Q_evap = 2000
		Q_cond = -3000
		COP1 = 5
		COP2 = 5
    end

    # ——— Property components ———
    @named cycle  = ThermodynamicCycle(fluid = fluid, refrigerant = refrigerant)
    @named bhe = Thermal.UTubeGrouted(depth=250.0, D_p1_o=32e-3, D_p1_i=28.24e-3, 
                                  D_p2_o=32e-3, D_p2_i=28.24e-3,
                                  D_b=0.115, λ_pipe_1=0.4, λ_pipe_2=0.4, 
                                  λ_grout=1.5, λ_soil=2.0,
                                  T_wall=273.15 + 10, s=0.075, χ=0.23)

    # Create the base model with a standard pump
    @named pump = Hydraulics.DynamicPropBasicPump(
        H_shutoff = 12,
        Q_max = 13/3600,
        fluid = fluid
    )



    @named phx = Hydraulics.DynamicPropPlateHeatExchanger(
        fluid = "MEA"
    )
    
    # ——— Governing equations ———
	
    eqs = [

        T_i_b ~ cycle.T_i_b,
        T_o_b ~ cycle.T_o_b,
        h_i_r ~ cycle.h_i_r,
        h_o_r ~ cycle.h_o_r,
        
		log_ΔTSF ~ cycle.log_ΔTSF,
        log_ΔTSH ~ cycle.log_ΔTSH,
        log_ΔTSC ~ cycle.log_ΔTSC,
		
        γṁ_sf ~ cycle.γṁ_sf,
        γṁ_ref ~ cycle.γṁ_ref,
        logp_l ~ cycle.logp_l,
        logp_h ~ cycle.logp_h,
        W_comp ~ cycle.W_comp,
        P_elec_pump ~ pump.base.P_elec,
        Q_evap ~ cycle.Q_evap,
        Q_cond ~ cycle.Q_cond,
        COP1 ~ cycle.COP1,
        COP2 ~ -Q_cond/(W_comp+P_elec_pump),
		
		bhe.T_2 ~ cycle.T_i_b,
		cycle.T_o_b ~ bhe.T_1,
		exp(cycle.γṁ_sf) ~ bhe.ṁ,
		pump.props.T ~ bhe.T_1,
		pump.base.ṁ ~ bhe.ṁ,
		phx.base.ṁ ~ bhe.ṁ,
		pump.base.Δp ~ 1e5,#phx.phx.Δp,
		phx.props.T ~ bhe.T_1,
    ]

    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name = name)
    return ModelingToolkit.compose(sys,
        [cycle, bhe, pump, phx], name = name)
end

# ╔═╡ 80885b7b-9ab2-4010-a954-ed4790c277f0
"""
Test routine for the GSHP model
"""
function test_GSHP()
    
    @named hp = GSHP(χ=0.280);
    eqs = [
        
        hp.cycle.evaporator.b ~ 1.5e-3,
        hp.cycle.evaporator.np ~ 20,
        hp.cycle.evaporator.W ~ 0.10,
        hp.cycle.evaporator.ϕ ~ 1.18,
        hp.cycle.evaporator.β ~ π/6,
        hp.cycle.evaporator.th ~ 0.5e-3,
        hp.cycle.evaporator.λ_pl ~ 16,
        hp.cycle.evaporator.λ_pitch ~ 4e-3,
        hp.cycle.evaporator.χ ~ 0.28,
        hp.cycle.evaporator.L ~ 0.32,

		
		hp.phx.L ~ 0.3,
		hp.phx.W ~ 0.10,
		hp.phx.b ~ 1.5e-3,
		hp.phx.np ~ 20,
		hp.phx.ϕ ~ 1.18,
		hp.phx.β ~ π/6,
		
        hp.cycle.f ~ 50,
        hp.cycle.log_ΔTSF ~ log(3.5),
        hp.cycle.log_ΔTSH ~ log(7),
        hp.cycle.log_ΔTSC ~ log(5),
        hp.cycle.logp_h ~ log(14e5),
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [hp])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])
    solver = NewtonRaphson(; 
        linesearch = BackTracking(; maxstep=20), 
        autodiff = AutoFiniteDiff()
    )
	
	solver = RobustMultiNewton(
        autodiff = AutoFiniteDiff())
	FastShortcutNLLSPolyalg
    sol = NonlinearSolve.solve(prob, solver; abstol = 1e-6)

    # Pick out the symbols you care about:
    vars = (
      # refrigerant side enthalpies
      h_in  = hp.cycle.inlet_evap.h,
      h_out = hp.cycle.outlet_comp.h,
      # borehole fluid temperatures
      T_in  = hp.cycle.T_i_b,
      T_out = hp.cycle.T_o_b,
      # mass‐flows
      m_sf  = hp.cycle.γṁ_sf,
      m_ref = hp.cycle.γṁ_ref,
      # pressures
      p_low = hp.cycle.logp_l,
      p_high= hp.cycle.logp_h,
      # loads & COP
      Q_ev  = hp.cycle.Q_evap,
      Q_co  = hp.cycle.Q_cond,
      COP1  = hp.cycle.COP1,
      COP2  = hp.COP2,
    )

    println("=== GSHP Steady‑State Results (retcode=$(sol.retcode)) ===\n")

    for (name, sym) in pairs(vars)
        val = sol[sym]                    # raw value
        # some sensible unit conversions:
        display = begin
          if name in (:h_in, :h_out)
            round(val/1e3, digits=2) , "kJ/kg"
          elseif name in (:T_in, :T_out)
            round(val - 273.15, digits=2) , "°C"
          elseif name in (:p_low, :p_high)
            round(exp(sol[sym]) / 1e5, digits=3) , "bar"
          elseif name in (:m_sf, :m_ref)
            round(exp(val), digits=4) , "kg/s"
          elseif name in (:Q_ev, :Q_co)
            round(val/1e3, digits=2) , "kW"
          else
            round(val, digits=3)
          end
        end
        println(rpad(string(name),12), ": ", display...)
    end
	println(hp.phx.Δp)
    return sol.retcode
end

# ╔═╡ 3ad72686-9a35-4486-a71c-180cba7d46a2
# ╠═╡ disabled = true
#=╠═╡
test_GSHP()
  ╠═╡ =#

# ╔═╡ 241265f4-f66f-4158-855f-89f002f23d4e
"""

"""
function GSHP_stripped(; name,
                     refrigerant = "R134A",
                     fluid       = "MEA",
                     χ           = 0.28,)
    # ——— Parameters ———
    ps = @parameters begin
    end

    # ——— State variables ———
    sts = @variables begin
		T_i_b = 9+273.15
		T_o_b = 9+273.15
		h_i_r = 200e3
		h_o_r = 200e3
		log_ΔTSF = 3
		log_ΔTSH = 7
		log_ΔTSC = 5
		γṁ_sf = 0.1
		γṁ_ref = 0.1
		#T_cond = 45
		#p_cond = 2e5
		logp_l = 1e5
		logp_h = 10e5

		W_comp = 1000
		Q_evap = 2000
		Q_cond = -3000
		COP1 = 5
    end

    # ——— Property components ———
    @named cycle  = ThermodynamicCycle(fluid = fluid, refrigerant = refrigerant)
    @named bhe = Thermal.UTubeGrouted(depth=250.0, D_p1_o=32e-3, D_p1_i=28.24e-3, 
                                  D_p2_o=32e-3, D_p2_i=28.24e-3,
                                  D_b=0.115, λ_pipe_1=0.4, λ_pipe_2=0.4, 
                                  λ_grout=1.5, λ_soil=2.0,
                                  T_wall=273.15 + 10, s=0.075, χ=0.23)

    
    # ——— Governing equations ———
	
    eqs = [

        T_i_b ~ cycle.T_i_b,
        T_o_b ~ cycle.T_o_b,
        h_i_r ~ cycle.h_i_r,
        h_o_r ~ cycle.h_o_r,
        
		log_ΔTSF ~ cycle.log_ΔTSF,
        log_ΔTSH ~ cycle.log_ΔTSH,
        log_ΔTSC ~ cycle.log_ΔTSC,
		
        γṁ_sf ~ cycle.γṁ_sf,
        γṁ_ref ~ cycle.γṁ_ref,
        logp_l ~ cycle.logp_l,
        logp_h ~ cycle.logp_h,
        W_comp ~ cycle.W_comp,
        Q_evap ~ cycle.Q_evap,
        Q_cond ~ cycle.Q_cond,
        COP1 ~ cycle.COP1,
		
		bhe.T_2 ~ cycle.T_i_b,
		cycle.T_o_b ~ bhe.T_1,
		exp(cycle.γṁ_sf) ~ bhe.ṁ,
    ]

    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name = name)
    return ModelingToolkit.compose(sys,
        [cycle, bhe], name = name)
end

# ╔═╡ 9aa6b70a-9464-4f48-a340-61030ca9310b
# ╠═╡ disabled = true
#=╠═╡
# 1. Get stripped solution
sol_test_GSHP = test_GSHP_stripped()

# 2. Build thermal system
  ╠═╡ =#

# ╔═╡ 5f7f03fb-20e2-4bd2-8bb7-f09a8c2d94e6
#=╠═╡
sol_test_GSHP.u
  ╠═╡ =#

# ╔═╡ 2a4c04e6-875e-453a-9be5-6194a5530c06
function GSHP_thermal(; name,
                     refrigerant = "R134A",
                     fluid       = "MEA",)
    # ——— Parameters ———


	
    ps = @parameters begin
        # Heat exchanger parameters
        b = 1.5e-3
        np = 20
        W = 0.10
        ϕ = 1.18
        β = π/6
        th = 0.5e-3
        λ_pl = 16
        λ_pitch = 4e-3
        χ = 0.28
        L = 0.32
        
        # BHE parameters
        L_BHE = 250
        T_wall = 273.15 + 10
        D_p1_o = 32e-3
        D_p1_i = 28.24e-3
        D_p2_o = 32e-3
        D_p2_i = 28.24e-3
        D_b = 0.115
        λ_pipe_1 = 0.4
        λ_pipe_2 = 0.4
        λ_grout = 1.5
        λ_soil = 2.0
        α = 0.5
        ε = 1.5e-6
        U_radius = 0.064
        Re_low = 2300
        Re_high = 4000
        g = 9.81

		
        # Cycle parameters
        f_ = 50
        log_ΔTSF_ = log(3.5)
        log_ΔTSH_ = log(7)
        log_ΔTSC_ = log(5)
        logp_h_ = log(14e5)
    end

    # ——— State variables ———
    sts = @variables begin
        # Original cycle variables
        T_i_b = 9 + 273.15
        T_o_b = 9 + 273.15
        h_i_r = 200e3
        h_o_r = 200e3
        γṁ_sf = 0.1
        γṁ_ref = 0.1
        logp_l = 1e5
        
        # Performance metrics
        W_comp = 1000
        Q_evap = 2000
        Q_cond = -3000
        COP1 = 5
        # Cycle parameters
        f = 50
        log_ΔTSF = log(3.5)
        log_ΔTSH = log(7)
        log_ΔTSC = log(5)
        logp_h = log(14e5)
        
    end

    # ——— Property components ———
    @named cycle = ThermodynamicCycle(fluid=fluid, refrigerant=refrigerant)
    @named bhe_thermal = Thermal.UTubeGrouted(
        depth = L_BHE,
        D_p1_o = D_p1_o,
        D_p1_i = D_p1_i,
        D_p2_o = D_p2_o,
        D_p2_i = D_p2_i,
        D_b = D_b,
        λ_pipe_1 = λ_pipe_1,
        λ_pipe_2 = λ_pipe_2,
        λ_grout = λ_grout,
        λ_soil = λ_soil,
        T_wall = T_wall,
        α = α,
        χ = χ
    )

    
    # ——— Governing equations ———
    eqs = [
        # ===== ORIGINAL THERMODYNAMIC CYCLE =====
        # Heat exchanger geometry
        cycle.evaporator.b ~ b,
        cycle.evaporator.np ~ np,
        cycle.evaporator.W ~ W,
        cycle.evaporator.ϕ ~ ϕ,
        cycle.evaporator.β ~ β,
        cycle.evaporator.th ~ th,
        cycle.evaporator.λ_pl ~ λ_pl,
        cycle.evaporator.λ_pitch ~ λ_pitch,
        cycle.evaporator.χ ~ χ,
        cycle.evaporator.L ~ L,
        
        
        # Cycle connections
        cycle.f ~ f,
        T_i_b ~ cycle.T_i_b,
        T_o_b ~ cycle.T_o_b,
        h_i_r ~ cycle.h_i_r,
        h_o_r ~ cycle.h_o_r,
        log_ΔTSF ~ cycle.log_ΔTSF,
        log_ΔTSH ~ cycle.log_ΔTSH,
        log_ΔTSC ~ cycle.log_ΔTSC,
        γṁ_sf ~ cycle.γṁ_sf,
        γṁ_ref ~ cycle.γṁ_ref,
        logp_l ~ cycle.logp_l,
        logp_h ~ cycle.logp_h,
        W_comp ~ cycle.W_comp,
        
        # Performance metrics
        Q_evap ~ cycle.Q_evap,
        Q_cond ~ cycle.Q_cond,
        COP1 ~ cycle.COP1,
        
        # Thermal BHE connections (original)
        bhe_thermal.T_2 ~ cycle.T_i_b,
        cycle.T_o_b ~ bhe_thermal.T_1,
        bhe_thermal.ṁ ~ exp(γṁ_sf),

		
        f ~ f_,
        log_ΔTSF ~ log_ΔTSF_,
        log_ΔTSH ~ log_ΔTSH_,
        log_ΔTSC ~ log_ΔTSC_,
        logp_h ~ logp_h_,
        
    ]

    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [cycle, bhe_thermal])
end

# ╔═╡ 22e21ab4-f5a9-4ccd-8fdf-280279e12096
function GSHP_thermal_msf(; name,
                          refrigerant = "R134A",
                          fluid       = "MEA")

    # ——— Parameters ———
    ps = @parameters begin
        # Heat exchanger parameters
        b = 1.5e-3
        np = 20
        W = 0.10
        ϕ = 1.18
        β = π/6
        th = 0.5e-3
        λ_pl = 16
        λ_pitch = 4e-3
        χ = 0.28
        L = 0.32

        # BHE parameters
        L_BHE = 250
        T_wall = 273.15 + 10
        D_p1_o = 32e-3
        D_p1_i = 28.24e-3
        D_p2_o = 32e-3
        D_p2_i = 28.24e-3
        D_b = 0.115
        λ_pipe_1 = 0.4
        λ_pipe_2 = 0.4
        λ_grout = 1.5
        λ_soil = 2.0
        α = 0.5
        ε = 1.5e-6
        U_radius = 0.064
        Re_low = 2300
        Re_high = 4000
        g = 9.81

        # Cycle parameters
        f_ = 50
        log_ΔTSF_ = log(3.5)
        log_ΔTSH_ = log(7)
        log_ΔTSC_ = log(5)
        logp_h_ = log(14e5)

        # NEW: make source-fluid mass flow a parameter (kg/s)
        ṁ_sf_ = 0.5
    end

    # ——— State variables ———
    sts = @variables begin
        # Original cycle variables
        T_i_b = 9 + 273.15
        T_o_b = 9 + 273.15
        h_i_r = 200e3
        h_o_r = 200e3
        γṁ_sf = 0.1          # log(ṁ_sf); kept as a state to connect with the cycle
        γṁ_ref = 0.1
        logp_l = 1e5

        # Performance metrics
        W_comp = 1000
        Q_evap = 2000
        Q_cond = -3000
        COP1 = 5

        # Cycle parameters
        f = 50
        log_ΔTSF = log(3.5)
        log_ΔTSH = log(7)
        log_ΔTSC = log(5)
        logp_h = log(14e5)
    end

    # ——— Property components ———
    @named cycle = ThermodynamicCycle(fluid=fluid, refrigerant=refrigerant)
    @named bhe_thermal = Thermal.UTubeGrouted(
        depth      = L_BHE,
        D_p1_o     = D_p1_o,
        D_p1_i     = D_p1_i,
        D_p2_o     = D_p2_o,
        D_p2_i     = D_p2_i,
        D_b        = D_b,
        λ_pipe_1   = λ_pipe_1,
        λ_pipe_2   = λ_pipe_2,
        λ_grout    = λ_grout,
        λ_soil     = λ_soil,
        T_wall     = T_wall,
        α          = α,
        χ          = χ
    )

    # ——— Governing equations ———
    eqs = [
        # ===== THERMODYNAMIC CYCLE =====
        # Heat exchanger geometry
        cycle.evaporator.b       ~ b,
        cycle.evaporator.np      ~ np,
        cycle.evaporator.W       ~ W,
        cycle.evaporator.ϕ       ~ ϕ,
        cycle.evaporator.β       ~ β,
        cycle.evaporator.th      ~ th,
        cycle.evaporator.λ_pl    ~ λ_pl,
        cycle.evaporator.λ_pitch ~ λ_pitch,
        cycle.evaporator.χ       ~ χ,
        cycle.evaporator.L       ~ L,

        # Connections
        cycle.f      ~ f,
        T_i_b        ~ cycle.T_i_b,
        T_o_b        ~ cycle.T_o_b,
        h_i_r        ~ cycle.h_i_r,
        h_o_r        ~ cycle.h_o_r,
        log_ΔTSF     ~ cycle.log_ΔTSF,
        log_ΔTSH     ~ cycle.log_ΔTSH,
        log_ΔTSC     ~ cycle.log_ΔTSC,
        γṁ_sf       ~ cycle.γṁ_sf,
        γṁ_ref      ~ cycle.γṁ_ref,
        logp_l       ~ cycle.logp_l,
        logp_h       ~ cycle.logp_h,
        W_comp       ~ cycle.W_comp,

        # Performance metrics
        Q_evap ~ cycle.Q_evap,
        Q_cond ~ cycle.Q_cond,
        COP1   ~ cycle.COP1,

        # ===== BHE THERMAL COUPLING =====
        bhe_thermal.T_2 ~ cycle.T_i_b,
        cycle.T_o_b     ~ bhe_thermal.T_1,

        # CHANGE: drive BHE with parameterized mass flow directly
        bhe_thermal.ṁ  ~ ṁ_sf_,

        # And bind the cycle’s log-massflow to this parameter
        γṁ_sf          ~ log(ṁ_sf_),

        # Fixed cycle parameters
        f          ~ f_,
        log_ΔTSF   ~ log_ΔTSF_,
        log_ΔTSH   ~ log_ΔTSH_,
        log_ΔTSC   ~ log_ΔTSC_,
        logp_h     ~ logp_h_
    ]

    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [cycle, bhe_thermal])
end


# ╔═╡ 8d3a070d-682b-47f4-bfd6-9feb53ce9ccc
# ╠═╡ disabled = true
#=╠═╡
begin    
	
	solver1 = RobustMultiNewton(
        autodiff = AutoFiniteDiff())
	
	# 1a) Build the system
	@named sys = GSHP_thermal()               # uses your defaults dict internally
	
	# 1b) Simplify/compile
	simp1 = structural_simplify(sys)
	
	# 1d) Create the NonlinearProblem, passing p0 as the 3rd argument
	prob1 = NonlinearProblem(simp1, nothing, [])
	
	# 1e) Solve it
	solth = solve(prob1, solver_robust, abstol = 1e-6)
end
  ╠═╡ =#

# ╔═╡ 1ac95f91-934a-4dc7-97a8-f6cc13ecb864
#=╠═╡
solth.retcode
  ╠═╡ =#

# ╔═╡ df4eb258-955d-4f9e-bfdd-6d202c5b8d90
# ╠═╡ disabled = true
#=╠═╡
p01 = [
	  #:b       => 1.5e-3,
	  #:np      => 20.0,
	  #:W       => 0.10,
	  #:ϕ       => 1.18,
	  #:β       => π/6,
	  #:th      => 0.5e-3,
	  #:λ_pl    => 16.0,
	  #:λ_pitch => 4e-3,
	  #:χ       => 0.28,
	  :log_ΔTSF_       => log(5.55),
	  #:f       => 45.0
	]
  ╠═╡ =#

# ╔═╡ 79b5e151-be4a-4d38-ae46-c975b4bf557e
#=╠═╡
# 1d) Create the NonlinearProblem, passing p0 as the 3rd argument
	prob01 = remake(prob1, u0 = solth.u, p=p01)
  ╠═╡ =#

# ╔═╡ 04ec6d0d-d812-464e-b244-8458207283da
#=╠═╡
solth01 = solve(prob01, solver_robust, abstol = 1e-6)
  ╠═╡ =#

# ╔═╡ a0e018cb-a75a-4844-8716-98492b29b6d9
#=╠═╡
solth01.retcode
  ╠═╡ =#

# ╔═╡ b9a1f8e0-7ca4-4435-b2fa-741b0a25da82
"""
Create a complete geothermal loop with a pump and BHE.
Implements parameter inheritance for easy remaking.
"""
function GSHP_hydraulics(; name::Symbol, fluid::String="MEA")
    
    # Define the loop-level parameters
    pars = @parameters begin

        # Heat exchanger parameters
        b = 1.5e-3
        np = 20
        W = 0.10
        ϕ = 1.18
        β = π/6
        χ = 0.28
        L = 0.32
        
        # BHE parameters
        L_BHE = 250
        T_wall = 273.15 + 10
        D_p1_o = 32e-3
        D_p1_i = 28.24e-3
        D_p2_o = 32e-3
        D_p2_i = 28.24e-3
        ε = 1.5e-6
        U_radius = 0.064
        Re_low = 2300
        Re_high = 4000
        g = 9.81

		T_leg1 = 300 
		T_leg2 = 300
		T_bend = 300 

		ṁ = 0.2

	Q_max = 1/1000
	H_shutoff = 12
    end
    
    # Create the pump
    @named pump = Hydraulics.DynamicPropBasicPump(
        H_shutoff = H_shutoff,
        Q_max = Q_max,
        fluid = fluid,
		N_max = 4980/60,
        χ = χ
    )
    
	@named bhe = Hydraulics.BHE_Upipe_Hydraulic_Complete(
		depth = L_BHE,
		Di = D_p1_i,
		fluid = fluid,
		χ = χ
	)
    
	@named phx = Hydraulics.DynamicPropPlateHeatExchanger(
		fluid = fluid,
		χ = χ,
		g = g,
		Re_low = Re_low,
		Re_high = Re_high,
	)


    # Create system-level variables
    vars = @variables begin
		ṁ
    end
    
	eqs = [
		# Connect temperatures
		bhe.T_leg1 ~ T_leg1,
		bhe.T_bend ~ T_bend,
		bhe.T_leg2 ~ T_leg2,

		bhe.ṁ ~ ṁ,
		
		phx.L ~ L,
        phx.W ~ W,
        phx.b ~ b,
        phx.np ~ np,
        phx.ϕ ~ ϕ,
        phx.β ~ β,
		phx.props.T ~ T_leg2,  # 20°C
        
		phx.base.ṁ ~ ṁ,

		pump.props.T ~ T_leg2,
        
		pump.base.ṁ ~ ṁ,
        pump.base.Δp ~ bhe.Δp + phx.base.Δp,
    ]

    
    parent_system = NonlinearSystem(eqs, vars, pars, name=name)
    
    # Compose the complete system
    return compose(parent_system, [pump, bhe, phx])
end

# ╔═╡ fbc3e326-5684-4029-ba49-6c0f065d151b
function solve_with_fallback(prob, solvers; abstol=1e-6, max_retries=1, u0=nothing)
    # Use provided initial guess or problem default
    current_prob = isnothing(u0) ? prob : remake(prob; u0=u0)
    
    last_error = nothing
    last_solver = nothing
    
    for solver in solvers
        try
            sol = solve(current_prob, solver; abstol)
            if SciMLBase.successful_retcode(sol)
                return sol
            else
                @warn "Solver $(typeof(solver)) returned non-success: $(sol.retcode)"
            end
        catch e
            @warn "Solver $(typeof(solver)) failed: $(sprint(showerror, e))"
            last_error = e
            last_solver = typeof(solver)
        end
    end
    
    # Final attempt with increased iterations
    @warn "Trying final fallback with increased iterations"
    try
        return solve(current_prob, last(solvers); abstol=abstol, maxiters=200)
    catch
        error("Complete failure after all attempts")
    end
end

# 2. Parameter update with solution propagation

# ╔═╡ 7ef1715e-2fd3-4deb-853d-233496735474
function solve_with_new_param(prob_base, sol_prev, pars_new, solvers)
	if isnothing(sol_prev)
    	prob_new = remake(prob_base; p=pars_new)
	else
		prob_new = remake(prob_base; u0=sol_prev.u, p=pars_new)
	end
    return solve_with_fallback(prob_new, solvers)
end

# 3. Thermal problem solver

# ╔═╡ 88624925-6092-44c8-9539-28527aa43a7c
function solve_thermal_problem(prob_thermal, previous_sol, pars, solvers)
    return solve_with_new_param(prob_thermal, previous_sol, pars, solvers)
end

# 4. Hydraulic problem solver

# ╔═╡ 526d0dad-939f-48af-a09f-c1cae30d6770
function solve_hydraulic_problem(prob_hydraulic, previous_sol, params, thermal_results, sys_thermal, solvers)
    # Extract needed values from thermal solution
    ṁ = exp(thermal_results[sys_thermal.γṁ_sf])
    T_in = thermal_results[sys_thermal.cycle.T_i_b]
    T_out = thermal_results[sys_thermal.cycle.T_o_b]
    
    combined_params = [
        params...,
        :ṁ => ṁ,
        :T_leg1 => T_out,  # Outlet from heat pump to BHE
        :T_leg2 => T_in,   # Inlet to heat pump from BHE
        :T_bend => (T_in + T_out)/2
    ]
    
    return solve_with_new_param(prob_hydraulic, previous_sol, combined_params, solvers)
end

# 5. Integrated system simulation

# ╔═╡ f6aed83e-02a5-4d62-874e-55f333266b42
function simulate_gshp_system(thermal_prob, hydraulic_prob, sys_thermal, 
                              pars, previous_sol_thermal, previous_sol_hydraulic, 
                              solvers)
    # Solve thermal system
    sol_thermal = solve_thermal_problem(thermal_prob, previous_sol_thermal, pars, solvers)
    
    # Solve hydraulic system
    sol_hydraulic = solve_hydraulic_problem(
        hydraulic_prob, 
        previous_sol_hydraulic, 
        pars, 
        sol_thermal, 
        sys_thermal, 
        solvers
    )
    
    return (sol_thermal, sol_hydraulic)
end

# 6. Parameter sweep driver

# ╔═╡ 34d568ce-25d5-4a0b-b1db-0ca3b8e8607a
function parameter_sweep(thermal_prob, hydraulic_prob, sys_thermal, 
                         base_params, param_sym, values, solvers; previous_sol_thermal=nothing, previous_sol_hydraulic=nothing)
    results_thermal = []
    results_hydraulic = []
    
    # Solve first point
    sol_thermal_prev = solve_thermal_problem(
        thermal_prob, 
        previous_sol_thermal,  # No previous solution
        base_params, 
        solvers
    )
    
    sol_hydraulic_prev = solve_hydraulic_problem(
        hydraulic_prob, 
        previous_sol_hydraulic,  # No previous solution
        base_params, 
        sol_thermal_prev, 
        sys_thermal, 
        solvers
    )
    
    push!(results_thermal, sol_thermal_prev)
    push!(results_hydraulic, sol_hydraulic_prev)
    
    # Solve subsequent points
    for i in 2:length(values)
        # Update parameter for this point
        current_params = [
            base_params...,
            param_sym => values[i]
        ]
        
        # Solve both systems
        sol_thermal, sol_hydraulic = simulate_gshp_system(
            thermal_prob,
            hydraulic_prob,
            sys_thermal,
            current_params,
            results_thermal[i-1],
            results_hydraulic[i-1],
            solvers
        )
        
        push!(results_thermal, sol_thermal)
        push!(results_hydraulic, sol_hydraulic)
    end
    
    return results_thermal, results_hydraulic
end

# 7. Helper function to create problems

# ╔═╡ 76cc998b-9a7f-4b1c-8403-c16bccbda74e
function prepare_problems(pars)
    # Build thermal system
    @named sys_thermal = GSHP_thermal()
    sys_thermal_simp = structural_simplify(sys_thermal)
    thermal_prob = NonlinearProblem(sys_thermal_simp, [], pars; jac = false)
    
    # Build hydraulic system
    @named sys_hydraulic = GSHP_hydraulics()
    sys_hydraulic_simp = structural_simplify(sys_hydraulic)
    hydraulic_prob = NonlinearProblem(sys_hydraulic_simp, [], pars; jac = false)
    
    return thermal_prob, hydraulic_prob, sys_thermal, sys_hydraulic
end


# ╔═╡ 5bf0b015-e0f5-4e49-baa8-c34ee5d75a73
"""
    continuation_solve(
      thermal_prob, hydraulic_prob, sys_thermal,
      base_params::Dict{Symbol,Real},
      param::Symbol,
      target_vals::AbstractVector{<:Real},
      solver;
      abstol = 1e-6
    ) -> (sol_th, sol_hy)

Start at `base_params[param] = first(target_vals)`, solve thermal+hydraulic from scratch,
then loop over the rest of `target_vals`, each time

  • setting `base_params[param] = v`  
  • calling `solve(...; p=base_params, u0=previous_sol.u, abstol=abstol)`  

and at the end returning the final `(sol_th, sol_hy)`.  
"""
function continuation_solve(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params::Dict{Symbol,Real},
    param::Symbol,
    target_vals::AbstractVector{<:Real},
    solver;
    abstol = 1e-6
)
    # 1) Kick off at the first value (solve from scratch)
    base_params[param] = first(target_vals)
    sol_th = NonlinearSolve.solve(thermal_prob, solver;
                                  p = base_params,
                                  abstol = abstol)
    sol_hy = NonlinearSolve.solve(hydraulic_prob, solver;
                                  p = base_params,
                                  abstol = abstol)

    # 2) March toward the end, re-using prev sol as u0
    for v in Iterators.rest(target_vals)
        base_params[param] = v
        sol_th = NonlinearSolve.solve(thermal_prob, solver;
                                      p  = base_params,
                                      u0 = sol_th.u,    # <-- warm start
                                      abstol = abstol)
        sol_hy = NonlinearSolve.solve(hydraulic_prob, solver;
                                      p  = base_params,
                                      u0 = sol_hy.u,    # <-- warm start
                                      abstol = abstol)
    end

    # 3) return only the terminal solutions
    return sol_th, sol_hy
end


# ╔═╡ d7090f5b-50f6-4198-9237-8447ae442789
solver_fast = FastShortcutNLLSPolyalg(autodiff=AutoFiniteDiff());

# ╔═╡ fdddfe24-f560-4a52-9756-ee1d1d6946d3
solver_robust = RobustMultiNewton(autodiff=AutoFiniteDiff());

# ╔═╡ 52e4fb9f-6e5e-4c7e-b34e-9f0bace4db63
"""
Test routine for the GSHP model
"""
function test_GSHP_stripped()
    
    @named hp = GSHP_stripped(χ=0.280);
    eqs = [
        
        hp.cycle.evaporator.b ~ 1.5e-3,
        hp.cycle.evaporator.np ~ 20,
        hp.cycle.evaporator.W ~ 0.10,
        hp.cycle.evaporator.ϕ ~ 1.18,
        hp.cycle.evaporator.β ~ π/6,
        hp.cycle.evaporator.th ~ 0.5e-3,
        hp.cycle.evaporator.λ_pl ~ 16,
        hp.cycle.evaporator.λ_pitch ~ 4e-3,
        hp.cycle.evaporator.χ ~ 0.28,
        hp.cycle.evaporator.L ~ 0.32,

		
        hp.cycle.f ~ 50,
        hp.cycle.log_ΔTSF ~ log(3.5),
        hp.cycle.log_ΔTSH ~ log(7),
        hp.cycle.log_ΔTSC ~ log(5),
        hp.cycle.logp_h ~ log(14e5),
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [hp])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])
    solver = NewtonRaphson(; 
        linesearch = BackTracking(; maxstep=20), 
        autodiff = AutoFiniteDiff()
    )
	
    sol = NonlinearSolve.solve(prob, solver_robust; abstol = 1e-6)

    return sol
end

# ╔═╡ 104492fc-1c63-4571-80e4-2df35b6abdcf
solver_fallback = NewtonRaphson(autodiff=AutoForwardDiff());

# ╔═╡ 0ab5ba9b-ee1d-4e7b-b16d-29ab59563ff2
solvers = [solver_fast, solver_robust, solver_fallback];

# ╔═╡ 58a5d982-5d07-4a63-950a-ca5639cff43f
base_params = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => 0.07,
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
    :L => 0.60,
	
    :χ => 0.28,
    :L_BHE => 150.0,
    :f_ => 50.0,
	:log_ΔTSF_ => log(3),
]

# ╔═╡ a4b78e57-3f4a-4bd1-a4fd-344aed782596
thermal_prob, hydraulic_prob, sys_thermal, sys_hydraulic = prepare_problems(base_params);

# ╔═╡ 89dc461f-71b0-4b91-ab47-a70d2929fec2
base_params_cold = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => 0.07,
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
    :L => 0.60,
	:T_wall => 273.15 + 5,
    :χ => 0.28,
    :L_BHE => 150.0,
    :f_ => 50.0,
	:log_ΔTSF_ => log(3),
]

# ╔═╡ ceb2b429-90ae-4bc5-b346-e65e5727d441
sol_start_thermal = NonlinearSolve.solve(thermal_prob, solver_robust; abstol = 1e-6)

# ╔═╡ 80898f01-7945-48d9-abd5-15de3c76d465
f_params = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => 0.07,
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
    :χ => 0.28,
    :L => 0.60,
    :L_BHE => 150.0,
    :f_ => 25.0,
	:log_ΔTSF_ => log(3),
]

# ╔═╡ b3bb8c8b-3825-4a21-94f7-a2d0d00fbb09
solve_thermal_problem(thermal_prob, sol_start_thermal, f_params, solvers)

# ╔═╡ a496bd64-5a9a-4b3f-bb8f-1a0567943d1c
solve_thermal_problem(thermal_prob, sol_start_thermal, base_params_cold, solvers)

# ╔═╡ 6ac1ecbf-4833-449f-9933-c782a8566fb0
sol_start_hydraulic = NonlinearSolve.solve(hydraulic_prob, solver_robust; abstol = 1e-6)

# ╔═╡ 8b6e0e2e-ecc5-475e-bd31-1f23893dcb38
W_vals = 0.07:0.02:0.14

# ╔═╡ 892ef213-3229-4d91-9ba7-742127f366c8
# ╠═╡ show_logs = false
results_thermal_init, results_hydraulic_init = parameter_sweep(
    thermal_prob,
    hydraulic_prob,
    sys_thermal,
    base_params,
    :W,
    W_vals,
    solvers,
	 previous_sol_thermal=sol_start_thermal, previous_sol_hydraulic=sol_start_hydraulic
)

# ╔═╡ 0ddc1362-ee46-452b-a5d7-6d03df36e31d
params_np = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => W_vals[end],
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
    :χ => 0.28,
    :L => 0.6,
    :L_BHE => 150.0,
    :f_ => 50.0
]

# ╔═╡ 30d6e371-c921-4fc2-8d8f-a6c5d7c34385
np_vals = 20:5:20

# ╔═╡ 19556fe1-b85d-4da1-ac62-a56b5a72618b
# ╠═╡ show_logs = false
results_thermal_init_np, results_hydraulic_init_np = parameter_sweep(
    thermal_prob,
    hydraulic_prob,
    sys_thermal,
    params_np,
    :np,
    np_vals,
    solvers,
	 previous_sol_thermal=results_thermal_init[end], previous_sol_hydraulic=results_hydraulic_init[end]
)

# ╔═╡ af94b9d0-44ca-42e4-9e8c-96f5864bc80d
params_b = [
    :b => 1.5e-3,
    :np => np_vals[end],
    :W => W_vals[end],
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
    :χ => 0.28,
    :L => 0.6,
    :L_BHE => 150.0,
    :f_ => 50.0
]

# ╔═╡ 8a4bc20a-195f-49f5-8379-ca6be9bdb1a1
b_vals = 1.5e-3:.25e-3:3e-3

# ╔═╡ a208968e-8f31-4a82-b8f3-9db4a8d27e8b
params_ΔT = [
    :b => b_vals[end],
    :np => np_vals[end],
    :W => W_vals[end],
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
    :χ => 0.28,
    :L => 0.6,
    :L_BHE => 150.0,
    :f_ => 50.0
]

# ╔═╡ b48efffc-d19c-4b66-936b-6fd9c9de5a88
res_ΔT = 0.1

# ╔═╡ 37d9352d-2357-4f36-b651-5dc73246fd42
ΔT_SF_vals = log.(3.5:res_ΔT:8)

# ╔═╡ 4b47374a-f1a3-402c-a4ba-bae8edcedc01
length(ΔT_SF_vals)

# ╔═╡ bfe4471c-7033-4ca5-beb1-530a8be5b757
ΔT_SF_vals_down = log.(3.5:-res_ΔT:1)

# ╔═╡ 670cd198-239a-4240-be5d-2def55383a36
# ╠═╡ show_logs = false
results_thermal_init_b, results_hydraulic_init_b = parameter_sweep(
    thermal_prob,
    hydraulic_prob,
    sys_thermal,
    params_b,
    :b,
    b_vals,
    solvers,
	 previous_sol_thermal=results_thermal_init[end], previous_sol_hydraulic=results_hydraulic_init[end]
)

# ╔═╡ 0f6b83f7-167a-4e3f-9fd5-d4f4695fbaea
# ╠═╡ show_logs = false
results_thermal, results_hydraulic = parameter_sweep(
    thermal_prob,
    hydraulic_prob,
    sys_thermal,
    params_ΔT,
    :log_ΔTSF_,
    ΔT_SF_vals,
    solvers,
	 previous_sol_thermal=results_thermal_init_b[end], previous_sol_hydraulic=results_hydraulic_init_b[end]
)

# ╔═╡ 798e29fe-16e5-4f68-baad-9edc88509245
# ╠═╡ show_logs = false
results_thermal_down, results_hydraulic_down = parameter_sweep(
    thermal_prob,
    hydraulic_prob,
    sys_thermal,
    params_ΔT,
    :log_ΔTSF_,
    ΔT_SF_vals_down,
    solvers,
	 previous_sol_thermal=results_thermal_init_b[end], previous_sol_hydraulic=results_hydraulic_init_b[end]
)

# ╔═╡ 05bbe96e-023c-413b-a347-d37e74fad7fe
res_thermal_overall = [results_thermal; results_thermal_down];

# ╔═╡ fecd7f3d-bf0b-464e-a1ad-238e508b0f38
res_hydraulic_overall = [results_hydraulic; results_hydraulic_down];

# ╔═╡ 0fe18091-e8f6-42d4-9961-5783452ddefd
#=╠═╡
function process_results(results_thermal, results_hydraulic, sys_thermal, sys_hydraulic)
    data = []  # Store rows as NamedTuples
    
    for (therm, hyd) in zip(results_thermal, results_hydraulic)
        # Thermal results
        ṁ_sf = exp(therm[sys_thermal.γṁ_sf])
        ṁ_ref = exp(therm[sys_thermal.γṁ_ref])
        T_in = therm[sys_thermal.cycle.T_i_b] .- 273.15
        T_mid = therm[sys_thermal.cycle.evaporator.T_mid_1] .- 273.15
        T_out = therm[sys_thermal.cycle.T_o_b] .- 273.15
		T_bot = therm[sys_thermal.bhe_thermal.thermal.T_b] .-273.15
        T_comp_in = therm[sys_thermal.cycle.evaporator.T_out_2] .- 273.15
        T_wall = 10.0
        #T_wall = therm[sys_thermal.bhe.T_wall] .- 273.15
		ΔT = T_in - T_out

		fr_sh = zip(therm[sys_thermal.cycle.evaporator.uz_fr_sh])
		
        T_evap = therm[sys_thermal.cycle.evaporator.sh.propref.Tsat] .- 273.15
        Q_evap = therm[sys_thermal.Q_evap]
        Q_cond = -therm[sys_thermal.Q_cond]
        W_comp = therm[sys_thermal.W_comp]
        COP_1 = therm[sys_thermal.COP1]

		U_EV = therm[sys_thermal.cycle.evaporator.ev.base.U]
		U_SH = therm[sys_thermal.cycle.evaporator.sh.base.U]
		fr_EV = (1-zap(therm[sys_thermal.cycle.evaporator.uz_fr_sh]))
		fr_SH = (zap(therm[sys_thermal.cycle.evaporator.uz_fr_sh]))
		
        # Hydraulic results
        Δp_bhe = hyd[sys_hydraulic.bhe.Δp]/1e5
        Δp_phx = hyd[sys_hydraulic.phx.base.Δp]/1e5
        Δp_tot = hyd[sys_hydraulic.pump.base.Δp]/1e5
        P_pump = hyd[sys_hydraulic.pump.base.P_elec]
        N_pump = hyd[sys_hydraulic.pump.base.N]
        # System performance
        COP_2 = Q_cond / (W_comp + P_pump)



		#temporary fix

        # Calculate Reynolds numbers from flow properties
        # Extract properties from hydraulic solution
        # For PHX
        ρ_phx = hyd[sys_hydraulic.phx.props.ρ]
        μ_phx = hyd[sys_hydraulic.phx.props.μ]
        ṁ_phx = hyd[sys_hydraulic.phx.base.ṁ]
        
        # PHX geometry (you'll need to extract these or use fixed values)
        b_phx = hyd[sys_hydraulic.phx.b]  # or use a fixed value like 1.5e-3
        W_phx = hyd[sys_hydraulic.phx.W]  # or use a fixed value like 0.07
        np_phx = hyd[sys_hydraulic.phx.np]  # or use a fixed value like 20
        ϕ_phx = hyd[sys_hydraulic.phx.ϕ]  # or use a fixed value like 1.18
        
        # Calculate PHX Reynolds number
        n_channels = (np_phx - 1) / 2
        A_flow_phx = n_channels * W_phx * (2*b_phx)
        v_phx = abs(ṁ_phx) / (ρ_phx * A_flow_phx)
        dh_phx = 2*b_phx / ϕ_phx
        Re_phx = ρ_phx * v_phx * dh_phx / μ_phx
        
        # For BHE leg 1
        ρ_bhe = hyd[sys_hydraulic.bhe.leg1.props.ρ]  # or similar path
        μ_bhe = hyd[sys_hydraulic.bhe.leg1.props.μ]
        ṁ_bhe = hyd[sys_hydraulic.bhe.ṁ]
        # BHE geometry (extract or use fixed values)
        D_i = hyd[sys_hydraulic.bhe.leg1.Di]  # inner diameter or extract from system
        A_pipe = π * (D_i/2)^2
        v_bhe = abs(ṁ_bhe) / (ρ_bhe * A_pipe)
        Re_bhe_leg1 = ρ_bhe * v_bhe * D_i / μ_bhe
        
		# For BHE leg 2
        ρ_bhe = hyd[sys_hydraulic.bhe.leg2.props.ρ]  # or similar path
        μ_bhe = hyd[sys_hydraulic.bhe.leg2.props.μ]
        ṁ_bhe = hyd[sys_hydraulic.bhe.ṁ]
        # BHE geometry (extract or use fixed values)
        D_i = hyd[sys_hydraulic.bhe.leg2.Di]  # inner diameter or extract from system
        A_pipe = π * (D_i/2)^2
        v_bhe = abs(ṁ_bhe) / (ρ_bhe * A_pipe)
        Re_bhe_leg2 = ρ_bhe * v_bhe * D_i / μ_bhe

		η_bh = 1 - (T_wall - T_in)/(T_wall - T_out)
		η_bh = therm[sys_thermal.bhe_thermal.ϵ]
		R_1b = therm[sys_thermal.bhe_thermal.utube_res.R_t1_b]
		R_2b = therm[sys_thermal.bhe_thermal.utube_res.R_t2_b]
		R_12 = therm[sys_thermal.bhe_thermal.utube_res.R_t1_t2]
		R_f1b = therm[sys_thermal.bhe_thermal.thermal.R_f1_b]
		R_f2b = therm[sys_thermal.bhe_thermal.thermal.R_f2_b]
		R_f1f2 = therm[sys_thermal.bhe_thermal.thermal.R_f1_f2]
		
        push!(data, (;
            ṁ_sf, ṁ_ref, T_wall, T_in, T_out, T_bot, T_evap, T_mid, T_comp_in, ΔT, Q_evap, Q_cond, W_comp, Δp_tot, Δp_phx, Δp_bhe, P_pump, N_pump, Re_bhe_leg1, Re_bhe_leg2, Re_phx, COP_1, COP_2, U_EV, U_SH, fr_EV, fr_SH, η_bh, R_1b, R_2b, R_12, R_f1b, R_f2b, R_f1f2
        ))
    end
	df = sort(DataFrame(data), :ΔT)
	max_Q = maximum(df.Q_cond)
	df.P_aux_heating = max_Q .- df.Q_cond
	df.COP_3 = (df.Q_cond .+ df.P_aux_heating) ./ (df.W_comp .+ df.P_pump .+ df.P_aux_heating)
    
    return df  # Convert to DataFrame
end
  ╠═╡ =#

# ╔═╡ 7d25aec4-c657-410a-8ae4-03195cefbee7
#=╠═╡
"""
Create a series of smoothed plots with optional y-axis zeroing.

# Arguments
- `df::DataFrame`: Input data with columns `:ṁ`, `:T_in`, `:T_out`, `:Q_evap`, `:W_comp`, `:P_pump`, `:Δp_tot`, `:COP_thermal`, `:COP_overall`.
- `smooth_degree::Int`: Degree of smoothing (LOESS).
- `y_zero::Bool`: If `true`, set the lower y-limit of COP and Δp plots to 0.

# Returns
A vector of four CairoMakie Figure objects in order:
1. Brine Temperature Profile
2. System Efficiency Comparison
3. Energy Flow Components
4. Hydraulic Resistance
"""
function create_plots(df_::DataFrame; smooth_degree::Int=1, y_zero::Bool=false)
    set_aog_theme!()

    df = transform(df_,
        :T_in  => ByRow(x -> x - 273.15) => :T_in_C,
        :T_out => ByRow(x -> x - 273.15) => :T_out_C
    )
    set_aog_theme!()

    # 1. Temperature Profile
    plt_temp = data(df) * mapping(
        :ṁ,
        :T_in_C => "Inlet Temperature (°C)",
        color = :T_out_C => "Outlet Temperature (°C)"
    ) * visual(Lines)
    temp_fig = draw(plt_temp; axis=(
        xlabel = L"Mass Flow Rate $\dot{m}$ (kg/s)",
        ylabel = "Temperature (°C)",
        title = "Brine Temperature Profile"
    ))

    # 2. COP Comparison
    layer_th = data(df) * mapping(:ṁ, :COP_1 => "Thermal COP")# * smooth(; degree=smooth_degree)
    layer_ov = data(df) * mapping(:ṁ, :COP_2 => "Overall COP")# * smooth(; degree=smooth_degree)
    plt_cop = layer_th + layer_ov
    cop_fig = draw(plt_cop; axis=(
        xlabel = L"Mass Flow Rate $\dot{m}$ (kg/s)",
        ylabel = "Coefficient of Performance",
        title = "System Efficiency Comparison",
        limits = (y = (y_zero ? 0 : nothing, nothing))
    ))

    # 3. Power Components
    power_long = stack(df, [:Q_evap, :Q_cond, :W_comp, :P_pump], variable_name=:component, value_name=:power)
    plt_power = data(power_long) * mapping(
        :ṁ,
        :power => "Power (W)",
        color = :component => nonnumeric
    )# * smooth(; degree=smooth_degree)
    power_fig = draw(plt_power; axis=(
        xlabel = L"Mass Flow Rate $\dot{m}$ (kg/s)",
        ylabel = "Power (W)",
        title = "Energy Flow Components"
    ))

    # 4. Pressure Drop
    plt_pressure = data(df) * mapping(:ṁ, :Δp_tot => "Pressure Drop (Pa)")# * smooth(; degree=smooth_degree)
    pressure_fig = draw(plt_pressure; axis=(
        xlabel = L"Mass Flow Rate $\dot{m}$ (kg/s)",
        ylabel = L"Pressure Drop $\Delta p$ (Pa)",
        title = "Hydraulic Resistance",
        limits = (y = (y_zero ? 0 : nothing, nothing))
    ))

    return [temp_fig, cop_fig, power_fig, pressure_fig]
end
  ╠═╡ =#

# ╔═╡ f100b275-2eb2-4d10-93e2-8954f314ae9e

# Merge extra style fields into your discovered transitions
_apply_transition_styles(base::Dict{String,NamedTuple},
                         styles::Dict{String,NamedTuple}) =
    Dict(k => (; v..., get(styles, k, (;))...) for (k, v) in base)


# ╔═╡ 058c533c-104d-4ae9-80e3-f356787d0d13
#=╠═╡
# Compute y(xt) by 1D linear interpolation on (xcol, ycol) from a *wide* df
# Assumes df is sorted by xcol (you already do sort!(df2, :ΔT))
function _y_at(df::DataFrame, xcol::Symbol, ycol::Symbol, xt::Real)
    x = df[!, xcol] |> collect
    y = df[!, ycol] |> collect
    itp = interpolate((x,), y, Gridded(Linear()))
    return itp(xt)
end

# Draw vertical guides and put text on a chosen curve at the transition ΔT
  ╠═╡ =#

# ╔═╡ 4899159f-1c86-478c-b874-ef8db1808bc5
# Accept Makie.Axis, AoG AxisEntries, or arrays and return the Makie.Axis
_as_axis(axlike) =
    axlike isa Makie.Axis                ? axlike :
    (hasproperty(axlike, :axis)          ? getproperty(axlike, :axis) :
     axlike isa AbstractArray            ? _as_axis(first(axlike)) :
     error("Couldn't extract Makie.Axis from $(typeof(axlike))"))

# Auto left/right alignment & small x offset near plot edges

# ╔═╡ c8044a03-837b-411c-8919-7d47f83111c3
# robust axis ranges across Makie versions
function _axis_xrange(ax::Makie.Axis)
    fr = ax.finallimits[]; (fr.origin[1], fr.origin[1] + fr.widths[1])
end

# ╔═╡ 5575055e-9136-4b4e-925a-ddf82ca1156c
# Horizontal auto-align, but force vertical anchor at :bottom / :top
_auto_align_offset_bottom(ax::Makie.Axis, x; margin=0.92, xpad=6, ypad=6) = begin
    xmin, xmax = _axis_xrange(ax)
    if !(isfinite(xmin) && isfinite(xmax)) || xmin == xmax
        ((:left, :bottom), (xpad, ypad))
    elseif x > xmin + margin * (xmax - xmin)
        ((:right, :bottom), (-xpad, ypad))
    else
        ((:left,  :bottom), ( xpad, ypad))
    end
end

# ╔═╡ f8b2fee6-eb35-493c-ac30-331176d9a795
_auto_align_offset_top(ax::Makie.Axis, x; margin=0.92, xpad=6, ypad=6) = begin
    xmin, xmax = _axis_xrange(ax)
    if !(isfinite(xmin) && isfinite(xmax)) || xmin == xmax
        ((:left, :top), (xpad, -ypad))
    elseif x > xmin + margin * (xmax - xmin)
        ((:right, :top), (-xpad, -ypad))
    else
        ((:left,  :top), ( xpad, -ypad))
    end
end

# ╔═╡ b220d228-e3cb-4dd2-843e-03db76f46fe9
function _axis_yrange(ax::Makie.Axis)
    fr = ax.finallimits[]; (fr.origin[2], fr.origin[2] + fr.widths[2])
end

# ╔═╡ 4412cbfc-d470-4003-a2bd-606df1be7022
_px_to_data(ax::Makie.Axis; dx_px::Real=0, dy_px::Real=0) = begin
    xmin, xmax = _axis_xrange(ax); ymin, ymax = _axis_yrange(ax)
    Δx = xmax - xmin; Δy = ymax - ymin
    pa = ax.scene.px_area[]
    sx = Δx / pa.widths[1]
    sy = Δy / pa.widths[2]
    (dx_px * sx, dy_px * sy)
end

# ╔═╡ d7384891-c5c0-4614-99cf-0f1cf58f759c


# auto left/right alignment based on where x sits in the axis range
function _auto_align_offset(ax::Makie.Axis, x; margin=0.92, pad=6)
    xmin, xmax = _axis_xrange(ax)
    if !(isfinite(xmin) && isfinite(xmax)) || xmin == xmax
        return ((:left, :center), (pad, 0))  # conservative fallback
    end
    if x > xmin + margin * (xmax - xmin)
        ((:right, :center), (-pad, 0))
    else
        ((:left,  :center), ( pad, 0))
    end
end


# ╔═╡ 9b2f80c6-18d0-4b10-827a-b9a665f7f1c4
#=╠═╡
function _label_transitions!(ax::Makie.Axis, df::DataFrame;
        xcol::Symbol=:ΔT, ycol::Union{Nothing,Symbol}=nothing,
        transitions,                               # Dict or Vector of NamedTuples
        where::Symbol = :bottom,                   # default placement (per-transition override allowed)
        rotation_deg::Real = 90,
        fontsize::Real = 12, textcolor=:gray10,    # default text style
        vcolor=:gray40, vstyle=:dash, vwidth=1,    # default line style
        # data-space padding fractions
        x_data_pad_frac::Real = 0.015,
        y_data_pad_frac::Real = 0.02,
        # extra pixel pads to keep labels clear of axes regardless of scale
        x_px_pad::Real = 6, y_px_pad::Real = 10,
        # switch from right→left alignment near the right edge
        edge_margin::Real = 0.92
    )

    angle = rotation_deg * π / 180
    xmin, xmax = _axis_xrange(ax); ymin, ymax = _axis_yrange(ax)
    Δx = xmax - xmin; Δy = ymax - ymin

    # base data-space pads
    xpad_d = x_data_pad_frac * Δx
    ypad_d = y_data_pad_frac * Δy

    # pixel→data pads (robust lift so text never hides behind x-axis)
    xpad_px_d, ypad_px_d = _px_to_data(ax; dx_px=x_px_pad, dy_px=y_px_pad)

    iter = transitions isa AbstractDict ? values(transitions) : transitions
    for t in iter
        # required fields
        xt    = t[:dT]
        label = t[:label]

        # per-transition overrides (fall back to the call-site defaults)
        t_where        = haskey(t, :where)         ? t[:where]         : where
        t_rot          = haskey(t, :rotation_deg)  ? t[:rotation_deg]  : rotation_deg
        t_size         = haskey(t, :fontsize)      ? t[:fontsize]      : fontsize
        t_tcolor       = haskey(t, :textcolor)     ? t[:textcolor]     : textcolor
        t_vcolor       = haskey(t, :vcolor)        ? t[:vcolor]        : vcolor
        t_vstyle       = haskey(t, :vstyle)        ? t[:vstyle]        : vstyle
        t_vwidth       = haskey(t, :vwidth)        ? t[:vwidth]        : vwidth
        t_xpad_frac    = haskey(t, :x_data_pad_frac) ? t[:x_data_pad_frac] : x_data_pad_frac
        t_ypad_frac    = haskey(t, :y_data_pad_frac) ? t[:y_data_pad_frac] : y_data_pad_frac
        t_xpx          = haskey(t, :x_px_pad)      ? t[:x_px_pad]      : x_px_pad
        t_ypx          = haskey(t, :y_px_pad)      ? t[:y_px_pad]      : y_px_pad
        t_z_text       = haskey(t, :z_text)        ? t[:z_text]        : 100
        t_z_line       = haskey(t, :z_line)        ? t[:z_line]        : 0

        # recompute pads if per-transition changed
        xpad_d_t = t_xpad_frac * Δx + (_px_to_data(ax; dx_px=t_xpx, dy_px=0)[1])
        ypad_d_t = t_ypad_frac * Δy + (_px_to_data(ax; dx_px=0, dy_px=t_ypx)[2])

        # y placement
        yt = if t_where === :curve && ycol !== nothing
            _y_at(df, xcol, ycol, xt) + ypad_d_t
        elseif t_where === :top
            ymax - ypad_d_t
        else
            ymin + ypad_d_t               # <- lifted enough to clear the x-axis
        end

        # side & alignment
        if xt > xmin + edge_margin*Δx
            xlbl  = xt + xpad_d_t
            align = ( :left,  t_where === :top ? :top : :bottom )
        else
            xlbl  = xt - xpad_d_t
            align = ( :right, t_where === :top ? :top : :bottom )
        end
		
		# 1) Guide line: keep using vlines! (no z)
		vlines!(ax, [xt]; color=t_vcolor, linestyle=t_vstyle, linewidth=t_vwidth)
		
		# 2) Label: drop z=..., use overdraw=true so it draws above decorations
		text!(ax, label;
		      position=(xlbl, yt),
		      align=align,
		      rotation=t_rot * π / 180,
		      fontsize=t_size,
		      color=t_tcolor,
		      space=:data,
		      overdraw=true)   # <- instead of z=t_z_text

    end
    ax
end

  ╠═╡ =#

# ╔═╡ 86429d81-3516-488f-be3d-df2d83bdbae8

# helper to make filesystem-friendly filenames
_slug(s::AbstractString) = begin
    t = lowercase(String(s))
    t = replace(t, r"[^\p{L}\p{N}]+" => "_")
    t = replace(t, r"_+" => "_")
    strip(t, '_')
end

# ╔═╡ 4982df2c-016a-4136-a074-b39c1befffb7
#=╠═╡

function create_plots_by_dT(df::DataFrame;
    smooth_degree::Int=1,
    y_zero::Bool=false,
    re_thresholds::Dict{String, Tuple{Float64, Float64}} = Dict(
        "Re_bhe_leg1" => (2300, 4000),
        "Re_bhe_leg2" => (2300, 4000),
        "Re_phx"      => (100, 2000)
    ),
    # saving options
    save_dir::Union{Nothing, AbstractString}=nothing,
    save_prefix::AbstractString = "",
    save_formats::Vector{Symbol} = [:svg],
)
    df2 = copy(df)
    set_aog_theme!()
    sort!(df2, :ΔT)

    # --- transitions (unchanged) ---
    transition_points = Dict{String, NamedTuple}()
    re_cols   = ["Re_bhe_leg1", "Re_bhe_leg2", "Re_phx"]
    re_labels = Dict(
        "Re_bhe_leg1" => "BHE Leg1",
        "Re_bhe_leg2" => "BHE Leg2",
        "Re_phx"      => "Evaporator"
    )
    for col in re_cols
        col ∈ names(df2) || continue   # now calls DataFrames.names
        x = df2.ΔT; y = df2[!, col]
        length(x) < 2 && error("Need ≥2 points for interpolation in $col")
        itp = interpolate((x,), y, Gridded(Linear()))
        for thresh in get(re_thresholds, col, ())
            f(dt) = itp(dt) - thresh
            idx = findfirst(i -> f(x[i]) * f(x[i+1]) ≤ 0, 1:length(x)-1)
            isnothing(idx) && continue
            a, b = x[idx], x[idx+1]
            root = find_zero(f, (a, b))
            key = "$col-$(thresh)"
            transition_points[key] = (dT = root,
                                      label = "$(re_labels[col]) Re=$(thresh)",
                                      threshold = thresh)
        end
    end

	# After you compute `transition_points` in create_plots_by_dT:
	transition_styles = Dict(
	    "Re_bhe_leg1-2300" => (
	        label    = "Leg 1: Re=2300 (L→T)",
	        vcolor   = :crimson, vstyle=:dash,  vwidth=2.0,
	        textcolor=:crimson, fontsize=11, where=:bottom, y_px_pad=14,
	    ),
	    "Re_bhe_leg1-4000" => (
	        label    = "Leg 1: Re=4000 (T→Tu)",
	        vcolor   = :crimson, vstyle=:dot,   vwidth=2.0,
	        textcolor=:crimson, fontsize=11, where=:bottom, y_px_pad=14,
	    ),
	    "Re_bhe_leg2-2300" => (
	        label    = "Leg 2: Re=2300",
	        vcolor   = :blue,    vstyle=:dashdot, vwidth=2.0,
	        textcolor=:blue,     where=:top, y_px_pad=12,
	    ),
	    "Re_bhe_leg2-4000" => (
	        label    = "Leg 2: Re=4000",
	        vcolor   = :blue,    vstyle=:dashdot, vwidth=2.0,
	        textcolor=:blue,     where=:top, y_px_pad=12,
	    ),
	)
	
	transition_points = _apply_transition_styles(transition_points, transition_styles)

	
    figures   = Vector{Figure}()
    fig_names = String[]   # <-- renamed (no shadowing)


    # 3) Temperature vs ΔT ------------------------------------------------------
    temp_long = stack(df2, [:T_wall, :T_in, :T_mid, :T_out, :T_bot, :T_evap, :T_comp_in],
                      variable_name=:temp_type, value_name=:temperature)
    map_temp = Dict(
        "T_wall"=>"Wall Temperature", "T_in"=>"BHE -> PHX Temperature",
        "T_mid"=>"Mid Temperature", "T_bot"=>"BHE bottom Temperature",
        "T_out"=>"PHX -> BHE Temperature", "T_evap"=>"Evaporation Temperature",
        "T_comp_in"=>"Inlet Compressor Temperature"
    )
    temp_long.temp_type .= map(s -> map_temp[s], temp_long.temp_type)

    temp_plt = data(temp_long) *
               mapping(:ΔT, :temperature, color=:temp_type=>nonnumeric) *
               visual(Lines)
    fig_temp, ax_temp = draw(temp_plt; axis=(; title="Temperature vs ΔT"),
                                        legend=(position=:right,))
    _label_transitions!(_as_axis(ax_temp), df2; ycol=:T_out,
        transitions=transition_points, textcolor=:gray10, fontsize=12)
    push!(figures, fig_temp); push!(fig_names, "temperature_vs_dT")

    # 4) COP vs ΔT --------------------------------------------------------------
    cop_long = stack(df2, [:COP_1, :COP_2, :COP_3],
                     variable_name=:cop_type, value_name=:cop_value)
    map_cop = Dict("COP_1"=>"COP 1", "COP_2"=>"COP 2", "COP_3"=>"COP 3")
    cop_long.cop_type .= map(s -> map_cop[s], cop_long.cop_type)

    base_cop = data(cop_long) *
               mapping(:ΔT=>"ΔT (K)", :cop_value=>"COP", color=:cop_type=>nonnumeric) *
               visual(Lines, linewidth=2)
    fig_cop, ax_cop = draw(base_cop;
        axis=(xlabel="ΔT (K)", ylabel="COP", title="COP vs ΔT",
              limits=(y=(y_zero ? 0 : nothing, nothing))),
        legend=(position=:right,))
    _label_transitions!(_as_axis(ax_cop), df2; ycol=:COP_3,
        transitions=transition_points)
    push!(figures, fig_cop);  push!(fig_names, "cop_vs_dT")

    # 5) Power vs ΔT ------------------------------------------------------------
    power_long = stack(df2, [:Q_evap, :Q_cond, :W_comp, :P_pump, :P_aux_heating],
                       variable_name=:component, value_name=:power)
    base_power = data(power_long) *
                 mapping(:ΔT=>"ΔT (K)", :power=>"Power (W)", color=:component=>nonnumeric) *
                 visual(Lines)
    fig_power, ax_power = draw(base_power;
        axis=(xlabel="ΔT (K)", ylabel="Power (W)", title="Energy Components vs ΔT"),
        legend=(position=:right,))
    _label_transitions!(_as_axis(ax_power), df2; ycol=:Q_cond,
        transitions=transition_points)
    push!(figures, fig_power); push!(fig_names, "energy_components_vs_dT")

    # 6) Pressure Drop vs ΔT ----------------------------------------------------
	base_pres = data(df2) *
    mapping(
        :ΔT      => "ΔT (K)",
        :Δp_tot  => (x -> x .* 1e5) => "Pressure (Pa)"   # bar → Pa
    ) *
    visual(Lines)

	
	fig_pres, ax_pres = draw(base_pres;
	    axis=(; xlabel="ΔT (K)",
	            ylabel=L"Pressure\ (Pa)",
	            title="Hydraulic Resistance vs ΔT",
	            yscale = log10,                 # <— log10 axis
	            limits=(y=(y_zero ? 1.0 : nothing, nothing))),  # avoid 0 on log scale
	    legend=(position=:right,))
	
	_label_transitions!(_as_axis(ax_pres), df2; transitions=transition_points, where=:bottom)
	push!(figures, fig_pres);  push!(fig_names, "pressure_drop_vs_dT")



    # 7) U-values vs ΔT ---------------------------------------------------------
    u_long = stack(df2, [:U_EV, :U_SH], variable_name=:u_type, value_name=:u_value)
    map_u = Dict("U_EV"=>"U Evaporator", "U_SH"=>"U Superheat")
    u_long.u_type .= map(s->map_u[s], u_long.u_type)
    u_layer = data(u_long) *
              mapping(:ΔT=>"ΔT (K)", :u_value=>"U (W/m²K)", color=:u_type=>nonnumeric) *
              visual(Lines)
    fig_u, ax_u = draw(u_layer;
        axis=(xlabel="ΔT (K)", ylabel="U (W/m²K)", title="Heat Transfer Coefficients vs ΔT"),
        legend=(position=:right,))
    _label_transitions!(_as_axis(ax_u), df2; ycol=:U_EV,
        transitions=transition_points)
    push!(figures, fig_u);     push!(fig_names, "u_values_vs_dT")

    # 8) Section Lengths vs ΔT --------------------------------------------------
    l_long = stack(df2, [:fr_EV, :fr_SH], variable_name=:l_type, value_name=:length)
    map_l = Dict("fr_EV"=>"Evaporator Length", "fr_SH"=>"Superheat Length")
    l_long.l_type .= map(s->map_l[s], l_long.l_type)
    l_layer = data(l_long) *
              mapping(:ΔT=>"ΔT (K)", :length=>"Length (m)", color=:l_type=>nonnumeric) *
              visual(Lines)
    fig_l1, ax_l1 = draw(l_layer;
        axis=(xlabel="ΔT (K)", ylabel="Length (m)", title="Section Lengths vs ΔT"),
        legend=(position=:right,))
    _label_transitions!(_as_axis(ax_l1), df2; ycol=:fr_EV,
        transitions=transition_points)
    push!(figures, fig_l1);    push!(fig_names, "section_lengths_vs_dT")

    # 9) Efficiency vs ΔT -------------------------------------------------------
    l2_long = stack(df2, [:η_bh], variable_name=:l_type, value_name=:length)
    map_l2 = Dict("η_bh"=>"BH efficiency")
    l2_long.l_type .= map(s->map_l2[s], l2_long.l_type)
    l2_layer = data(l2_long) *
               mapping(:ΔT=>"ΔT (K)", :length=>"Length (m)", color=:l_type=>nonnumeric) *
               visual(Lines)
    fig_l2, ax_l2 = draw(l2_layer;
        axis=(xlabel="ΔT (K)", ylabel="Length (m)", title="Section Lengths vs ΔT"),
        legend=(position=:right,))
    _label_transitions!(_as_axis(ax_l2), df2; ycol=:η_bh,
        transitions=transition_points)
    push!(figures, fig_l2);    push!(fig_names, "bh_efficiency_vs_dT")

    # 10) Thermal resistances vs ΔT --------------------------------------------
    r_long = stack(df2, [:R_1b, :R_2b, :R_12, :R_f1b, :R_f2b, :R_f1f2],
                   variable_name=:r_type, value_name=:r_val)
    map_r = Dict("R_1b"=>"R pipe 1 wall", "R_2b"=>"R pipe 2 wall",
                 "R_12"=>"R pipe 1 ↔ 2", "R_f1b"=>"R fluid 1 ↔ wall",
                 "R_f2b"=>"R fluid 2 ↔ wall", "R_f1f2"=>"R fluid 1 ↔ 2")
    r_long.r_type .= map(s->map_r[s], r_long.r_type)
    r_layer = data(r_long) *
              mapping(:ΔT=>"ΔT (K)", :r_val=>"R (K/W)", color=:r_type=>nonnumeric) *
              visual(Lines)
    fig_r, ax_r = draw(r_layer;
        axis=(xlabel="ΔT (K)", ylabel="R (K/W)", title="Thermal Resistances vs ΔT"),
        legend=(position=:right,))
    # anchor on one available resistance for y (fallback if missing)
    anchor_y = :R_1b ∈ names(df2) ? :R_1b : first(Symbol.(names(df2, r"^R_")))
    _label_transitions!(_as_axis(ax_r), df2; ycol=anchor_y,
        transitions=transition_points)
    push!(figures, fig_r);     push!(fig_names, "thermal_resistances_vs_dT")

    # --- save if requested ---
    if save_dir !== nothing
        mkpath(save_dir)
        for (fig, nm) in zip(figures, fig_names)
            base = isempty(save_prefix) ? _slug(nm) : "$( _slug(save_prefix) )_$( _slug(nm) )"
            for fmt in save_formats
                path = joinpath(save_dir, "$(base).$(String(fmt))")
                CairoMakie.save(path, fig)
            end
        end
    end


    return figures
end
  ╠═╡ =#

# ╔═╡ c18c39fe-6896-433e-90d4-4a8149ae9b5b
#=╠═╡
df = process_results(res_thermal_overall, res_hydraulic_overall, sys_thermal, sys_hydraulic);
  ╠═╡ =#

# ╔═╡ ce0ce754-c26a-4b8d-81d9-8234a165a0c6
new_thresh = Dict("Re_bhe_leg1"=>(2300.,4000.), "Re_bhe_leg2"=>(2300.,4000.))

# ╔═╡ bed552fa-1a51-492b-8552-81186c96c6a1
# ╠═╡ show_logs = false
#=╠═╡
plot_A_ΔT, plot_B_ΔT, plot_C_ΔT, plot_D_ΔT, plot_E_ΔT, plot_F_ΔT, plot_G_ΔT , plot_H_ΔT = create_plots_by_dT(df; re_thresholds=new_thresh)
  ╠═╡ =#

# ╔═╡ 473f4d2b-92a3-495c-975d-38738f714fb9
#=╠═╡
plot_A_ΔT
  ╠═╡ =#

# ╔═╡ c2755bcf-548b-47ff-be58-b6cbfda319f7
#=╠═╡
plot_B_ΔT
  ╠═╡ =#

# ╔═╡ afba0093-bdb2-4ee9-95da-0f4d1933d22c
#=╠═╡
plot_C_ΔT
  ╠═╡ =#

# ╔═╡ bdbf26be-8a8b-4c4d-be59-a67ea68226ef
#=╠═╡
plot_D_ΔT
  ╠═╡ =#

# ╔═╡ 85978ad5-c80a-4eed-b04f-72b57afcc13e
#=╠═╡
plot_E_ΔT
  ╠═╡ =#

# ╔═╡ 2ce9d158-4156-47db-97bc-efa0e18618b3
#=╠═╡
plot_F_ΔT
  ╠═╡ =#

# ╔═╡ aa8498ba-14aa-421d-abb8-67ac1286bd41
#=╠═╡
plot_G_ΔT
  ╠═╡ =#

# ╔═╡ 0761ab79-e4d9-43c2-b772-fe6df8f333b7
#=╠═╡
plot_H_ΔT
  ╠═╡ =#

# ╔═╡ b92cc8df-5884-4129-99fe-d22038b4eb0c
"""
    param_pair_range(sym::Symbol, a::Real, b::Real, N::Integer)
Return a Vector of N linearly-spaced `sym => value` pairs from `a` to `b`.
"""
function param_pair_range(sym::Symbol, a::Real, b::Real, N::Integer)
    values = range(a, stop=b, length=N)
    return [ sym => v for v in values ]
end


# ╔═╡ 1e116b39-5019-4e87-80eb-13b1a9d9ae58
# ╠═╡ show_logs = false
thermal_prob2, hydraulic_prob2, sys_thermal2, sys_hydraulic2 = prepare_problems(base_params);

# ╔═╡ fb29ede8-6f49-4e14-a6d9-d8caee90e5ba
Npoints = 45

# ╔═╡ 870d81bf-2254-427b-8471-3c9d23b13a72
f_values = collect(range(50, 10, length = Npoints))

# ╔═╡ 24b94de6-1d03-4225-8fdd-95cac5beb33c
250+708+7634-1400

# ╔═╡ 977b3df3-0fa4-452f-884b-94b24ec54cf3
# high-contrast, CVD-safe accents (Okabe–Ito-ish)
const OI_VERMILLION = RGBf(0.835, 0.365, 0.000)  # #D55E00

# ╔═╡ 1320fac5-2cf2-484d-a01e-98ecd559588d
const OI_MAGENTA    = RGBf(0.800, 0.475, 0.655)  # #CC79A7

# ╔═╡ a3fe0b30-007d-415d-af04-f5712e3d0489
const CHARCOAL      = RGBf(0.100, 0.100, 0.100)

# ╔═╡ 298f0a4b-d6e7-4772-9862-54895ba7aab8


# ╔═╡ ea33f723-8f17-4000-b483-eb8f1f9354f3
md" # Comparison BH efficiency"

# ╔═╡ dc7866d3-68e8-40de-b82d-26717fff53d9
md" # Comparison T outlet BH"

# ╔═╡ 21a3e541-8ad9-41de-995a-e2b43221d15a
md" # Comparison T evap"

# ╔═╡ 712bf240-09c1-4afb-be0d-5d76f9e5bbb1
md" # Comparison COP3"

# ╔═╡ 5ce8337c-aae0-476f-914b-854a054b26db
md" # Comparison ridge position"

# ╔═╡ a966eea4-f3e4-46f9-8230-7ab19eb7f905


# ╔═╡ e2504f4a-6606-4051-b91e-189bbd218edf


# ╔═╡ fdfb7456-6d8e-44c8-a883-1d1faeb96ede
md" # Comparison pump efficiency"

# ╔═╡ 1a9a9e58-41ac-45fb-9891-2e2e531621e5


# ╔═╡ 8c67c3c7-1452-42cf-a7f6-0a40b5bdd071
md" # Other explorations"

# ╔═╡ 4e1204d9-eab8-4890-a957-92341c8edfc4
md" # 2D Parametric Simulations"

# ╔═╡ 987290e0-8b32-481a-85a0-442b075463cb
md" ## Base"

# ╔═╡ 8ddb2096-a2b0-4369-a188-5d972af91428
md" ## Variant with colder wall"

# ╔═╡ 44e95eb8-ef96-4832-9168-df5c504b0569
md" ## Variant with water"

# ╔═╡ 8a1007ab-a556-4670-a2a6-18292d721c89
base_params_water = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => 0.07,
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
	:λ_grout => 1.5,
	:λ_soil => 3.5,
    :χ => 0.15,
    :L => 0.60,
    :L_BHE => 150.0,
    :f_ => 50.0,
	:log_ΔTSF_ => log(3),
]

# ╔═╡ 0d69490f-83ee-46fc-b034-8370e8a66361
md" ## Variant with bigger pipes (40mm OD)"

# ╔═╡ 6d9c073c-e05b-45ec-a4b4-1b47c43cc5ba
base_params_pipe = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => 0.07,
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
	:λ_grout => 1.5,
	:λ_soil => 3.5,
    :χ => 0.15,
    :L => 0.60,
    :L_BHE => 150.0,
    :f_ => 50.0,
	:log_ΔTSF_ => log(3),
	:D_p1_o => 40e-3,
	:D_p1_i => 28.24e-3/32e-3*40e-3,
	:D_p2_o => 40e-3,
	:D_p2_i => 28.24e-3/32e-3*40e-3,
	:D_b => 0.115,

		
]

# ╔═╡ 6b0ca393-a200-4ee0-b71e-ea3f3b68c3f4
md" ## Variant with longer BH (40mm OD)"

# ╔═╡ 880e6614-2556-4d3a-9408-f59b45ba94fe
base_params_long = [
    :b => 1.5e-3,
    :np => 20.0,
    :W => 0.07,
    :ϕ => 1.18,
    :β => π/3,
    :th => 0.5e-3,
    :λ_pl => 16.0,
    :λ_pitch => 4e-3,
	:λ_grout => 1.5,
	:λ_soil => 3.5,
    :χ => 0.15,
    :L => 0.60,
    :L_BHE => 200.0,
    :f_ => 50.0,
	:log_ΔTSF_ => log(3),
	:D_b => 0.115,

		
]

# ╔═╡ c5875432-f3f0-4495-a0f0-853d87d73353
md" ## Looking at the wall temperature"

# ╔═╡ acbfba23-bbb5-4d9e-93b7-1120677332bd
log_ΔTSF_values = log.(collect(range(0.5, 8.0, length = Npoints)))

# ╔═╡ b0b5c52c-a795-40a0-9c8d-f678010c454a
log_msf_values = log.(collect(range(0.1, 2.0, length = Npoints)))

# ╔═╡ 08d88b40-980e-4d17-91d1-2b1a06180bdd
T_wall_values = 273.15 .+ collect(range(5, 10, length = Npoints))

# ╔═╡ f0e43670-bbd9-487a-b8b6-8b5d6804039c
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
sols_thermal_wall, sols_hydraulic_wall = multithreaded_grid_sweep(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params, :log_ΔTSF_, log_ΔTSF_values, :T_wall, T_wall_values, solvers;
    previous_sol_thermal=nothing,
    previous_sol_hydraulic=nothing
)
  ╠═╡ =#

# ╔═╡ d7d4f78c-a845-42ae-86ac-1fbca7c3cfcb
#=╠═╡
df_grid_wall = process_grid_results(
    log_ΔTSF_values, T_wall_values,
    sols_thermal_wall, sols_hydraulic_wall,
    sys_thermal2, sys_hydraulic2;
    name1 = :log_ΔTSF_,
    name2 = :T_wall,
	target_Q_cond = 5000
);
  ╠═╡ =#

# ╔═╡ 09952ea8-ef5d-414f-ac15-478ad8ea6343
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:η_bh, levels = 0.0:0.05:1.0,
    overlays=[
        (z=:N_pump,       levels=83 .*(0.3:0.1:1),       color=:black, linewidth=1,   labels=false),
        (z=:Re_bhe_leg1,  levels=[2300, 4000],     color=:red,   linewidth=1.2, linestyle=:dash, labels=false),
        (z=:Re_bhe_leg2,  levels=[2300, 4000],     color=:blue,  linewidth=1.2, linestyle=:dash, labels=false),
    ]
)
  ╠═╡ =#

# ╔═╡ 16868428-1547-4b8d-a900-5e5120f38088
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:fr_EV, levels = 0:0.1:1.0,
    overlays=[
        (z=:N_pump,       levels=83 .*(0.3:0.1:1),       color=:black, linewidth=1,   labels=true),
        (z=:Re_bhe_leg1,  levels=[2300, 4000],     color=:red,   linewidth=1.2, linestyle=:dash, labels=true),
        (z=:Re_bhe_leg2,  levels=[2300, 4000],     color=:blue,  linewidth=1.2, linestyle=:dash, labels=true),
    ]
)
  ╠═╡ =#

# ╔═╡ 62b6a403-53ac-42c1-9fff-055a9b7a478b
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:COP_1,
    levels=2.0:0.1:4.0,
    overlays=[
        (z=:N_pump,      levels=83 .* (0.3:0.1:1.0), color=:black, linewidth=1,   labels=true),
        (z=:Re_bhe_leg1, levels=[2300, 4000],         color=:red,   linewidth=1.2, linestyle=:dash, labels=true),
        (z=:Re_bhe_leg2, levels=[2300, 4000],         color=:blue,  linewidth=1.2, linestyle=:dash, labels=true),
    ]
)

  ╠═╡ =#

# ╔═╡ 849d55d9-3e76-4ccb-8635-c8b70376d0b1
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:COP_2,
    levels=2.0:0.1:4.0,
    overlays=[
        (z=:N_pump,      levels=83 .* (0.3:0.1:1.0), color=:black, linewidth=1,   labels=true),
        (z=:Re_bhe_leg1, levels=[2300, 4000],         color=:red,   linewidth=1.2, linestyle=:dash, labels=true),
        (z=:Re_bhe_leg2, levels=[2300, 4000],         color=:blue,  linewidth=1.2, linestyle=:dash, labels=true),
    ]
)

  ╠═╡ =#

# ╔═╡ 34ed86a4-20e9-4b2b-9344-7ec6dfaab815
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:COP_3,
    levels=2.0:0.1:4.0,
    overlays=[
        (z=:N_pump,      levels=83 .* (0.3:0.1:1.0), color=:black, linewidth=1,   labels=true),
        (z=:Re_bhe_leg1, levels=[2300, 4000],         color=:red,   linewidth=1.2, linestyle=:dash, labels=true),
        (z=:Re_bhe_leg2, levels=[2300, 4000],         color=:blue,  linewidth=1.2, linestyle=:dash, labels=true),
    ]
)

  ╠═╡ =#

# ╔═╡ 41fe5c48-3b28-4d41-8969-e3d59a84cbfa
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:Q_cond,
    overlays=ovs_plots_2
)

  ╠═╡ =#

# ╔═╡ 54809c43-f7c3-493b-b06a-a13bf0a44bc6

ovs_plots_2 = [
        (z=:N_pump,color=:black, linewidth=1,   labels=true),
        (z=:Re_bhe_leg1, levels=[2300, 4000],         color=:red,   linewidth=1.2, linestyle=:dash, labels=true),
        (z=:Re_bhe_leg2, levels=[2300, 4000],         color=:blue,  linewidth=1.2, linestyle=:dash, labels=true),
    ]


# ╔═╡ bf3a0c4a-5fea-465d-a3ca-e742b74dc78f
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:Δp_tot,
    levels=range(0.1^0.5, 5^0.5, length = 15).^2, xlabel = "ΔT Secondary Fluid [K]", ylabel = "Borehole wall temperature [K]", zlabel = "Pressure drop [bar]", title = "",
    overlays=ovs_plots_2
)

  ╠═╡ =#

# ╔═╡ c8c35ed3-bd6c-40ac-8049-b603b003f603
#=╠═╡
plot_contourf(
    df_grid_wall;  xcol=:log_ΔTSF_, ycol=:T_wall, xlogbase=:log, zcol=:COP_2,
    levels=2.0:0.1:4.0,
    overlays=ovs_plots_3
)

  ╠═╡ =#

# ╔═╡ 33c1996a-aade-400a-9b2f-259db12c0835
#=╠═╡
df_grid_wall|> names
  ╠═╡ =#

# ╔═╡ 79b17009-14ee-422e-8eb6-9e81c270b602
#=╠═╡
begin
	# --- 1) Pack the SR grid with the same keys as df_grid_wall -------------------
	df_sr_grid_wall = DataFrame()
	df_sr_grid_wall[!, xcol] = Xgrid_ord[!, xcol]
	df_sr_grid_wall[!, ycol] = Xgrid_ord[!, ycol]
	for t in targets
	    df_sr_grid_wall[!, t] = Ygrid_hat[!, t]         # unitless during math is OK
	end
	
	# --- 2) Join and build error columns (unit-safe) ------------------------------
	true_tbl = select(df_grid_wall, vcat([xcol, ycol], targets))
	pred_tbl = select(df_sr_grid_wall, vcat([xcol, ycol], targets))
	
	rename!(true_tbl, Dict(t => Symbol(t, "_true") for t in targets))
	rename!(pred_tbl, Dict(t => Symbol(t, "_sr")   for t in targets))
	
	cmp = innerjoin(true_tbl, pred_tbl, on=[xcol, ycol])
	
	for t in targets
	    y  = ustrip.(cmp[!, Symbol(t, "_true")])
	    ŷ  = ustrip.(cmp[!, Symbol(t, "_sr")])
	    cmp[!, Symbol(t, "_Δ")]    = ŷ .- y
	    cmp[!, Symbol(t, "_rel%")] = 100 .* (cmp[!, Symbol(t, "_Δ")] ./ max.(abs.(y), eps()))
	end
	
	# --- 3) Summary table of metrics ---------------------------------------------
	function grid_metrics(cmp, targets)
	    rows = NamedTuple[]
	    for t in targets
	        y  = ustrip.(cmp[!, Symbol(t, "_true")])
	        ŷ  = ustrip.(cmp[!, Symbol(t, "_sr")])
	        e  = ŷ .- y
	        push!(rows, (target=t,
	                     RMSE = sqrt(mean(e.^2)),
	                     MAE  = mean(abs.(e)),
	                     MAPE = 100*mean(abs.(e) ./ max.(abs.(y), eps())),
	                     R2   = 1 - sum(e.^2) / max(sum((y .- mean(y)).^2), eps())))
	    end
	    DataFrame(rows)
	end
	
	metrics_tbl = grid_metrics(cmp, targets)  # ← tidy numeric comparison
	
	# --- 4) A few “shape” comparisons --------------------------------------------
	# Peak shift between true and SR fields
	function peak_shift(cmp, xcol, ycol, z::Symbol)
	    yt = ustrip.(cmp[!, Symbol(z, "_true")]); it = argmax(yt)
	    yp = ustrip.(cmp[!, Symbol(z, "_sr")]);   ip = argmax(yp)
	    (; Δx = ustrip(cmp[ip, xcol] - cmp[it, xcol]),
	       Δy = ustrip(cmp[ip, ycol] - cmp[it, ycol]),
	       max_true = yt[it], max_sr = yp[ip])
	end
	peak_Q_cond = peak_shift(cmp, xcol, ycol, :Q_cond)
	
	# Dice coefficient for superlevel sets (iso-contour agreement)
	function dice_at(cmp, z::Symbol, τ::Real)
	    A = ustrip.(cmp[!, Symbol(z, "_true")]) .>= τ
	    B = ustrip.(cmp[!, Symbol(z, "_sr")])   .>= τ
	    2*sum(A .& B) / max(sum(A) + sum(B), 1)
	end
	dice_Q5000 = dice_at(cmp, :Q_cond, 5_000.0)
	
	# Gradient-direction alignment (cosine similarity of ∇Z)
	function grad_alignment(cmp, xcol, ycol, z::Symbol, nx::Int, ny::Int)
	    Zt = reshape(ustrip.(cmp[!, Symbol(z, "_true")]), ny, nx)
	    Zp = reshape(ustrip.(cmp[!, Symbol(z, "_sr")]),   ny, nx)
	    ∂x(A) = @views 0.5 .* (A[:,3:end] .- A[:,1:end-2])
	    ∂y(A) = @views 0.5 .* (A[3:end,:] .- A[1:end-2,:])
	    gx_t, gy_t = ∂x(Zt)[2:end-1,:], ∂y(Zt)[:,2:end-1]
	    gx_p, gy_p = ∂x(Zp)[2:end-1,:], ∂y(Zp)[:,2:end-1]
	    num = gx_t .* gx_p .+ gy_t .* gy_p
	    den = sqrt.(gx_t.^2 .+ gy_t.^2) .* sqrt.(gx_p.^2 .+ gy_p.^2) .+ eps()
	    mean(num ./ den)
	end
	cos_grad_Q = grad_alignment(cmp, xcol, ycol, :Q_cond, length(xvals), length(yvals))
	
end
  ╠═╡ =#

# ╔═╡ 94138b2d-0fd6-494e-9fd1-ef5aecadaabf
md" # Plotting functions"

# ╔═╡ 2097d037-7f79-4fe2-b683-b74b77d98752
function plot_pred_vs_actual_long(df_long; markersize=2, width=1200, row_height=240, band_pct=0.10)
    targets = collect(unique(df_long.target))
    n = length(targets)
    fig = Figure(resolution=(width, max(200, row_height*n)))
    colsize!(fig.layout, 1, Relative(1.0))

    for (i, t) in enumerate(targets)
        df_t = df_long[df_long.target .== t, :]

        # let AoG handle plotting of unitful columns
        ax = Axis(fig[i, 1]; title=string(t), aspect=1)
        draw!(ax, data(df_t) * mapping(:actual, :predicted) *
                  visual(Scatter; markersize=markersize, strokewidth=0))

        # compute numeric limits from data (unit-safe)
        a = ustrip.(df_t.actual)
        p = ustrip.(df_t.predicted)
        lo = min(minimum(a), minimum(p))
        hi = max(maximum(a), maximum(p))
        pad = (hi - lo == 0) ? 1e-6 : 0.02*(hi - lo)
        lo -= pad;  hi += pad

        xlims!(ax, lo, hi)   # Makie’s documented API for setting limits
        ylims!(ax, lo, hi)   # (no xlimits/ylimits field access)  ⟵

        # 1:1 and ±band overlays (numeric)
        lines!(ax, [lo, hi], [lo, hi], linewidth=2)
        lines!(ax, [lo, hi], [(1+band_pct)*lo, (1+band_pct)*hi], linestyle=:dash)
        lines!(ax, [lo, hi], [(1-band_pct)*lo, (1-band_pct)*hi], linestyle=:dash)
    end

    resize_to_layout!(fig)
    fig
end

# ╔═╡ bb0f32e5-f72f-41d1-af47-ed1cbb571bb1
#=╠═╡

# ------------- grid builder (Makie/Contour expects a rectangular grid) -------------
function _grid_from_df(df::DataFrame, xcol::Symbol, ycol::Symbol, zcol::Symbol;
                       invx = identity, agg = first)
    df2 = combine(groupby(df, [ycol, xcol]), zcol => agg => :_z)

    xs = sort(unique(df2[!, xcol]))
    ys = sort(unique(df2[!, ycol]))
    x2i = Dict(x => i for (i,x) in enumerate(xs))
    y2j = Dict(y => j for (j,y) in enumerate(ys))

    Z = fill(NaN, length(xs), length(ys))  # size(Z) == (nx, ny)
    @inbounds for r in eachrow(df2)
        Z[x2i[r[xcol]], y2j[r[ycol]]] = Float64(r._z)
    end

    keep_x = vec(any(isfinite, Z; dims=2))
    keep_y = vec(any(isfinite, Z; dims=1))
    Z = Z[keep_x, keep_y]
    xs = xs[keep_x]; ys = ys[keep_y]

    X = Float64.(invx.(xs))
    Y = Float64.(ys)
    @assert all(diff(X) .> 0) && all(diff(Y) .> 0)
    return X, Y, Z
end
  ╠═╡ =#

# ╔═╡ 3ec7080c-2972-484b-805e-3d9fe0a7d9d1
# choose reasonable default levels for iso-lines
function _default_iso_levels(Z; n::Int=8)
    vals = vec(Z)
    m = isfinite.(vals)
    if !any(m)
        return n               # fallback: let Makie auto-pick n levels
    end
    zmin, zmax = extrema(vals[m])
    if zmin == zmax
        return [zmin]          # constant field: one contour at that value
    else
        return collect(range(zmin, zmax; length=n))
    end
end


# ╔═╡ 7148e421-b7e8-4f3c-a0a9-da3d14dd32e5
# map selected numeric levels to custom strings (with tolerance)
_level_labelformatter(labels::AbstractDict; tol=1e-9) = let tbl = Dict(float(k)=>string(v) for (k,v) in labels)
    l -> begin
        # try exact hit first
        if haskey(tbl, float(l))
            return tbl[float(l)]
        end
        # tolerant match (useful when levels came from computations)
        for (k,v) in tbl
            if isapprox(k, l; atol=tol, rtol=0)
                return v
            end
        end
        # fallback: print the numeric value if not mapped
        @sprintf("%.3g", l)
    end
end

# ╔═╡ c21e1161-99f5-4d15-adb8-7f7c669eadb0

ovs_plots = [
    # Reynolds thresholds — leg 1 (red) with halo
    (z=:Re_bhe_leg1, levels=[2300], color=:red, linewidth=1.3, linestyle=:dash, labels=true,
     #halo=true, halocolor=:white, halowidth=1.2,
     labelformatter=_level_labelformatter(Dict(2300 => "onset of turbulence - down"))),

    (z=:Re_bhe_leg1, levels=[4000], color=:red, linewidth=1.3, linestyle=:dot, labels=true,
     labelformatter=_level_labelformatter(Dict(4000 => "full turbulence - down"))),

    # Reynolds thresholds — leg 2 (magenta instead of blue) with halo
    (z=:Re_bhe_leg2, levels=[2300], color=:red, linewidth=1.3, linestyle=:dash, labels=true,
     labelformatter=_level_labelformatter(Dict(2300 => "onset of turbulence - up"))),

    (z=:Re_bhe_leg2, levels=[4000], color=:red, linewidth=1.3, linestyle=:dot, labels=true,
     labelformatter=_level_labelformatter(Dict(4000 => "full turbulence - up"))),


    # f_pump guides — dense, unobtrusive (no halo)
    (z=:r_N_pump, levels=(0.3:0.1:1)[2:end-1], color=CHARCOAL, linewidth=.5, labels=false),

    # f_pump endpoints — highlighted w/ labels + halo
    (z=:r_N_pump, levels=[1.0], color=CHARCOAL, linewidth=2.0, labels=true,
     labelformatter=_level_labelformatter(Dict(1.0 => "f = 1.0"))),

    (z=:r_N_pump, levels=[0.3], color=CHARCOAL, linewidth=2.0, labels=true,
     labelformatter=_level_labelformatter(Dict(0.3 => "f = 0.3"))),

    (z=:P_pump_P_max, levels=[1.0], color=CHARCOAL, linewidth=2.5, labels=true, linestyle=:dashdot, labelformatter=_level_labelformatter(Dict(1.0 => "P pump max"))),
]


# ╔═╡ 72ec6c2e-498a-47c1-a15f-31d1e5cd14c7
ovs_plots

# ╔═╡ 7ceca66f-f75a-4ba5-aa44-1f2b6939da77

ovs_plots_3 = [
    (z=:Re_bhe_leg1, levels=[2300], color=:red, linewidth=1.3, linestyle=:dash, labels=true,
     labelformatter=_level_labelformatter(Dict(2300 => "onset of turbulence - down"))),

    (z=:Re_bhe_leg1, levels=[4000], color=:red, linewidth=1.3, linestyle=:dot, labels=true,
     labelformatter=_level_labelformatter(Dict(4000 => "full turbulence - down"))),

    (z=:Re_bhe_leg2, levels=[2300], color=:red, linewidth=1.3, linestyle=:dash, labels=true,
     labelformatter=_level_labelformatter(Dict(2300 => "onset of turbulence - up"))),

    (z=:Re_bhe_leg2, levels=[4000], color=:red, linewidth=1.3, linestyle=:dot, labels=true,
     labelformatter=_level_labelformatter(Dict(4000 => "full turbulence - up"))),


    # f_pump guides — dense, unobtrusive (no halo)
    (z=:r_N_pump, levels=(0.3:0.1:1)[2:end-1], color=CHARCOAL, linewidth=.5, labels=false),

    # f_pump endpoints — highlighted w/ labels + halo
    (z=:r_N_pump, levels=[1.0], color=CHARCOAL, linewidth=2.0, labels=true,
     labelformatter=_level_labelformatter(Dict(1.0 => "f = 1.0"))),
]


# ╔═╡ 23f9f68a-9e3c-4bbe-8abe-1b851206162f
#=╠═╡
function _add_isolines!(ax::Makie.Axis, df::DataFrame, xcol::Symbol, ycol::Symbol, zcol::Symbol;
                        xlogbase::Symbol = :log10,
                        levels = nothing,
                        color::Union{Symbol,RGBAf,RGBf} = :black,
                        linestyle::Symbol = :solid,
                        linewidth::Real = 1.5,
                        labels::Bool = false,
                        labelformatter::Function = l -> @sprintf("%.3g", l),
                        labelsize::Int = 20,
                        halo::Bool = false,                # <— NEW
                        halowidth::Real = 1.5,             # <— NEW
                        halocolor::Union{Symbol,RGBAf,RGBf} = :white)  # <— NEW

    invx = xlogbase === :log10 ? exp10 : xlogbase === :log ? exp : identity
    X, Y, Z = _grid_from_df(df, xcol, ycol, zcol; invx=invx)

    levs =
        levels === nothing ? _default_iso_levels(Z) :
        (levels isa AbstractVector && isempty(levels)) ? _default_iso_levels(Z) :
        (levels isa AbstractVector ? Float64.(collect(levels)) : levels)

    # 1) optional halo pass (no labels)
    if halo
        contour!(ax, X, Y, Z;
            levels = levs, color = halocolor, linewidth = linewidth + halowidth,
            linestyle = linestyle, labels = false, overdraw = false)
    end

    # 2) colored line with labels
    contour!(ax, X, Y, Z;
        levels = levs, color = color, linestyle = linestyle, linewidth = linewidth,
        labels = labels, labelformatter = labelformatter, labelsize = labelsize,
        labelcolor = color, overdraw = false)

    return ax
end

  ╠═╡ =#

# ╔═╡ 7520ab74-810f-431a-b24e-724e24595b13


_resolve_label(sym::Symbol, fallback::String, map::Dict{Symbol,String},
               kw::Union{Nothing,AbstractString}) =
    kw === nothing ? get(map, sym, fallback) : String(kw)

# ╔═╡ 1d8f15e6-0462-4f86-9a97-3364e3e73bc0


# --- helpers ---
_normalize_labelmap(x)::Dict{Symbol,String} = begin
    out = Dict{Symbol,String}()
    for (k,v) in pairs(x)          # works for Dict and NamedTuple
        out[Symbol(k)] = string(v)
    end
    out
end

# ╔═╡ 3afecaa8-1c77-4b26-b534-2c419ddfb027
#=╠═╡
# pick a ridge curve by optimizing z over the *other* axis
# wrt = :x  -> for each x, choose y that extremizes z  (y = y*(x))   [old behavior]
# wrt = :y  -> for each y, choose x that extremizes z  (x = x*(y))   [what you want]
function _ridge(df::DataFrame, xcol::Symbol, ycol::Symbol, zcol::Symbol;
                invx = identity, mode::Symbol = :max, wrt::Symbol = :x)
    X, Y, Z = _grid_from_df(df, xcol, ycol, zcol; invx=invx)  # size(Z) == (length(X), length(Y))

    xs = Float64[]; ys = Float64[]
    if wrt === :x
        # iterate over x-slices -> y*(x)
        for i in eachindex(X)
            best = (mode === :max) ? -Inf : Inf
            bestj = 0
            col = @view Z[i, :]
            @inbounds for j in eachindex(Y)
                v = col[j]
                if isfinite(v) && ((mode === :max && v > best) || (mode === :min && v < best))
                    best = v; bestj = j
                end
            end
            bestj != 0 && (push!(xs, X[i]); push!(ys, Y[bestj]))
        end
    elseif wrt === :y
        # iterate over y-slices -> x*(y)
        for j in eachindex(Y)
            best = (mode === :max) ? -Inf : Inf
            besti = 0
            row = @view Z[:, j]
            @inbounds for i in eachindex(X)
                v = row[i]
                if isfinite(v) && ((mode === :max && v > best) || (mode === :min && v < best))
                    best = v; besti = i
                end
            end
            besti != 0 && (push!(xs, X[besti]); push!(ys, Y[j]))
        end
    else
        error("zridge_wrt must be :x or :y")
    end
    return xs, ys
end

  ╠═╡ =#

# ╔═╡ 2164a2d3-59a0-4573-bd12-b62431f46c04
# Dense max/min over x for each fixed y using Interpolations.jl (linear, irregular grid OK)
function _ridge_resample_wrt_y(X::AbstractVector, Y::AbstractVector, Z::AbstractMatrix;
                               mode::Symbol=:max, nsamples::Int=401)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        sel = isfinite.(Z[:, j])
        if count(sel) ≥ 2
            xj = X[sel]; zj = Z[sel, j]
            itp = interpolate((xj,), zj, Gridded(Linear()))  # linear on irregular grid
            xd = range(first(xj), last(xj), length=nsamples)
            vals = itp.(xd)
            k = (mode === :max) ? argmax(vals) : argmin(vals)
            push!(xs, xd[k]); push!(ys, Y[j])
        end
    end
    xs, ys
end

# Same, but per-slice smoothing spline with Dierckx (cubic), then dense evaluate
# Requires: using Dierckx

# ╔═╡ 318a7ca3-1d6a-4a32-9fc7-c2a0072c92d3
function _ridge_spline_wrt_y(X::AbstractVector, Y::AbstractVector, Z::AbstractMatrix;
                             mode::Symbol=:max, nsamples::Int=401, s::Real=0.0)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        sel = isfinite.(Z[:, j])
        if count(sel) ≥ 2
            xj = X[sel]; zj = Z[sel, j]
            spl = Dierckx.Spline1D(xj, zj; k=3, s=s)     # smoothing s≥0
            xd  = range(first(xj), last(xj), length=nsamples)
            vals = spl.(xd)
            k = (mode === :max) ? argmax(vals) : argmin(vals)
            push!(xs, xd[k]); push!(ys, Y[j])
        end
    end
    xs, ys
end

# Optional post-smoothing of the ridge path x*(y)
# Requires: using Dierckx

# ╔═╡ c1333a2c-28a8-4a31-a049-8479e43f6c0e
# Vertex (x at extremum) of quadratic through three points (x1,y1),(x2,y2),(x3,y3)
# Works for non-uniform spacing. Returns NaN if the fit is degenerate.
function _parabolic_vertex(x1,y1,x2,y2,x3,y3)
    denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
    if denom == 0
        return NaN
    end
    a = (x1*(y3 - y2) + x2*(y1 - y3) + x3*(y2 - y1)) / denom
    b = (x1^2*(y2 - y3) + x2^2*(y3 - y1) + x3^2*(y1 - y2)) / denom
    if a == 0
        return NaN
    end
    return -b / (2a)  # x at extremum
end

# Parabolic refinement for wrt=:y  (x* as a smooth function of y)

# ╔═╡ 1177eea2-1d56-4090-8b70-a636f404b4fb
function _ridge_parabolic_wrt_y_mask_interp(X::AbstractVector, Y::AbstractVector,
                                            Z::AbstractMatrix, W::AbstractMatrix;
                                            mode::Symbol=:max, wmin::Real=-Inf, wmax::Real=Inf)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        mj = isfinite.(Z[:,j]) .& isfinite.(W[:,j])
        idxs = findall(mj)
        length(idxs) < 2 && continue

        xi, zi, wi = X[idxs], Z[idxs, j], W[idxs, j]
        w_itp = interpolate((xi,), wi, Gridded(Linear()))

        # pick best grid point
        k = mode === :max ? argmax(zi) : argmin(zi)
        x̂ = xi[k]

        # parabolic refinement if neighbors exist
        if 1 < k < length(xi)
            x1,y1 = xi[k-1], zi[k-1]
            x2,y2 = xi[k],   zi[k]
            x3,y3 = xi[k+1], zi[k+1]
            xv = _parabolic_vertex(x1,y1,x2,y2,x3,y3)
            if isfinite(xv) && (min(x1,x3)-eps() ≤ xv ≤ max(x1,x3)+eps())
                curv = (y1 - 2y2 + y3)
                if (mode === :max && curv < 0) || (mode === :min && curv > 0)
                    x̂ = xv
                end
            end
        end

        # interpolated limiter check
        ŵ = w_itp(x̂)
        if isfinite(ŵ) && (wmin ≤ ŵ ≤ wmax)
            push!(xs, x̂); push!(ys, Y[j])
        end
    end
    xs, ys
end


# ╔═╡ 1d782694-8f7d-4d1f-a1bc-b7646ac01b23
function _ridge_parabolic_wrt_y(X::AbstractVector, Y::AbstractVector, Z::AbstractMatrix;
                                mode::Symbol=:max)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        zcol = view(Z, :, j)
        m = isfinite.(zcol)
        n = count(m)
        n < 2 && continue
        xi = X[m]; zi = zcol[m]

        k = (mode === :max) ? argmax(zi) : argmin(zi)
        x̂ = xi[k]  # fallback = grid point

        if 1 < k < length(xi)  # have neighbors
            x1,y1 = xi[k-1], zi[k-1]
            x2,y2 = xi[k],   zi[k]
            x3,y3 = xi[k+1], zi[k+1]
            xv = _parabolic_vertex(x1,y1,x2,y2,x3,y3)
            if isfinite(xv) && min(x1,x3) - eps() ≤ xv ≤ max(x1,x3) + eps()
                # ensure it's the right extremum shape
                # For maxima we expect concave (a<0); for minima convex (a>0).
                # We approximatively check by central second diff sign.
                curv = (y1 - 2y2 + y3)
                if (mode === :max && curv < 0) || (mode === :min && curv > 0)
                    x̂ = xv
                end
            end
        end
        push!(xs, x̂); push!(ys, Y[j])
    end
    xs, ys
end

# Optional post-smoothing of the ridge path: x(y) -> smooth cubic spline
# (Requires: using Dierckx)

# ╔═╡ a35e92ea-bee1-4117-9ccf-f2016be29dc3
function _smooth_path(xs::AbstractVector, ys::AbstractVector; s::Real=0.0)
    spl = Dierckx.Spline1D(ys, xs; k=3, s=s)
    return spl.(ys), ys
end

# ╔═╡ 17aa8588-7089-4247-aaad-0ca54e1d888a
# Split ridge into continuous segments by jump/Gap/slope rules
function _split_discontinuities(xs::AbstractVector, ys::AbstractVector;
        xjump_frac::Real = 0.08,      # break if |Δx| > xjump_frac * (xmax-xmin)
        ygap_factor::Real = 3.0,      # break if |Δy| > ygap_factor * median(Δy)
        slope_max::Real = Inf)        # break if |Δx/Δy| > slope_max

    n = length(xs)
    n ≥ 2 || return Vector{UnitRange{Int}}()

    # basic scales
    xr = maximum(xs) - minimum(xs)
    dy = diff(ys)
    dx = diff(xs)
    med_dy = median(abs.(dy[isfinite.(dy)]))
    med_dy = (isfinite(med_dy) && med_dy > 0) ? med_dy : 1.0

    segs = Vector{UnitRange{Int}}()
    s = 1
    for k in 1:n-1
        jump_x  = xr > 0 ? abs(dx[k]) > xjump_frac * xr : false
        gap_y   = abs(dy[k]) > ygap_factor * med_dy
        steep   = (slope_max < Inf) ? (abs(dx[k]) / max(abs(dy[k]), eps()) > slope_max) : false
        if jump_x || gap_y || steep || !isfinite(dx[k]) || !isfinite(dy[k])
            if k >= s+1
                push!(segs, s:k)
            end
            s = k + 1
        end
    end
    if n >= s+1
        push!(segs, s:n)
    end
    segs
end

# Optional: discard ridge points too close to left/right edges

# ╔═╡ 286269bd-f4db-4690-81f4-b4c45edeb24d
function _exclude_edge_points!(xs::Vector{Float64}, ys::Vector{Float64}; edge_frac=0.0)
    edge_frac ≤ 0 && return xs, ys
    xmin, xmax = minimum(xs), maximum(xs)
    pad = edge_frac * (xmax - xmin)
    keep = @. (xs > xmin + pad) & (xs < xmax - pad) & isfinite(xs) & isfinite(ys)
    return xs[keep], ys[keep]
end

# ╔═╡ 34522b07-e28b-4fa6-b2bf-2c7c73895aa2
# contiguous True-runs helper (if you don't already have it)
function _true_runs(b::AbstractVector{Bool})
    runs = Vector{UnitRange{Int}}()
    i = 1
    while i <= length(b)
        if b[i]
            s = i
            while i <= length(b) && b[i]
                i += 1
            end
            push!(runs, s:(i-1))
        else
            i += 1
        end
    end
    return runs
end

# ╔═╡ 599e27a4-1aba-4de4-a76d-fb09511de993
function _mask_intervals_spline(xmin::Real, xmax::Real, sw::Dierckx.Spline1D;
                                wmin::Real, wmax::Real, nsamples::Int=401)
    xd = range(xmin, xmax; length=nsamples)
    w  = sw.(xd)
    keep = (w .>= wmin) .& (w .<= wmax)
    runs = _true_runs(keep)
    isempty(runs) && return Tuple{Float64,Float64}[]

    intervals = Tuple{Float64,Float64}[]
    for r in runs
        # crude ends
        a0, b0 = xd[first(r)], xd[last(r)]
        # refine left edge if we have a sample just before the run
        al = a0
        if first(r) > 1
            xl = xd[first(r)-1]
            # which boundary did we cross?
            wl = sw(xl)
            if wl < wmin
                al = try find_zero(x -> sw(x) - wmin, (xl, a0), Bisection()) catch; a0 end
            elseif wl > wmax
                al = try find_zero(x -> sw(x) - wmax, (xl, a0), Bisection()) catch; a0 end
            end
        end
        # refine right edge if we have a sample just after the run
        br = b0
        if last(r) < length(xd)
            xr = xd[last(r)+1]
            wr = sw(xr)
            if wr < wmin
                br = try find_zero(x -> sw(x) - wmin, (b0, xr), Bisection()) catch; b0 end
            elseif wr > wmax
                br = try find_zero(x -> sw(x) - wmax, (b0, xr), Bisection()) catch; b0 end
            end
        end
        (al < br) && push!(intervals, (al, br))
    end
    intervals
end

# Smooth constrained ridge: for each fixed y, maximize/minimize z(x) s.t. w(x)∈[wmin,wmax].
# z and w are smoothed with cubic splines on irregular x; optimization is Brent on each interval.

# ╔═╡ 4b5d8e48-6775-4ec6-ab28-3e96b9954af3
function _ridge_spline_opt_wrt_y_mask(
    X::AbstractVector, Y::AbstractVector, Z::AbstractMatrix, W::AbstractMatrix;
    mode::Symbol=:max, s_z::Real=0.0, s_w::Real=0.0, nsamples_mask::Int=401
)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        sel = isfinite.(Z[:,j]) .& isfinite.(W[:,j])
        count(sel) ≥ 2 || continue
        xj, zj, wj = X[sel], Z[sel,j], W[sel,j]

        # smooth splines (cubic; s>0 to denoise)
        sz = Dierckx.Spline1D(xj, zj; k=3, s=s_z)
        sw = Dierckx.Spline1D(xj, wj; k=3, s=s_w)

        # feasible intervals from mask spline
        intervals = _mask_intervals_spline(first(xj), last(xj), sw, zridge_min, zridge_max;
                                           nsamples=nsamples_mask)
        isempty(intervals) && continue

        best_v = (mode === :max) ? -Inf : Inf
        best_x = NaN
        for (a,b) in intervals
            # Brent 1D (bounded)
            obj = mode === :max ? (x->-sz(x)) : (x->sz(x))
            res = Optim.optimize(obj, a, b, Brent())
            x̂  = Optim.minimizer(res)
            v̂  = sz(x̂)
            if (mode === :max && v̂ > best_v) || (mode === :min && v̂ < best_v)
                best_v = v̂; best_x = x̂
            end
        end
        isfinite(best_x) && (push!(xs, best_x); push!(ys, Y[j]))
    end
    xs, ys
end

# ╔═╡ c9aa37be-53b4-4f23-9898-b5179f7117f2
# smooth, interpolated-mask ridge: maximize/minimize z(x) for each fixed y,
# but only where w(x) lies in [wmin, wmax].
function _ridge_resample_wrt_y_mask_interp(X::AbstractVector, Y::AbstractVector,
                                           Z::AbstractMatrix, W::AbstractMatrix;
                                           mode::Symbol=:max, nsamples::Int=401,
                                           wmin::Real=-Inf, wmax::Real=Inf)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        m = isfinite.(Z[:, j]) .& isfinite.(W[:, j])
        if count(m) ≥ 2
            xj, zj, wj = X[m], Z[m, j], W[m, j]

            # linear on irregular grids
            z_itp = interpolate((xj,), zj, Gridded(Linear()))
            w_itp = interpolate((xj,), wj, Gridded(Linear()))

            xd   = range(first(xj), last(xj), length=nsamples)
            wd   = w_itp.(xd)
            keep = (wd .>= wmin) .& (wd .<= wmax)

            runs = _true_runs(keep)
            isempty(runs) && continue

            best_val = mode === :max ? -Inf : Inf
            best_x   = NaN
            for r in runs
                # need ≥2 points in a sub-interval to be meaningful
                if length(r) < 2; continue; end
                vals = z_itp.(xd[r])
                k    = mode === :max ? argmax(vals) : argmin(vals)
                v    = vals[k]
                if (mode === :max && v > best_val) || (mode === :min && v < best_val)
                    best_val = v
                    best_x   = xd[r[k]]
                end
            end
            isfinite(best_x) && (push!(xs, best_x); push!(ys, Y[j]))
        end
    end
    xs, ys
end


# ╔═╡ 9d0ffbe1-deff-4de9-a0d7-9e7c8df5bfaf
function _ridge_parabolic_resample_wrt_y_mask_interp(
    X::AbstractVector, Y::AbstractVector, Z::AbstractMatrix, W::AbstractMatrix;
    mode::Symbol=:max, nsamples::Int=601, wmin::Real=-Inf, wmax::Real=Inf
)
    xs = Float64[]; ys = Float64[]
    for j in eachindex(Y)
        # finite pairs only
        sel = isfinite.(Z[:, j]) .& isfinite.(W[:, j])
        count(sel) ≥ 2 || continue

        xj = X[sel]; zj = Z[sel, j]; wj = W[sel, j]

        z_itp = interpolate((xj,), zj, Gridded(Linear()))
        w_itp = interpolate((xj,), wj, Gridded(Linear()))

        xd   = range(first(xj), last(xj), length=nsamples)
        keep = (w_itp.(xd) .>= wmin) .& (w_itp.(xd) .<= wmax)
        runs = _true_runs(keep)
        isempty(runs) && continue

        best_val = mode === :max ? -Inf : Inf
        best_x   = NaN

        for r in runs
            # sample z on this allowed sub-interval
            vals = z_itp.(xd[r])
            k    = mode === :max ? argmax(vals) : argmin(vals)
            xi   = xd[r[k]]; zi = vals[k]

            # parabolic refine if we have neighbors inside the run
            if 1 < k < length(r)
                x1, x2, x3 = xd[r[k-1]], xd[r[k]], xd[r[k+1]]
                y1, y2, y3 = z_itp(x1), z_itp(x2), z_itp(x3)
                xv = _parabolic_vertex(x1,y1,x2,y2,x3,y3)
                # keep refinement only if inside the allowed sub-interval AND still allowed by w
                if isfinite(xv) && x1 - eps() ≤ xv ≤ x3 + eps()
                    wv = w_itp(xv)
                    if isfinite(wv) && (wmin ≤ wv ≤ wmax)
                        xi = xv; zi = z_itp(xv)
                    end
                end
            end

            if (mode === :max && zi > best_val) || (mode === :min && zi < best_val)
                best_val = zi
                best_x   = xi
            end
        end

        isfinite(best_x) && (push!(xs, best_x); push!(ys, Y[j]))
    end
    return xs, ys
end

# ╔═╡ bf82a17d-f524-4aa7-ab17-eadcc045d6cb
#=╠═╡
function plot_contourf(df::DataFrame;
    xcol::Symbol=:log_ΔTSF_, ycol::Symbol=:f_, zcol::Symbol=:COP_3,
    xlogbase::Symbol=:log10,
    levels=10,
    palette=:viridis, lowclip=nothing, highclip=nothing,
    labels_for_lines::Bool=false,
    overlays::Vector = Vector{Any}(),
    return_axis::Bool=false,

    # labeling options
    xlabel::Union{Nothing,AbstractString}=nothing,
    ylabel::Union{Nothing,AbstractString}=nothing,
    zlabel::Union{Nothing,AbstractString}=nothing,
    title::Union{Nothing,AbstractString}=nothing,
    labelmap::Union{AbstractDict,NamedTuple}=Dict(),

    # ridge (extremum curve of z)
    zridge_enable::Bool = false,
    zridge_wrt::Symbol  = :x,        # :x | :y
    zridge_mode::Symbol = :max,      # :max | :min
    zridge_color::Union{Symbol,RGBAf,RGBf} = :red,
    zridge_linestyle::Symbol = :solid,
    zridge_linewidth::Real = 2.0,
    zridge_halo::Bool = false,
    zridge_halowidth::Real = 1.5,
    zridge_halocolor::Union{Symbol,RGBAf,RGBf} = :white,
    zridge_labels::Bool = false,
    zridge_label::Union{Nothing,AbstractString} = nothing,
    zridge_labelsize::Int = 14,
    zridge_method::Symbol = :resample,
    zridge_nsamples::Int = 401,
    zridge_spline_s::Real = 0.0,
    zridge_postsmooth::Bool = false,
    zridge_post_s::Real = 0.0,
    zridge_xjump_frac::Real = 0.08,
    zridge_ygap_factor::Real = 3.0,
    zridge_slope_max::Real = Inf,
    zridge_exclude_edge_frac::Real = 0.0,
	
    # ridge domain limit via *interpolated* limiter column (e.g. :r_N_pump)
    zridge_limit_col::Union{Nothing,Symbol} = nothing,
    zridge_min::Real = -Inf,
    zridge_max::Real =  Inf,
    zridge_limit_interp::Bool = true,   # use the interpolated-mask variant


    # base rendering mode
    base_style::Symbol = :filled,    # :filled | :lines | :none
    show_colorbar::Bool = true,
    line_width::Real = 1.2,          # for :lines
    line_color::Union{Nothing,Symbol,RGBAf,RGBf} = nothing  # single color if set
)
    invx = xlogbase === :log10 ? exp10 : xlogbase === :log ? exp : identity
    xlab_default = (invx === identity) ? string(xcol) : "ΔTSF"

    lm = _normalize_labelmap(labelmap)
    xname = _resolve_label(xcol, xlab_default, lm, xlabel)
    yname = _resolve_label(ycol, string(ycol), lm, ylabel)
    zname = _resolve_label(zcol, string(zcol), lm, zlabel)
    plt_title = title === nothing ? zname : String(title)

    fig = Figure()
    ax  = Axis(fig[1,1]; xlabel=xname, ylabel=yname, title=plt_title)
    if hasproperty(ax, :clip)
        ax.clip = true
    elseif hasproperty(ax, :scene) && hasproperty(ax.scene, :clip)
        ax.scene.clip = true
    end

    # grid + helpers
    Xg, Yg, Zg = _grid_from_df(df, xcol, ycol, zcol; invx=invx)
    zmin, zmax = extrema(vec(Zg))
    levs = levels isa Integer ? levels : collect(levels)
    bfmt(x) = @sprintf("%.4g", x)

    if base_style === :filled
        # filled bands via AoG
        pal = (lowclip !== nothing || highclip !== nothing) ?
              clipped(from_continuous(palette); low=lowclip, high=highclip) :
              from_continuous(palette)
        base  = data(df) * mapping(xcol => invx => xname, ycol, zcol)
        fc    = levels isa Integer ? filled_contours(bands=levels) : filled_contours(levels=collect(levels))
        specs = base * fc
        sc    = scales(Color=(; palette=pal))
        draw!(ax, specs, sc)

        # outline boundaries
        boundaries = levels isa Integer ?
                     collect(range(zmin, zmax; length=levels+1)) :
                     collect(sort(levels))
        contour!(ax, Xg, Yg, Zg; levels=boundaries, color=:white, linewidth=.3,
                 labels=labels_for_lines, labelformatter=bfmt, labelsize=20)

        if show_colorbar
            boundaries = sort!(boundaries)
            nbands = max(length(boundaries)-1, 1)
            cm = cgrad(palette, nbands; categorical=true)
            Colorbar(fig[1,2];
                colormap=cm,
                limits=(first(boundaries), last(boundaries)),
                ticks=(boundaries, map(bfmt, boundaries)),
                label=zname, vertical=true, lowclip=lowclip, highclip=highclip)
        end

    elseif base_style === :lines
        # contour lines only
        if line_color === nothing
            cm = cgrad(palette)
            contour!(ax, Xg, Yg, Zg;
                levels=levs, colormap=cm, colorrange=(zmin, zmax),
                linewidth=line_width,
                labels=labels_for_lines, labelformatter=bfmt, labelsize=18)
            if show_colorbar
                chosen = levels isa Integer ? collect(range(zmin, zmax; length=levels)) : Float64.(levs)
                Colorbar(fig[1,2];
                    colormap=cm, limits=(zmin, zmax),
                    ticks=(chosen, map(bfmt, chosen)),
                    label=zname, vertical=true)
            end
        else
            contour!(ax, Xg, Yg, Zg;
                levels=levs, color=line_color, linewidth=line_width,
                labels=labels_for_lines, labelformatter=bfmt, labelsize=18)
            # (no colorbar by default)
        end

    elseif base_style === :none
        # draw nothing; overlays/ridge only
        # (no colorbar)
    else
        error("base_style must be :filled, :lines, or :none")
    end

    # overlays
    for ov in overlays
        zsym = get(ov, :z, nothing); zsym === nothing && error("overlay missing :z")
        _add_isolines!(ax, df, xcol, ycol, zsym;
            xlogbase        = xlogbase,
            levels          = get(ov, :levels, []),
            color           = get(ov, :color, :black),
            linestyle       = get(ov, :linestyle, :solid),
            linewidth       = get(ov, :linewidth, 1.5),
            labels          = get(ov, :labels, false),
            labelformatter  = get(ov, :labelformatter, l -> @sprintf("%.3g", l)),
            labelsize       = get(ov, :labelsize, 10),
            halo            = get(ov, :halo, false),
            halowidth       = get(ov, :halowidth, 1.5),
            halocolor       = get(ov, :halocolor, :white))
    end

	if zridge_enable
	    Xg, Yg, Zg = _grid_from_df(df, xcol, ycol, zcol; invx=invx)
	
	    Wg = nothing
	    if zridge_limit_col !== nothing
	        _, _, Wg = _grid_from_df(df, xcol, ycol, zridge_limit_col; invx=invx)
	    end
	
	    if zridge_wrt === :y
	        if zridge_method === :resample
	            if zridge_limit_interp && Wg !== nothing
	                xs, ys = _ridge_resample_wrt_y_mask_interp(
	                    Xg, Yg, Zg, Wg; mode=zridge_mode,
	                    nsamples=zridge_nsamples,
	                    wmin=zridge_min, wmax=zridge_max)
	            else
	                # your previous (discrete-mask) path:
	                allow = Wg === nothing ? nothing :
	                        @. isfinite(Wg) & (Wg >= zridge_min) & (Wg <= zridge_max)
	                xs, ys = _ridge_resample_wrt_y(
	                    Xg, Yg, Zg; mode=zridge_mode,
	                    nsamples=zridge_nsamples, allow=allow)
	            end
				elseif zridge_method === :parabolic
			    if zridge_limit_interp && Wg !== nothing
			        xs, ys = _ridge_parabolic_resample_wrt_y_mask_interp(
			            Xg, Yg, Zg, Wg; mode=zridge_mode,
			            nsamples=zridge_nsamples,             # e.g. 601
			            wmin=zridge_min, wmax=zridge_max
			        )
			    else
			        # (fallbacks unchanged)
			        allow = Wg === nothing ? nothing :
			                @. isfinite(Wg) & (Wg >= zridge_min) & (Wg <= zridge_max)
			        xs, ys = _ridge_parabolic_wrt_y(Xg, Yg, Zg; mode=zridge_mode)  # or a discrete-mask variant if you keep one
			    end

	        else
	            error("zridge_method must be :resample or :parabolic for wrt=:y")
	        end
	    else
	        # (wrt=:x variant unchanged; add a similar mask if you ever need it)
	        xs, ys = _ridge(df, xcol, ycol, zcol; invx=invx, mode=zridge_mode, wrt=:x)
	    end
	    # 2) optional path smoothing of x*(y)
	    if zridge_postsmooth && !isempty(xs)
	        xs, ys = _smooth_path(xs, ys; s=zridge_post_s)  # requires Dierckx
	    end
	
	    # 3) drop edge-hugging points (optional)
	    xs, ys = _exclude_edge_points!(collect(xs), collect(ys);
	                                   edge_frac = zridge_exclude_edge_frac)
	
	    # 4) split into continuous segments (jump/gap/slope rules)
	    segs = _split_discontinuities(xs, ys;
	        xjump_frac = zridge_xjump_frac,
	        ygap_factor = zridge_ygap_factor,
	        slope_max = zridge_slope_max)
	
	    # 5) draw each continuous segment (with optional halo)
	    for r in segs
	        if zridge_halo
	            lines!(ax, xs[r], ys[r]; color=zridge_halocolor,
	                   linewidth = zridge_linewidth + zridge_halowidth,
	                   linestyle = zridge_linestyle)
	        end
	        lines!(ax, xs[r], ys[r]; color=zridge_color,
	               linewidth = zridge_linewidth, linestyle = zridge_linestyle)
	    end
	
	    # 6) optional label at the end of last segment
	    if zridge_labels && !isempty(segs)
	        r = last(segs)
	        lbl = zridge_label === nothing ?
	            string(zridge_mode == :max ? "max " : "min ", zname, " wrt ",
	                   (zridge_wrt == :x ? xname : yname)) :
	            String(zridge_label)
	        text!(ax, lbl; position=(xs[last(r)], ys[last(r)]),
	              align=(:left, :center), offset=(8, 0),
	              color=zridge_color, fontsize=zridge_labelsize)
	    end
	end

    resize_to_layout!(fig)
    return return_axis ? (fig, ax) : fig
end

  ╠═╡ =#

# ╔═╡ b570d453-43f2-4a12-a9d6-889752dd5aeb
#=╠═╡
function plot_heatmap(df::DataFrame;
                      xcol::Symbol = :log_ΔTSF_,
                      ycol::Symbol = :f_,
                      zcol::Symbol = :COP_3,
                      xlogbase::Symbol = :log10)  # use :log if natural log

    invx = xlogbase === :log10 ? exp10 :
           xlogbase === :log   ? exp   :
           error("xlogbase must be :log10 or :log")

    plt = data(df) *
          mapping(                     # <-- POSitional: x, y, z
              xcol => invx => "ΔTSF",  # log(ΔTSF) -> ΔTSF for plotting
              ycol,
              zcol
          ) *
          visual(Heatmap)

    fig, ax = draw(
        plt;
        axis = (
            xlabel = "ΔTSF",           # linear axis (default)
            ylabel = string(ycol),
            title  = string(zcol, " over ΔTSF vs ", ycol),
        )
        # colormap = :viridis,
        # colorbar = true
    )

    return fig
end

  ╠═╡ =#

# ╔═╡ 21cd24f0-4d8b-4298-89e1-7ba7817c2eef
#=╠═╡

"""
    process_grid_results(
      values1, values2,
      results_thermal, results_hydraulic,
      sys_thermal, sys_hydraulic;
      name1=:param1, name2=:param2
    ) -> DataFrame

Loop over the 2D sweep defined by `values1`×`values2`, extract all the same
metrics you did in 1D, and return a DataFrame with columns:

  name1, name2, ṁ_sf, ṁ_ref, T_in, …, COP_3

where `name1`/`name2` label the two sweep axes.
"""
function process_grid_results(
    values1::AbstractVector,
    values2::AbstractVector,
    results_thermal::AbstractMatrix,
    results_hydraulic::AbstractMatrix,
    sys_thermal,
    sys_hydraulic;
    name1::Symbol = :param1,
    name2::Symbol = :param2,
	target_Q_cond::Real = NaN,
)
    n1, n2 = length(values1), length(values2)
    @assert size(results_thermal) == (n1, n2)
    @assert size(results_hydraulic) == (n1, n2)

    rows = NamedTuple[]

    for i in 1:n1, j in 1:n2
        therm = results_thermal[i,j]
        hyd   = results_hydraulic[i,j]

        # — your thermal extractions —
        ṁ_sf      = exp( therm[sys_thermal.γṁ_sf] )
        ṁ_ref     = exp( therm[sys_thermal.γṁ_ref] )
        T_in       = therm[sys_thermal.cycle.T_i_b] .- 273.15
        T_mid      = therm[sys_thermal.cycle.evaporator.T_mid_1] .- 273.15
        T_out      = therm[sys_thermal.cycle.T_o_b] .- 273.15
        T_bot      = therm[sys_thermal.bhe_thermal.thermal.T_b] .- 273.15
        T_comp_in  = therm[sys_thermal.cycle.evaporator.T_out_2] .- 273.15
        T_wall     = 10.0
        ΔTSF         = T_in - T_out

        fr_SH      = zap(therm[sys_thermal.cycle.evaporator.uz_fr_sh])
        fr_EV      = 1 - fr_SH

        T_evap     = therm[sys_thermal.cycle.evaporator.sh.propref.Tsat] .- 273.15
        Q_evap     = therm[sys_thermal.Q_evap]
        Q_cond     = -therm[sys_thermal.Q_cond]
        W_comp     = therm[sys_thermal.W_comp]
        COP_1      = therm[sys_thermal.COP1]

        U_EV       = therm[sys_thermal.cycle.evaporator.ev.base.U]
        U_SH       = therm[sys_thermal.cycle.evaporator.sh.base.U]

        # — your hydraulic extractions —
        Δp_bhe     = hyd[sys_hydraulic.bhe.Δp]          / 1e5
        Δp_phx     = hyd[sys_hydraulic.phx.base.Δp]     / 1e5
        Δp_tot     = hyd[sys_hydraulic.pump.base.Δp]    / 1e5
        P_pump     = hyd[sys_hydraulic.pump.base.P_elec]
        N_pump     = hyd[sys_hydraulic.pump.base.N]
        r_N_pump     = hyd[sys_hydraulic.pump.base.N]/hyd[sys_hydraulic.pump.base.N_max]
		η_pump = hyd[sys_hydraulic.pump.base.η]

		P_pump_P_max = P_pump / 185#hyd[sys_hydraulic.pump.base.P_elec_max]
		
		pump_compressor = P_pump/W_comp

        COP_2      = Q_cond / (W_comp + P_pump)

        # — PHX Reynolds number —
        ρ_phx, μ_phx, ṁ_phx = hyd[sys_hydraulic.phx.props.ρ],
                              hyd[sys_hydraulic.phx.props.μ],
                              hyd[sys_hydraulic.phx.base.ṁ]
        b_phx, W_phx, np_phx, ϕ_phx = hyd[sys_hydraulic.phx.b],
                                     hyd[sys_hydraulic.phx.W],
                                     hyd[sys_hydraulic.phx.np],
                                     hyd[sys_hydraulic.phx.ϕ]
        nch       = (np_phx - 1)/2
        A_flow    = nch * W_phx * (2*b_phx)
        v_phx     = abs(ṁ_phx)/(ρ_phx*A_flow)
        dh_phx    = 2*b_phx/ϕ_phx
        Re_phx    = ρ_phx * v_phx * dh_phx / μ_phx

        # — BHE legs (replace L1 & L2 with your real field‐names) —
        ρ1 = hyd[sys_hydraulic.bhe.leg1.props.ρ]
        μ1 = hyd[sys_hydraulic.bhe.leg1.props.μ]
        D1 = hyd[sys_hydraulic.bhe.leg1.Di]
        Re_bhe_leg1 = abs(hyd[sys_hydraulic.bhe.ṁ]) * D1 / (μ1 * π*(D1/2)^2)

        ρ2 = hyd[sys_hydraulic.bhe.leg2.props.ρ]
        μ2 = hyd[sys_hydraulic.bhe.leg2.props.μ]
        D2 = hyd[sys_hydraulic.bhe.leg2.Di]
        Re_bhe_leg2 = abs(hyd[sys_hydraulic.bhe.ṁ]) * D2 / (μ2 * π*(D2/2)^2)

        # — BHE thermal & resistances —
        η_bh   = therm[sys_thermal.bhe_thermal.ϵ]
        R_1b    = therm[sys_thermal.bhe_thermal.utube_res.R_t1_b]
        R_2b    = therm[sys_thermal.bhe_thermal.utube_res.R_t2_b]
        R_12    = therm[sys_thermal.bhe_thermal.utube_res.R_t1_t2]
        R_f1b   = therm[sys_thermal.bhe_thermal.thermal.R_f1_b]
        R_f2b   = therm[sys_thermal.bhe_thermal.thermal.R_f2_b]
        R_f1f2  = therm[sys_thermal.bhe_thermal.thermal.R_f1_f2]

		# Try these corrected paths:
		ρ_l = therm[sys_thermal.cycle.evaporator.ev.propref.ρb]  # liquid density
		ρ_g = therm[sys_thermal.cycle.evaporator.ev.propref.ρd]  # gas density  
		σ = therm[sys_thermal.cycle.evaporator.ev.propref.σ]     # surface tension
		
		# Hydraulic diameter from evaporator geometry
		D_h = 2 * therm[sys_thermal.cycle.evaporator.ev.b] / therm[sys_thermal.cycle.evaporator.ev.ϕ]
		
		# Calculate Bond number
		Bd = (ρ_l - ρ_g) * 9.81 * D_h^2 / σ

		# scalar, not a 1-element vector
		f_measured = sys_thermal.cycle.compressor.f
		
		# build the payload without params
		payload = (;  ṁ_sf, ṁ_ref, T_in, T_mid, T_out, ΔTSF, T_bot, T_evap, T_comp_in,
		             Q_evap, Q_cond, W_comp, COP_1, U_EV, U_SH, fr_EV, fr_SH,
		             Δp_bhe, Δp_phx, Δp_tot, P_pump, N_pump, r_N_pump, η_pump, COP_2,
		             Re_phx, Re_bhe_leg1, Re_bhe_leg2,
		             η_bh, R_1b, R_2b, R_12, R_f1b, R_f2b, R_f1f2,
		             pump_compressor, P_pump_P_max, Bd)
		
		# include measured frequency only if it won’t collide with the axis name
		if :f_ ∉ (name1, name2)
		    payload = (; payload..., f_ = f_measured)
		else
		    # (optional) keep it under another name if you still want it
		    payload = (; payload..., f_meas = f_measured)
		end
		
		# now make the row: params first, then payload (no duplicates)
		row = (; (name1 => values1[i]), (name2 => values2[j]), payload...)
		push!(rows, row)

    end
    df = DataFrame(rows)

    # (optional but handy)
    sort!(df, [name2, :ΔTSF])   # keep rows grouped by f then ΔT
	
	if isnan(target_Q_cond)
    # Compute per-f max Q_cond
	    gdf = groupby(df, name2)
	    transform!(gdf, :Q_cond => maximum => :Qmax_same_f)
	
	    # Aux heater needed to meet the best point at this f
	    df.P_aux_heating = max.(0, df.Qmax_same_f .- df.Q_cond)

	else
			
	    # Aux heater needed to meet the best point at this f
	    df.P_aux_heating = max.(0, target_Q_cond .- df.Q_cond)

	end
    # COP_3 normalized within each f-group
    df.COP_3 = (df.Q_cond .+ df.P_aux_heating) ./ (df.W_comp .+ df.P_pump .+ df.P_aux_heating)

	g = groupby(df, name2)
	
	# per-f maximum and per-f normalized COP_3
	transform!(g,
	    :COP_3 => (v -> fill(maximum(skipmissing(v)), length(v))) => :COP_3_max_f,
	    :COP_3 => (v -> v ./ maximum(skipmissing(v)))             => :r_COP_3,
	)
	
	# optional: deficit from the local max
	df.p_r_COP_3 = df.r_COP_3 .- 1.0

	freq_col = :f_ ∈ names(df) ? :f_ : (:f_meas ∈ names(df) ? :f_meas : name2)
	df.r_fc_fp = df[!, freq_col] ./ df.N_pump

    return df
end


  ╠═╡ =#

# ╔═╡ fff91fa6-8ffb-458a-946d-244b5b13c3cd
function multithreaded_grid_sweep(thermal_prob, hydraulic_prob, sys_thermal, 
                                  base_params, param_sym1, values1, param_sym2, values2, solvers;
                                  previous_sol_thermal=nothing, previous_sol_hydraulic=nothing)
    n1 = length(values1)
    n2 = length(values2)
    
    # Preallocate solution arrays
    sols_thermal = Matrix{Any}(undef, n1, n2)
    sols_hydraulic = Matrix{Any}(undef, n1, n2)
    
    # Calculate Manhattan distances
    center_i = ceil(Int, n1/2)
    center_j = ceil(Int, n2/2)
    distances = [abs(i - center_i) + abs(j - center_j) for i in 1:n1, j in 1:n2]
    max_distance = maximum(distances)
    
    # Compute parameter ranges for normalization
    range1 = maximum(values1) - minimum(values1)
    range2 = maximum(values2) - minimum(values2)
    range1 = range1 == 0 ? 1.0 : range1  # Avoid division by zero
    range2 = range2 == 0 ? 1.0 : range2
    
    # Solve center point with staged parameters
    center_params = [
        base_params...,
        param_sym1 => values1[center_i],
        param_sym2 => values2[center_j]
    ]
    
    # Stage 1: Solve with base_params only if no previous solution
    if previous_sol_thermal === nothing
        staged_sol_thermal = solve_thermal_problem(
            thermal_prob, 
            nothing, 
            base_params,  # Only base_params
            solvers
        )
        staged_sol_hydraulic = solve_hydraulic_problem(
            hydraulic_prob, 
            nothing, 
            base_params,
            staged_sol_thermal, 
            sys_thermal, 
            solvers
        )
    else
        staged_sol_thermal = previous_sol_thermal
        staged_sol_hydraulic = previous_sol_hydraulic
    end

    # Stage 2: Solve with full parameters using staged solution as initial condition
    sols_thermal[center_i, center_j] = solve_thermal_problem(
        thermal_prob, 
        staged_sol_thermal,  # Use staged solution
        center_params,       # Full parameters with sweeps
        solvers
    )
    
    sols_hydraulic[center_i, center_j] = solve_hydraulic_problem(
        hydraulic_prob, 
        staged_sol_hydraulic, 
        center_params,
        sols_thermal[center_i, center_j], 
        sys_thermal, 
        solvers
    )

    # Create wavefronts (points grouped by distance from center)
    wavefronts = [Tuple{Int,Int}[] for _ in 0:max_distance]
    for i in 1:n1, j in 1:n2
        d = distances[i, j]
        if d > 0  # Skip center point
            push!(wavefronts[d], (i, j))
        end
    end
    
    # Progress bar
    total_points = n1 * n2 - 1  # Minus center point
    p = Progress(total_points, 1, "Grid sweep progress: ")
    
    # Process wavefronts in order
    for d in 1:max_distance
        points = wavefronts[d]
        isempty(points) && continue
        
        # Process all points in current wavefront
        for idx in 1:length(points)
            i, j = points[idx]
            
            # Find nearest solved neighbor in PARAMETER SPACE
            best_thermal = nothing
            best_hydraulic = nothing
            min_param_dist = Inf
            
            # Search all solved points for closest in parameter space
            for i2 in 1:n1, j2 in 1:n2
                # Only consider solved points
                if isassigned(sols_thermal, i2, j2) && sols_thermal[i2, j2] !== nothing
                    # Compute normalized parameter distance
                    param_dist = sqrt(
                        ((values1[i] - values1[i2]) / range1)^2 +
                        ((values2[j] - values2[j2]) / range2)^2
                    )
                    
                    # Update best solution if closer
                    if param_dist < min_param_dist
                        min_param_dist = param_dist
                        best_thermal = sols_thermal[i2, j2]
                        best_hydraulic = sols_hydraulic[i2, j2]
                    end
                end
            end
            
            # Fallback to center solution if no neighbors found
            if best_thermal === nothing
                best_thermal = sols_thermal[center_i, center_j]
                best_hydraulic = sols_hydraulic[center_i, center_j]
            end
            
            # Solve systems
            current_params = [
                base_params...,
                param_sym1 => values1[i],
                param_sym2 => values2[j]
            ]
            
            # Step solving for better convergence
            try
                sol_thermal = solve_thermal_problem(
                    thermal_prob, best_thermal, current_params, solvers
                )
                sol_hydraulic = solve_hydraulic_problem(
                    hydraulic_prob, best_hydraulic, current_params,
                    sol_thermal, sys_thermal, solvers
                )
            # Store solutions
            sols_thermal[i, j] = sol_thermal
            sols_hydraulic[i, j] = sol_hydraulic
            catch e
                # Fallback to center solution if direct solve fails
                @warn "Solve failed at ($i,$j), retrying with center solution"
                sol_thermal = solve_thermal_problem(
                    thermal_prob, sols_thermal[center_i, center_j], current_params, solvers
                )
                sol_hydraulic = solve_hydraulic_problem(
                    hydraulic_prob, sols_hydraulic[center_i, center_j], current_params,
                    sol_thermal, sys_thermal, solvers
                )
            # Store solutions
            sols_thermal[i, j] = sol_thermal
            sols_hydraulic[i, j] = sol_hydraulic
            end
            
            
            next!(p)
        end
    end
    
    return sols_thermal, sols_hydraulic
end

# ╔═╡ 1f4b84eb-7ca1-478e-8e4f-376ac069ea1e
# ╠═╡ show_logs = false
sols_thermal2, sols_hydraulic2 = multithreaded_grid_sweep(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params, :log_ΔTSF_, log_ΔTSF_values, :f_, f_values, solvers;
    previous_sol_thermal=nothing,
    previous_sol_hydraulic=nothing
)

# ╔═╡ 6a624b9b-1e84-4e98-b1a8-60fd8e95cb13
#=╠═╡
df_grid = process_grid_results(
    log_ΔTSF_values, f_values,
    sols_thermal2, sols_hydraulic2,
    sys_thermal2, sys_hydraulic2;
    name1 = :log_ΔTSF_,
    name2 = :f_
);
  ╠═╡ =#

# ╔═╡ 02497be0-7d2f-43d1-b365-d4a8d72f045e
#=╠═╡
plot_heatmap(df_grid; xlogbase = :log)
  ╠═╡ =#

# ╔═╡ fc5900c1-9704-4384-9a29-be545faf74d6
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :COP_3, levels=20, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "Base scenario")
  ╠═╡ =#

# ╔═╡ 8262db7e-b994-4afa-b823-661a0bb4949b
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :Re_phx, levels=20,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Reynolds PHX [-]", title = "Base scenario")
  ╠═╡ =#

# ╔═╡ eb8cd9d6-9237-4723-84d3-3f293d74d1e9
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :COP_1, levels=3.0:0.05:4.1,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 1 (Compressor only) [-]", title = "Base scenario")
  ╠═╡ =#

# ╔═╡ 0b7c9c92-e263-4a91-bbf9-d32cb0ac7024
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :COP_2, levels=3.0:0.05:4.1,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 2 (Including pump) [-]", title = "Base scenario")
  ╠═╡ =#

# ╔═╡ 6eab6576-10af-432d-abf6-42f424541e5b
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :T_evap,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation temperature [°C]", title = "Base scenario")
  ╠═╡ =#

# ╔═╡ 8dea2e5f-82b5-408e-843d-c25b7d265d48
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :ṁ_sf,)
  ╠═╡ =#

# ╔═╡ b153b909-a8b2-4e5d-abe7-a3c1bad81358
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :fr_EV,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation fraction [-]", title = "Base scenario")
  ╠═╡ =#

# ╔═╡ 2be85177-c03c-48de-8d4f-c2053fcf2a13
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :η_bh, levels = 20,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Borehole effectiveness [-]", title = "")
  ╠═╡ =#

# ╔═╡ a24d0373-b198-42f6-a6bb-c2ede7259244
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :ṁ_ref,levels = 20)

  ╠═╡ =#

# ╔═╡ 74455461-5277-4fee-a53b-edf34e9e26f4
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :ṁ_sf,levels = 20)

  ╠═╡ =#

# ╔═╡ dc49f979-5f43-4234-864d-acfef268f5f9
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :P_pump, levels = 20,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Pump power [W]", title = "")
  ╠═╡ =#

# ╔═╡ 50af6ba3-363d-4ca6-977a-7eb35919b27b
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :Q_evap,levels = 20)

  ╠═╡ =#

# ╔═╡ 88cbfed8-9de9-4421-8dd0-8f42b82234bd
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :Q_cond,levels = 20,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Condenser heat duty [W]", title = "")

  ╠═╡ =#

# ╔═╡ 9035962a-af30-4870-8ccb-c340892e81c3
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :r_fc_fp, levels = 20,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Pump frequency ratio [-]", title = "")
  ╠═╡ =#

# ╔═╡ b6a370b1-bfe7-4324-b5a6-6018ddf55e12
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :T_out,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Outlet temperature [°C]", title = "Base scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ d7b463f3-4a9d-4572-9ae6-eb5d4ea47531
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :T_in,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Inlet temperature [°C]", title = "Base scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 60d5e58e-237b-4f13-ba00-c50a5315ffd7
#=╠═╡
plot_contourf(df_grid; xlogbase = :log, zcol = :T_bot,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Bottom temperature [°C]", title = "Base scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ b9879207-d8b7-4465-b899-e7656e48e5fc
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:COP_2, levels=2.7:0.1:4.2,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 2 (Including pump) [-]", title = "",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 09b7509b-a308-49a4-966d-33d5ef7f9cb1
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:COP_3, levels=10,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 1945625f-482a-4eda-985d-7112b4cf432d
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:U_EV,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporator UA [W/K]", title = "",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 7386981b-7839-4588-aa68-5c4b3d461b06
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:U_SH,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Superheater UA [W/K]", title = "",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ daed3493-da2b-489c-b7ae-5bcd66bfc03b
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:η_bh, levels = 0.0:0.1:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Borehole effectiveness [-]", title = "Base scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 5e034021-1116-4d1d-8034-d582c8a9643f
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:T_in, levels = -0:0.5:10, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Outlet BH temperature [°C]", title = "Base Scenario",    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 1ce3b2c5-71ba-4a1f-ab48-961585e57943
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:T_evap, levels = -5:0.5:5, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation temperature [°C]", title = "Base Scenario",    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 16928040-8730-4f0c-a67c-50f9d1ae852b
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:COP_2, levels=2.7:0.1:4.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 05d571f0-c990-430e-971c-965ba527e814
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:COP_3, levels=2.7:0.1:4.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 9fa61b9f-1d33-444f-9c1e-d7ed70933000
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:p_r_COP_3, palette=:viridis, levels=-.2:.02:.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relative COP 3 deficit [-]", title = "Base scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 2f8b1672-934b-4b93-9e05-6169ba730cc8
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:η_pump, levels = 0.0:0.05:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Pump efficiency [-]", title = "Base scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 6ebdb8db-d5cc-4d02-a471-1e29a29ab42b
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:r_COP_3, levels = [0.0, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999,1.0],
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relative COP 3 [-]", title = ""
)
  ╠═╡ =#

# ╔═╡ 97b66208-6b79-40dd-bfd9-b710807a4e3f
#=╠═╡
plot_contourf(df_grid;
    xlogbase=:log, zcol=:COP_3, levels=2.7:0.1:4.1,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 4.0,
    zridge_labels = false)

  ╠═╡ =#

# ╔═╡ 960c408c-3040-4a00-b285-c01080cb1c02
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:Q_cond,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Condenser heat duty [W]", title = "",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 644901ea-7778-4e7d-a5cf-44e34df77827
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:R_f1f2,
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 57e985e5-023c-4389-8e08-2185bed0751a
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:fr_EV, levels=0.0:0.1:1.0, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation fraction [-]", title = "", 
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ bb90f897-17d3-4a9d-a462-fd6d3f4d98e7
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:pump_compressor, levels = 0.0:0.01:0.2,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Pump/Compressor power ratio [-]", title = "",
    overlays=ovs_plots_2
)
  ╠═╡ =#

# ╔═╡ 6acbadd7-94f6-4399-ac5f-3a2724ce641c
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:r_fc_fp, levels = 0.2:0.2:1.8,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relataive frequency", title = "",     base_style = :lines, show_colorbar = false, labels_for_lines=true,
)
  ╠═╡ =#

# ╔═╡ 29e4a45f-7a0b-425b-90d5-b2fc84f43d10
#=╠═╡
plot_contourf(
    df_grid; xlogbase=:log, zcol=:Δp_tot, levels = 0.2:0.2:1.8,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relataive frequency", title = "",     base_style = :lines, show_colorbar = true, labels_for_lines=true,
)
  ╠═╡ =#

# ╔═╡ 195ab453-b1f2-4632-9420-091f7f702ed7
# ╠═╡ show_logs = false
sols_thermal_cold, sols_hydraulic_cold = multithreaded_grid_sweep(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params_cold, :log_ΔTSF_, log_ΔTSF_values, :f_, f_values, solvers;
    previous_sol_thermal=nothing,
    previous_sol_hydraulic=nothing
)

# ╔═╡ 517540bb-548e-4b97-8e96-bd674f8a9d18
#=╠═╡
df_grid_cold = process_grid_results(
    log_ΔTSF_values, f_values,
    sols_thermal_cold, sols_hydraulic_cold,
    sys_thermal2, sys_hydraulic2;
    name1 = :log_ΔTSF_,
    name2 = :f_
);
  ╠═╡ =#

# ╔═╡ e086552e-9ce0-4353-a16b-395d2986d1a9
#=╠═╡
plot_contourf(
    df_grid_cold; xlogbase=:log, zcol=:η_bh, levels = 0.0:0.1:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Borehole effectiveness [-]", title = "Cold ground scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 4217415e-772c-44d2-8c3d-2e85cf01fb1e
#=╠═╡
plot_contourf(
    df_grid_cold; xlogbase=:log, zcol=:T_in, levels = -0:0.5:10, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Outlet BH temperature [°C]", title = "Cold ground scenario",    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 629aaec1-eca1-4b91-910d-5579b49c8df2
#=╠═╡
plot_contourf(
    df_grid_cold; xlogbase=:log, zcol=:T_evap, levels = -10:0.5:0, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation temperature [°C]", title = "Cold ground scenario",
    overlays=ovs_plots_2
)
  ╠═╡ =#

# ╔═╡ 6ae2052e-907c-49b9-bb51-e2e3b0c641bc
#=╠═╡
plot_contourf(
    df_grid_cold; xlogbase=:log, zcol=:COP_3, levels=2.5:0.1:4.0, palette=:viridis,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "Cold ground scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ a0332504-e854-4b16-bb9d-e5f4ca277f40
#=╠═╡
plot_contourf(
    df_grid_cold; xlogbase=:log, zcol=:p_r_COP_3, palette=:viridis, levels=-.2:.02:.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relative COP 3 deficit [-]", title = "Colder wall scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 0bc8ed9a-6d4d-4983-a204-9b54ef7a1fa9
#=╠═╡
plot_contourf(df_grid_cold;
    xlogbase=:log, zcol=:COP_3, levels=2.7:0.1:4.1,
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 4.0,
    zridge_labels = false)

  ╠═╡ =#

# ╔═╡ bba28e41-70b6-409a-a4ec-1a1246261e82
#=╠═╡
plot_contourf(df_grid_cold;
    xlogbase=:log, zcol=:COP_3, levels=2.7:0.1:4.1,
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :white,
    zridge_linestyle = :dot,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 4.0,
    zridge_labels = false)

  ╠═╡ =#

# ╔═╡ cf7f3e02-3ac6-495d-a6a8-81711ae0fe57
# ╠═╡ show_logs = false
sols_thermal_water, sols_hydraulic_water = multithreaded_grid_sweep(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params_water, :log_ΔTSF_, log_ΔTSF_values, :f_, f_values, solvers;
    previous_sol_thermal=nothing,
    previous_sol_hydraulic=nothing
)

# ╔═╡ 6011e1fa-1058-4bce-bb65-fcdf25951bfa
#=╠═╡
df_grid_water = process_grid_results(
    log_ΔTSF_values, f_values,
    sols_thermal_water, sols_hydraulic_water,
    sys_thermal2, sys_hydraulic2;
    name1 = :log_ΔTSF_,
    name2 = :f_
);
  ╠═╡ =#

# ╔═╡ b72f0888-4716-48ec-8a99-fc068a988abf
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:η_bh, levels = 0.0:0.1:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Borehole effectiveness [-]", title = "Pure water scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ e6a2deea-3658-42b3-8ef9-83e31be5a300
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:T_in, levels = -0:0.5:10, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Outlet BH temperature [°C]", title = "Pure Water Scenario",    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 5d5f3c09-164b-4a55-9bf1-b2bb90221465
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:T_evap, levels = -5:0.5:5, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation temperature [°C]", title = "Pure water scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ e174ba39-dd8c-4af9-88c1-caa504af2c7b
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:COP_2, levels=2.7:0.1:4.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "Pure water scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 3b86a38d-9b62-4db0-b943-a1b973eb9b3c
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:p_r_COP_3, palette=:viridis, levels=-.2:.02:.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relative COP 3 deficit [-]", title = "Pure Water scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
) |> x -> save("cop3_pure_water_scenario.svg", x)
  ╠═╡ =#

# ╔═╡ 9beedddb-8b5a-4811-98cb-9a26ecd77b52
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:η_pump, levels = 0.0:0.05:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Pump efficiency [-]", title = "Pure water scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ eac9124f-ffa0-4a9a-948d-3c6e5ec0c7a2
#=╠═╡
plot_contourf(
    df_grid_water; xlogbase=:log, zcol=:fr_EV, levels=0.0:0.1:1.0, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation fraction [-]", title = "Pure water scenario", 
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ f28598ab-9f11-4bd5-b5e4-6e4ae2aea350
# ╠═╡ show_logs = false
sols_thermal_pipe, sols_hydraulic_pipe = multithreaded_grid_sweep(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params_pipe, :log_ΔTSF_, log_ΔTSF_values, :f_, f_values, solvers;
    previous_sol_thermal=nothing,
    previous_sol_hydraulic=nothing
)

# ╔═╡ 941dcb04-88be-43ae-bc21-ae3ef9dd0db0
#=╠═╡
df_grid_pipe = process_grid_results(
    log_ΔTSF_values, f_values,
    sols_thermal_pipe, sols_hydraulic_pipe,
    sys_thermal2, sys_hydraulic2;
    name1 = :log_ΔTSF_,
    name2 = :f_
);
  ╠═╡ =#

# ╔═╡ dfc85639-a09b-44fa-99af-61706a192549
#=╠═╡
plot_contourf(
    df_grid_pipe; xlogbase=:log, zcol=:η_bh, levels = 0.0:0.1:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Borehole effectiveness [-]", title = "Wider pipe scenario",
    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 99e84662-a903-424d-b270-d665a52f479f
#=╠═╡
plot_contourf(
    df_grid_pipe; xlogbase=:log, zcol=:T_in, levels = -0:0.5:10, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Outlet BH temperature [°C]", title = "Wider Pipe Scenario",    overlays=ovs_plots
)
  ╠═╡ =#

# ╔═╡ 4175c8ba-346c-4afa-9267-ed7e689cff35
#=╠═╡
plot_contourf(
    df_grid_pipe; xlogbase=:log, zcol=:T_evap, levels = -5:0.5:5, 
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Evaporation temperature [°C]", title = "Wider pipe scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 387e3115-866a-4831-bb59-a18edf749375
#=╠═╡
plot_contourf(
    df_grid_pipe; xlogbase=:log, zcol=:COP_3, levels=3.0:0.1:4.0, palette=:viridis,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "Wider pipe scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ a7e7ede3-1caa-4483-91f3-c26393f2b45a
#=╠═╡
plot_contourf(
    df_grid_pipe; xlogbase=:log, zcol=:p_r_COP_3, palette=:viridis, levels=-.2:.02:.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Relative COP 3 deficit [-]", title = "Wider pipe scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
	
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
) |> x -> save("cop3_wider_pipe_scenario.svg", x)
  ╠═╡ =#

# ╔═╡ 5eed01c3-561f-490b-80d5-23db0846bfe3
#=╠═╡
plot_contourf(
    df_grid_pipe; xlogbase=:log, zcol=:η_pump, levels = 0.0:0.05:1.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "Pump efficiency [-]", title = "Wider pipe scenario"
)
  ╠═╡ =#

# ╔═╡ 8bc61318-43a1-447d-81a5-d2fed06799b2
#=╠═╡
plot_contourf(df_grid_pipe;
    xlogbase=:log, zcol=:COP_3, levels=2.7:0.1:4.1,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 4.0,
    zridge_labels = false,)
  ╠═╡ =#

# ╔═╡ 592a2cfb-9d8a-481e-a39d-3147dc76b643
# ╠═╡ show_logs = false
sols_thermal_long, sols_hydraulic_long = multithreaded_grid_sweep(
    thermal_prob, hydraulic_prob, sys_thermal,
    base_params_long, :log_ΔTSF_, log_ΔTSF_values, :f_, f_values, solvers;
    previous_sol_thermal=nothing,
    previous_sol_hydraulic=nothing
)

# ╔═╡ c2a86042-1d4b-498c-8370-0d38f73a2bfd
#=╠═╡
df_grid_long = process_grid_results(
    log_ΔTSF_values, f_values,
    sols_thermal_long, sols_hydraulic_long,
    sys_thermal2, sys_hydraulic2;
    name1 = :log_ΔTSF_,
    name2 = :f_
);
  ╠═╡ =#

# ╔═╡ d9a2618d-e923-4370-a20d-b8faf44bf06c
#=╠═╡
plot_contourf(
    df_grid_long; xlogbase=:log, zcol=:COP_2, levels=2.7:0.1:4.0,
    xlabel = "ΔT Secondary Fluid [K]", ylabel = "Compressor frequency [Hz]", 
    zlabel = "COP 3 [-]", title = "Pure water scenario",
    overlays=ovs_plots,
    zridge_enable = true,
    zridge_wrt    = :y,     # <-- maximize over x for each y (x = x*(y))
    zridge_mode   = :max,
    zridge_color  = :red,
    zridge_linestyle = :solid,
    zridge_method   = :parabolic,  # <- smoother than :resample
			  
    zridge_linewidth = 3.0,
    zridge_labels = false,
    # smooth limiter on pump speed
    zridge_limit_col   = :r_N_pump,
    zridge_min = 0.3, zridge_max = 1.0,
    zridge_limit_interp = true,
)
  ╠═╡ =#

# ╔═╡ 698de8c1-b7a4-477c-ba42-88e78cabddf2


# ╔═╡ dfdcded3-026e-4925-ac02-ba801488b267


# ╔═╡ 1abb42cf-af03-49a9-837b-902999a89744


# ╔═╡ 51d3286a-d89d-41f5-88d5-bb361e76ff09


# ╔═╡ 5e9e8e04-4aca-416d-a4b4-410b243fd722


# ╔═╡ 47ab6638-1135-452c-ae79-bb5682094a51
md" # Sampling and symbolic regression"

# ╔═╡ 23899c17-5414-405e-9bb8-7f505b45e2f6
begin
	const QMC = QuasiMonteCarlo

	# --- Canonical sampling over bounds (NamedTuple of (lo, hi)) ---
	# Returns Vector{NamedTuple} of length n with keys matching `bounds`.
	sample_canonical(bounds::NamedTuple; n::Int=512) = begin
        syms = collect(keys(bounds))
        lb   = Float64[first(bounds[k]) for k in syms]
        ub   = Float64[last(bounds[k])  for k in syms]

        X = QMC.sample(n, lb, ub, QMC.SobolSample())  # shape: d × n (cols are samples)

        # Build one NamedTuple per *column*/sample
        [NamedTuple{Tuple(syms)}(Tuple(col)) for col in eachcol(X)]
    end

	# --- Merge "vector of pairs" with overrides (last-write wins) ---
	merge_pairs(base::AbstractVector{<:Pair}, overrides::AbstractVector{<:Pair}) = begin
	    d = Dict(base); foreach(p -> (d[first(p)] = last(p)), overrides); collect(pairs(d))
	end

	 """
        make_mapper(; transforms=Dict(), renames=Dict(), drop=Symbol[], base_params=nothing)

    Build a mapper that:
      • passes through EVERY key in θ as `:key => value` by default
      • applies `transforms[key](value, θ, base_params)` when provided
      • otherwise renames via `renames[key] => value`
      • skips keys in `drop`

    `transforms[key]` may return a `Pair`, a `Vector{<:Pair}`, a `NamedTuple`,
    or a plain value (then the target name is taken from `renames[key]` or defaults to the same key).
    """
    function make_mapper(; transforms=Dict{Symbol,Function}(),
                           renames=Dict{Symbol,Symbol}(),
                           drop=Symbol[],
                           base_params=nothing)

        # Return a closure that accepts θ (and optionally base_params if you prefer calling it directly)
        return function (θ::NamedTuple)
            ov = Pair[]
            for s in propertynames(θ)  # iterate canonical keys  # :contentReference[oaicite:2]{index=2}
                s in drop && continue
                val = getproperty(θ, s)

                if haskey(transforms, s)
                    out = transforms[s](val, θ, base_params)
                    if out isa Pair
                        push!(ov, out)
                    elseif out isa AbstractVector{<:Pair}
                        append!(ov, out)
                    elseif out isa NamedTuple
                        append!(ov, collect(pairs(out)))
                    else
                        # plain value: map to renamed key if provided, else same key
                        target = get(renames, s, s)
                        push!(ov, target => out)
                    end
                else
                    target = get(renames, s, s)
                    push!(ov, target => val)  # default pass-through
                end
            end
            return ov
        end
    end

	make_params(base_params::AbstractVector{<:Pair}, θ::NamedTuple; mapper=default_mapper) = merge_pairs(base_params, mapper(θ))
end


# ╔═╡ 312118cb-b321-4d8a-b1d9-62411134a980
begin
	# Pick super-simple physical bounds to start (edit as you wish)
	canon_bounds = (
	    T_wall = ( -5.0+273.15, 15.0+273.15),    # K
	    L_BHE = ( 50, 300),    # m
		ΔTSF = (1, 8),
		f_ = (10, 50),
		
	)

	points = sample_canonical(canon_bounds; n=3600);

end;


# ╔═╡ e319da9f-acb5-43e4-830f-01822b9a19ab
"""
	run_cloud_sweep(thermal_prob, hydraulic_prob, sys_thermal,
					 base_params, points::Vector{<:NamedTuple}, solvers;
					 mapper = default_mapper)

Solve center first, then iterate points in increasing distance from the center,
seeding each solve with the *closest solved* neighbor so far (normalized Euclidean).

Returns (solT::Vector{Any}, solH::Vector{Any}, points::Vector{NamedTuple})
"""
function run_cloud_sweep(thermal_prob, hydraulic_prob, sys_thermal,
                         base_params::AbstractVector{<:Pair},
                         points::AbstractVector{<:NamedTuple},
                         solvers; mapper=default_mapper)

	@assert !isempty(points)
	syms = collect(keys(points[1]))

	# normalize each dim to [0,1] over the sampled set
	mat  = reduce(hcat, [[p[s] for s in syms] for p in points])'  # n×d
	mins = mapslices(minimum, mat; dims=1)[:] 
	rngs = mapslices(x->max(maximum(x)-minimum(x), 1e-12), mat; dims=1)[:]
	embed(v::AbstractVector) = (v .- mins) ./ rngs
	embed(θ::NamedTuple)     = embed([θ[s] for s in syms])

	# pick center as the sample closest to the coordinate-wise median
	med  = mapslices(median, mat; dims=1)[:]
	center_idx = argmin(sum(abs2, embed(mat[i,:]) .- embed(med)) for i in 1:size(mat,1))
	center     = points[center_idx]

	# Stage 1: staged/base → Stage 2: full center (like your 2D code)
	staged_T = solve_thermal_problem(thermal_prob, nothing, base_params, solvers)
	staged_H = solve_hydraulic_problem(hydraulic_prob, nothing, base_params,
									   staged_T, sys_thermal, solvers)
	center_params = make_params(base_params, center; mapper=mapper)
	solT_c = solve_thermal_problem(thermal_prob, staged_T, center_params, solvers)
	solH_c = solve_hydraulic_problem(hydraulic_prob, staged_H, center_params,
									 solT_c, sys_thermal, solvers)

	n = length(points)
	solT = Vector{Any}(undef, n); solH = Vector{Any}(undef, n)
	solved = Int[]; push!(solved, center_idx)
	solT[center_idx] = solT_c; solH[center_idx] = solH_c

	# order remaining by distance from center
	xc = embed(center)
	order = sortperm(1:n; by = i -> sum(abs2, embed(points[i]) .- xc))
	order = filter(!=(center_idx), order)

	for idx in order
		θ  = points[idx]; x = embed(θ)
		# find nearest solved neighbor (brute force is simplest & fine to start)
		nearest = argmin( sum(abs2, x .- embed(points[j])) for j in solved )
		seed_idx = solved[nearest]
		seedT, seedH = solT[seed_idx], solH[seed_idx]
		params = make_params(base_params, θ; mapper=mapper)

		# try with nearest seed; on failure, fallback to center seeds
		local st, sh
		try
			st = solve_thermal_problem(thermal_prob, seedT, params, solvers)
			sh = solve_hydraulic_problem(hydraulic_prob, seedH, params, st, sys_thermal, solvers)
		catch
			st = solve_thermal_problem(thermal_prob, solT_c, params, solvers)
			sh = solve_hydraulic_problem(hydraulic_prob, solH_c, params, st, sys_thermal, solvers)
		end
		solT[idx] = st; solH[idx] = sh
		push!(solved, idx)
	end

	return solT, solH, points
end

# ╔═╡ 31aa8e9b-9122-497f-866f-ef359f7604ad
#=╠═╡
function collect_cloud_results(points, results_thermal, results_hydraulic,
							  sys_thermal, sys_hydraulic; target_Q_cond::Real = NaN)
	rows = NamedTuple[]

	for i in eachindex(points)
		θ     = points[i]
		therm = results_thermal[i]
		hyd   = results_hydraulic[i]

		# --- thermal (unchanged from your code) ---
		ṁ_sf      = exp( therm[sys_thermal.γṁ_sf] )
		ṁ_ref     = exp( therm[sys_thermal.γṁ_ref] )
		T_in       = therm[sys_thermal.cycle.T_i_b] .- 273.15
		T_mid      = therm[sys_thermal.cycle.evaporator.T_mid_1] .- 273.15
		T_out      = therm[sys_thermal.cycle.T_o_b] .- 273.15
		T_bot      = therm[sys_thermal.bhe_thermal.thermal.T_b] .- 273.15
		T_comp_in  = therm[sys_thermal.cycle.evaporator.T_out_2] .- 273.15
		ΔT         = T_in - T_out
		

		fr_SH      = zap(therm[sys_thermal.cycle.evaporator.uz_fr_sh])
		fr_EV      = 1 - fr_SH
		T_evap     = therm[sys_thermal.cycle.evaporator.sh.propref.Tsat] .- 273.15
		Q_evap     = therm[sys_thermal.Q_evap]
		Q_cond     = -therm[sys_thermal.Q_cond]
		W_comp     = therm[sys_thermal.W_comp]
		COP_1      = therm[sys_thermal.COP1]
		U_EV       = therm[sys_thermal.cycle.evaporator.ev.base.U]
		U_SH       = therm[sys_thermal.cycle.evaporator.sh.base.U]

		# --- hydraulic (unchanged) ---
		Δp_bhe     = hyd[sys_hydraulic.bhe.Δp]          / 1e5
		Δp_phx     = hyd[sys_hydraulic.phx.base.Δp]     / 1e5
		Δp_tot     = hyd[sys_hydraulic.pump.base.Δp]    / 1e5
		P_pump     = hyd[sys_hydraulic.pump.base.P_elec]
		N_pump     = hyd[sys_hydraulic.pump.base.N]
		COP_2      = Q_cond / (W_comp + P_pump)

		# --- PHX Reynolds ---
		ρ_phx, μ_phx, ṁ_phx = hyd[sys_hydraulic.phx.props.ρ],
							  hyd[sys_hydraulic.phx.props.μ],
							  hyd[sys_hydraulic.phx.base.ṁ]
		b_phx, W_phx, np_phx, ϕ_phx = hyd[sys_hydraulic.phx.b],
									 hyd[sys_hydraulic.phx.W],
									 hyd[sys_hydraulic.phx.np],
									 hyd[sys_hydraulic.phx.ϕ]
		nch     = (np_phx - 1)/2
		A_flow  = nch * W_phx * (2*b_phx)
		v_phx   = abs(ṁ_phx)/(ρ_phx*A_flow)
		dh_phx  = 2*b_phx/ϕ_phx
		Re_phx  = ρ_phx * v_phx * dh_phx / μ_phx

		# --- BHE legs ---
		ρ1 = hyd[sys_hydraulic.bhe.leg1.props.ρ]
		μ1 = hyd[sys_hydraulic.bhe.leg1.props.μ]
		D1 = hyd[sys_hydraulic.bhe.leg1.Di]
		Re_bhe_leg1 = abs(hyd[sys_hydraulic.bhe.ṁ]) * D1 / (μ1 * π*(D1/2)^2)

		ρ2 = hyd[sys_hydraulic.bhe.leg2.props.ρ]
		μ2 = hyd[sys_hydraulic.bhe.leg2.props.μ]
		D2 = hyd[sys_hydraulic.bhe.leg2.Di]
		Re_bhe_leg2 = abs(hyd[sys_hydraulic.bhe.ṁ]) * D2 / (μ2 * π*(D2/2)^2)

		# --- BHE thermal & resistances ---
		η_bh   = therm[sys_thermal.bhe_thermal.ϵ]
		R_1b    = therm[sys_thermal.bhe_thermal.utube_res.R_t1_b]
		R_2b    = therm[sys_thermal.bhe_thermal.utube_res.R_t2_b]
		R_12    = therm[sys_thermal.bhe_thermal.utube_res.R_t1_t2]
		R_f1b   = therm[sys_thermal.bhe_thermal.thermal.R_f1_b]
		R_f2b   = therm[sys_thermal.bhe_thermal.thermal.R_f2_b]
		R_f1f2  = therm[sys_thermal.bhe_thermal.thermal.R_f1_f2]

		# row: include canonical params first so they become columns too
		row = (; θ...,
				ṁ_sf, ṁ_ref, T_in, T_mid, T_out, ΔT, T_bot, T_evap, T_comp_in,
				Q_evap, Q_cond, W_comp, COP_1, U_EV, U_SH, fr_EV, fr_SH,
				Δp_bhe, Δp_phx, Δp_tot, P_pump, N_pump, COP_2,
				Re_phx, Re_bhe_leg1, Re_bhe_leg2,
				η_bh, R_1b, R_2b, R_12, R_f1b, R_f2b, R_f1f2)
		push!(rows, row)
	end

	df = DataFrame(rows)

	# Aux heater and COP_3 (global target or none)
	if isnan(target_Q_cond)
		# group by last canonical param just to have a grouping (tweak if you like)
		lastkey = propertynames(points[1])[end]
		gdf = groupby(df, lastkey)
		transform!(gdf, :Q_cond => maximum => :Qmax_same_group)
		df.P_aux_heating = max.(0, df.Qmax_same_group .- df.Q_cond)
	else
		df.P_aux_heating = max.(0, target_Q_cond .- df.Q_cond)
	end
	df.COP_3 = (df.Q_cond .+ df.P_aux_heating) ./ (df.W_comp .+ df.P_pump .+ df.P_aux_heating)
	df
end

  ╠═╡ =#

# ╔═╡ 6ccf5d2b-1e32-4aa6-bea2-f1830e795280
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
begin
	# exceptions:
	transforms = Dict(
	    :ΔTSF      => (v, θ, bp) -> :log_ΔTSF_ => log(v),    # map ΔTSF → log_ΔTSF_
	)
	renames = Dict(
	    :log_ΔTSF_ => :log_ΔTSF_ # (no-op)
	)

	# build the mapper
	mapper = make_mapper(; transforms, renames)


	# Run the sweep (sequential & sturdy)
	solT, solH, pts = run_cloud_sweep(
		thermal_prob, hydraulic_prob, sys_thermal,
		base_params, points, solvers;
		mapper = mapper
	)

	# Pack into one DataFrame
	df_cloud = collect_cloud_results(pts, solT, solH, sys_thermal, sys_hydraulic);
end;
  ╠═╡ =#

# ╔═╡ 0e0ec7a4-31a2-4583-a8e5-0a44f48d2c59
#=╠═╡
df_cloud
  ╠═╡ =#

# ╔═╡ 794b5b4d-62d7-4dda-a379-a7f1d9af7573
md"## 1) Operators (physics-y, NaN-safe)"


# ╔═╡ 36f7dc30-ec0b-493b-bdd0-092c393342fa
begin
	square(x) = x*x
	gmean(x::T,y::T) where {T} = (x ≥ 0 && y ≥ 0) ? sqrt(x*y) : T(NaN)
	hmean(x::T,y::T) where {T} = (x + y) == 0 ? T(NaN) : (2x*y)/(x + y)
	lmean(x::T,y::T) where {T} = (x > 0 && y > 0) ? (x == y ? x : (x-y)/(log(x)-log(y))) : T(NaN)
	amean(x,y) = 0.5(x + y)
		
	const BIN_OPS = [+, -, *, /, lmean, amean]
	const UN_OPS  = [exp, log, tanh, sqrt, square]
end

# ╔═╡ e78e5fac-9201-4208-a643-3583eecae761
md"## 3) Choose columns + (optional) filter rows"

# ╔═╡ a8a0c668-be27-4e64-8a03-039c3bf07df8
#=╠═╡
begin
	ys_m = [:W_comp; :P_pump; :Q_cond; :ṁ_sf; :ṁ_ref]
	xs_m = [:T_wall; :L_BHE; :ΔTSF; :f_]
	
	# optional row filter (example from your earlier step)
	df_sub = copy(df_cloud)
end
  ╠═╡ =#

# ╔═╡ 92b07487-ebec-4026-989a-bf1bcdc9d1c9
#=╠═╡
begin
	# add units to X and y columns
	dfu = copy(df_sub)
	dfu.T_wall .= dfu.T_wall .* us"K"
	dfu.L_BHE  .= dfu.L_BHE  .* us"m"
	dfu.ΔTSF   .= dfu.ΔTSF   .* us"K"
	dfu.f_     .= dfu.f_     .* us"Hz"
	
	dfu.W_comp .= dfu.W_comp .* us"W"
	dfu.Q_cond .= dfu.Q_cond .* us"W"
	dfu.ṁ_sf  .= dfu.ṁ_sf  .* us"kg/s"
	dfu.ṁ_ref .= dfu.ṁ_ref .* us"kg/s"
end

  ╠═╡ =#

# ╔═╡ 3812a920-c1f8-4ecd-ae82-8d5ab03b3522
md"## 4) Fit SR, show best formulas"

# ╔═╡ 76482fe3-b1c0-4a40-a41a-a95362d181ee
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
begin
	model_mt = MultitargetSRRegressor(
	    niterations = 500,
	    binary_operators = [+, -, *, /, amean, lmean],
	    unary_operators  = [sqrt, erf, log, exp],
	    dimensional_constraint_penalty = 1e5,   # strongly discourage unit violations
	    seed = 0, turbo = true
	)

	
	sr_m = sr_fit_multi(dfu; ys=ys_m, x=xs_m, model=model_mt)
	r = sr_m.report
	print_best_equations(r, ys_m)
end
  ╠═╡ =#

# ╔═╡ dd40bbe0-6b5c-4bcc-a371-f6a956b7759b
md"## 5) Predictions on training set + long table"

# ╔═╡ 1a0f0c34-0fe2-44ef-b67b-8ecb51732f7f
# ╠═╡ show_logs = false
#=╠═╡
begin
	# make predictions with the best member per target
	Xtrain   = select(dfu, xs_m)
	Ytrain   = select(dfu, ys_m)
	best     = r.best_idx
	Yhat_tbl = sr_predict_multi(sr_m.mach, Xtrain; idx=best)
	
	# build tidy long table: :row, :target, :actual, :predicted
	dfA = select(Ytrain, ys_m); dfP = select(Yhat_tbl, ys_m)
	dfA.row = 1:nrow(dfA); dfP.row = 1:nrow(dfP)
	longA = stack(dfA, Not(:row); variable_name=:target, value_name=:actual)
	longP = stack(dfP, Not(:row); variable_name=:target, value_name=:predicted)
	df_long = innerjoin(longA, longP, on=[:row, :target])
end
  ╠═╡ =#

# ╔═╡ 65f2bffb-80eb-42f6-9685-bd719b0f61ed
#=╠═╡
eqs_per_target = r.equation_strings    # Vector{Vector{String}}
  ╠═╡ =#

# ╔═╡ e7f2b4bd-7b17-41d8-9c60-a41d88181456
#=╠═╡
best_idxs      = r.best_idx            # Vector{Int}
  ╠═╡ =#

# ╔═╡ 3a4edf4f-d73b-4dd0-a706-da7f14c4801b
#=╠═╡
targets        = ys_m
  ╠═╡ =#

# ╔═╡ 01ee0a4e-6d90-44ea-8119-bfd48892bff1
#=╠═╡
for (i, t) in pairs(targets)           # i is the target index (1..length(ys_m))
    idx = best_idxs[i]                 # pick the best for this target
    println("Best for $t: (idx=$idx) ", eqs_per_target[i][idx])
end
  ╠═╡ =#

# ╔═╡ 19f8c59d-f243-45db-837c-f664193b24a4
#=╠═╡
plot_pred_vs_actual_long(df_long)
  ╠═╡ =#

# ╔═╡ e2a509cf-f4eb-4db2-9cf9-b390e3c36cd8
md"## 2-D parametric study (grid → predict → contour)"


# ╔═╡ ce47fdb8-b5ef-4a17-b1a7-a47240576545
#=╠═╡

# 1) Make a grid (works with unitful df). NOTE: no Float64 coercion!
function make_param_grid(df::DataFrame, xs::Vector{Symbol};
        xcol::Symbol, ycol::Symbol, nx::Int=80, ny::Int=80,
        xscale::Symbol=:identity, xr=nothing, yr=nothing, fixed=:median)

    _ext(v) = extrema(skipmissing(v))

    # x values (don’t log unitful quantities; keep :identity unless dimensionless)
    xvals = xr === nothing ? begin
        xmin, xmax = _ext(df[!, xcol])
        xscale === :identity || @assert !(df[!, xcol] isa AbstractVector{<:DynamicQuantities.Quantity}) "Use xscale=:identity for unitful columns."
        collect(range(xmin, xmax, length=nx))
    end : collect(xr)

    yvals = yr === nothing ? begin
        ymin, ymax = _ext(df[!, ycol])
        collect(range(ymin, ymax, length=ny))
    end : collect(yr)

    # fixed values for other inputs (keep their units)
    fixed_map = Dict{Symbol, Any}()
    for c in xs
        (c == xcol || c == ycol) && continue
        v = fixed === :median  ? Statistics.median(skipmissing(df[!, c])) :
            fixed === :mean    ? Statistics.mean(skipmissing(df[!, c]))   :
            (fixed isa Dict && haskey(fixed, c)) ? fixed[c] :
            Statistics.median(skipmissing(df[!, c]))
        fixed_map[c] = v
    end

    # cartesian grid (unitful); order → reshape to (ny, nx)
    Xgrid = DataFrame()
    Xgrid[!, xcol] = repeat(xvals, inner=ny)    # x slow
    Xgrid[!, ycol] = repeat(yvals, outer=nx)    # y fast
    for c in xs
        (c == xcol || c == ycol) && continue
        Xgrid[!, c] = fill(fixed_map[c], nx*ny)
    end
    return Xgrid, xvals, yvals
end

  ╠═╡ =#

# ╔═╡ c6469779-0f69-4845-9183-719ecd3b7c26
md"## 2) SR helpers (fit/report/predict)"

# ╔═╡ 60a39d0c-f909-48b7-922f-eb88bc11b088
#=╠═╡
begin
	
# Build a simple SR model with sane defaults
sr_model(; niterations::Int=200,
         binary_ops = [+, -, *, /],
         unary_ops  = [cos, exp, log],
         seed::Union{Int,Nothing}=0,         # use `nothing` for non-deterministic runs
         deterministic::Bool=false) =

    SRRegressor(
        niterations        = niterations,
        binary_operators   = binary_ops,
        unary_operators    = unary_ops,
        seed               = seed,           # ← instead of rng=
        deterministic      = deterministic,  # (serial-only if true)
		turbo = true,
    )

	"""
	    sr_fit(df; y, x, model = sr_model(), drop_missings=true)

	Train SymbolicRegression (MLJ) on a DataFrame.
	Returns a NamedTuple with:
	  - mach::Machine
	  - report::Any        (contains the Pareto equations)
	  - Xnames::Vector{Symbol}
	  - yname::Symbol
	"""
	function sr_fit(df::DataFrame; y::Symbol, x::Vector{Symbol},
	                model = sr_model(), drop_missings::Bool=true)

	    data = select(df, vcat(x, [y]))
	    data = drop_missings ? dropmissing(data) : data

	    X = select(data, x)    # MLJ accepts any Tables.jl table
	    yv = data[!, y]

	    mach = machine(model, X, yv)
	    fit!(mach)  # trains and stores the Pareto front

	    r = report(mach)  # equations & metadata
	    return (; mach, report=r, Xnames=x, yname=y)
	end
end

  ╠═╡ =#

# ╔═╡ ae544fa4-9b4a-48a4-ade9-01881147fa8f
begin
	"""
	    sr_evaluate(mach; folds=5, measure=rms)

	Cross-validated performance (MLJ). Returns the PerformanceEvaluation object.
	"""
	sr_evaluate(mach; folds::Int=5, measure=rms) =
	    evaluate!(mach, resampling=CV(nfolds=folds), measure=measure, verbosity=0)

	"""
	    sr_predict_df(mach, Xnew; idx=nothing)

	Predict on a new table `Xnew`. If `idx` is given, uses that Pareto member.
	"""
	function sr_predict_df(mach, Xnew; idx::Union{Nothing,Int}=nothing)
	    idx === nothing ? predict(mach, Xnew) : predict(mach, (data = Xnew, idx = idx))
	end

	"""
	    sr_equations(mach)

	Return (eq_strings, best_idx), being robust to minor report format changes.
	"""
	function sr_equations(mach)
	    r = report(mach)
	    eqs =
	        hasproperty(r, :equation_strings) ? r.equation_strings :
	        hasproperty(r, :equations)        ? r.equations        :
	        String[]
	    best_idx = hasproperty(r, :best_idx) ? r.best_idx : (isempty(eqs) ? 1 : 1)
	    return (; eqs, best_idx, report=r)
	end

	"""
	    sr_symbolic(mach, Xnames; idx=nothing)

	Convert a chosen SR tree to a SymbolicUtils expression.
	"""
	function sr_symbolic(mach, Xnames::Vector{Symbol}; idx::Union{Nothing,Int}=nothing)
	    r = report(mach)
	    i = isnothing(idx) ? r.best_idx : idx
	    tree = r.equations[i]                      # ← not r.trees
	    try
	        # newer versions accept variable_names as a keyword
	        return node_to_symbolic(tree; variable_names=String.(Xnames))
	    catch
	        # fallback for versions expecting only the tree
	        return node_to_symbolic(tree)
	    end
	end

end


# ╔═╡ c5dfcd1d-29de-48ba-b292-2933f3d0f0e6
begin

	"""
	    sr_predict_multi(mach, Xnew; idx=nothing)
	
	Predict all targets. If `idx` is a Vector{Int}, it chooses a specific
	Pareto member *per* target (same length as number of targets).
	"""
	sr_predict_multi(mach, Xnew; idx=nothing) =
	    isnothing(idx) ? predict(mach, Xnew) :
	                     predict(mach, (data=Xnew, idx=idx))
	
	"""
	    sr_equations_multi(mach)
	
	Returns a vector of equation strings per target (and best indices).
	"""
	function sr_equations_multi(mach)
	    r = report(mach)
	    eqs =
	        hasproperty(r, :equation_strings) ? r.equation_strings :
	        hasproperty(r, :equations)        ? r.equations : Vector{Vector{String}}()
	    best =
	        hasproperty(r, :best_idx) ? r.best_idx :
	        (length(eqs)==0 ? Int[] : ones(Int, length(eqs)))
	    return (; eqs, best, report=r)
	end
	
	"""
	    sr_symbolic_multi(mach, Xnames; idxs=nothing)
	
	Convert chosen trees to SymbolicUtils expressions (one per target).
	If `idxs` omitted, uses each target’s best.
	"""
	function sr_symbolic_multi(mach, Xnames::Vector{Symbol}; idxs::Union{Nothing,Vector{Int}}=nothing)
	    r = report(mach)
	    eq_lists = hasproperty(r, :equations) ? r.equations : nothing
	    eq_lists === nothing && return nothing
	    ids = isnothing(idxs) ? r.best_idx : idxs
	    [node_to_symbolic(eq_lists[t][ids[t]]; variable_names=String.(Xnames))
	        for t in eachindex(eq_lists)]
	end
end


# ╔═╡ 580b4506-f37d-43c6-a75d-56086e03f9b9
#=╠═╡
begin
	sr_model_mt(; niterations::Int=300, seed::Int=0) = MultitargetSRRegressor(
	    niterations = niterations,
	    binary_operators = BIN_OPS,
	    unary_operators  = UN_OPS,
	    seed = seed,
	    turbo = true,
	)
	
	"""
	    sr_fit_multi(df; ys, x, model=sr_model_mt(), drop_missings=true)
	Train multitarget SR on `df` and return a NamedTuple.
	"""
	function sr_fit_multi(df::DataFrame; ys::Vector{Symbol}, x::Vector{Symbol},
	                      model = sr_model_mt(), drop_missings::Bool=true)
	    data = select(df, vcat(x, ys))
	    data = drop_missings ? dropmissing(data) : data
	    X = select(data, x); Y = select(data, ys)
	    mach = machine(model, X, Y); fit!(mach)
	    (; mach, report=report(mach), Xnames=x, Ynames=ys)
	end
	
	print_best_equations(r, ys) = for (i,t) in pairs(ys)
	    best = r.best_idx[i]
	    println("Best for $t  (idx=$(best))  =>  ", r.equation_strings[i][best])
	end
end
  ╠═╡ =#

# ╔═╡ 77e0847b-d179-4f91-8c2b-ca1e407895ba
"""
    sr_predict_on_grid(sr_m, Xgrid; idx=r.best_idx)

Use your trained multitarget SR `sr_m.mach` to predict on `Xgrid`.
`idx` can be a Vector{Int} (one best index per target) or a single Int.
Returns the prediction table (same columns as `ys_m`).
"""
sr_predict_on_grid(sr_m, Xgrid; idx) = sr_predict_multi(sr_m.mach, Xgrid; idx=idx)

# ╔═╡ 50b5bb13-afd8-4ad8-8ec5-8fd2bbd3b050
"""
    plot_sr_contour(xvals, yvals, Z; title="", xlabel="", ylabel="",
                    levels=10, colormap=:viridis, add_lines=true, labels=false)

Simple Makie contourf + optional isolines and colorbar.
Z must be a Matrix with size (length(yvals), length(xvals)).
"""
function plot_sr_contour(xvals, yvals, Z;
        title::AbstractString="", xlabel::AbstractString="", ylabel::AbstractString="",
        levels=10, colormap=:viridis, add_lines::Bool=true, labels::Bool=false)

    xu = ustrip.(xvals)          # -> Vector{Float64}
    yu = ustrip.(yvals)
    Zu = ustrip.(Z)              # -> Matrix{Float64}

    fig = Figure(resolution=(900, 650))
    ax  = Axis(fig[1,1]; title=title, xlabel=xlabel, ylabel=ylabel)

    contourf!(ax, xu, yu, Zu; levels=levels, colormap=colormap)
    #add_lines && contour!(ax, xu, yu, Zu; levels=levels, color=:black, linewidth=0.8, labels=labels)
    #Colorbar(fig[1,2], ax, label=title)
    fig
end

# ╔═╡ 39f35083-8fbd-4c1b-8e3b-796752dc8850
md"### Example: contour of predicted `Q_cond` over (ΔTSF, f_)"

# ╔═╡ c4355205-1929-48ec-9582-5b3d2b9255d3
#=╠═╡
begin
	# choose which two inputs to vary, and which target to visualize
	xcol, ycol  = :ΔTSF, :T_wall           # vary ΔTSF (x) and f_ (y)
	ztarget     = :ṁ_ref               # visualize predicted condenser heat
	
	# build the grid (hold others fixed at their medians)
	Xgrid, xvals, yvals = make_param_grid(dfu, xs_m;
	    xcol=xcol, ycol=ycol, nx=80, ny=80,
	    xscale=:identity,                # or :log10 to make x log-spaced
	    fixed=:median)                   # or fixed=Dict(:T_wall=>293.15, :L_BHE=>200.0, ...)
	
    # ⬇️ enforce column order before prediction
    Xgrid_ord = select(Xgrid, xs_m)

    Ygrid_hat = sr_predict_on_grid(sr_m, Xgrid_ord; idx=r.best_idx)

    ny, nx = length(yvals), length(xvals)
    Z = reshape(Ygrid_hat[!, ztarget], ny, nx)

    xl = string(xcol); yl = string(ycol); ttl = string(ztarget)*" (SR prediction)"

    fig_contour = plot_sr_contour(xvals, yvals, Z;
        title=ttl, xlabel=xl, ylabel=yl,
        levels=12, colormap=:viridis, add_lines=true, labels=false)
    fig_contour
end

  ╠═╡ =#

# ╔═╡ 106e37c3-8c25-4fae-b906-00b86831e58d


# ╔═╡ d6386d9d-fc6b-495c-aed3-9030073281a4


# ╔═╡ c0287358-04ff-49b5-926d-75c10a97da44


# ╔═╡ a2456a78-9ef7-4834-a535-a03c47ea735c


# ╔═╡ d56e6679-1745-4efe-8205-aab080a2c0b8
function test_pump_phx()
    # ——— Instantiate components ———
    @named pump = Hydraulics.DynamicPropBasicPump(
        H_shutoff = 12.0,      # m
        Q_max     = 13/3600,   # m³/s
        fluid     = "MEA"
    )

	@named phx = Hydraulics.DynamicPropPlateHeatExchanger(
        fluid = "MEA"
    )

	
	@named bhe = Hydraulics.BHE_Upipe_Hydraulic_Complete(depth = 100.0, Di = 0.032)

    # ——— Unknown state variable ———
    @variables ṁ        # mass flow [kg/s]

    # ——— Simple algebraic equations ———
    eqs = [
		
		phx.L ~ 0.38,
		phx.W ~ 0.07,
		phx.b ~ 1.6e-3,
		phx.np ~ 15,
		phx.ϕ ~ 1.118,
		phx.β ~ 30*π/180,
		
		phx.props.T ~ 293.15,  # 20°C
        # fix inlet temperature to 10 °C
        pump.props.T      ~ 273.15 + 10,
		bhe.T_leg1 ~ 303.15, 
		bhe.T_bend ~ 303.15, 
		bhe.T_leg2 ~ 303.15,

        # fix pump speed to 50 Hz
        #pump.pump.N       ~ 50.0,
        ṁ   ~ .8,

        # tie the flow through both devices to the same ṁ
        pump.base.ṁ      ~ ṁ,
        phx.base.ṁ        ~ ṁ,
        bhe.ṁ ~ ṁ,

		
        # enforce Δp_pump = Δp_phx
        pump.base.Δp      ~ phx.base.Δp + bhe.Δp
    ]

    # ——— Build & solve the system ———

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [ṁ], [], name = :test), [pump, phx, bhe])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])
    solver = NewtonRaphson(; 
        linesearch = BackTracking(; maxstep=20), 
        autodiff = AutoFiniteDiff()
    )
    sol = NonlinearSolve.solve(prob, solver; abstol = 1e-6)

    # ——— Report ———
    println("=== Pump–HX Steady–State ===")
    println("Converged: ", sol.retcode)
    println("Mass flow rate ṁ = ", round(sol[ṁ], digits=4), " kg/s")
    println("Pump Δp = ", round(sol[pump.base.Δp], digits=2), " Pa")
    println("HX   Δp = ", round(sol[phx.base.Δp], digits=2), " Pa")
    println("HX   Δp = ", round(sol[bhe.Δp], digits=2), " Pa")
    println("Power = ", round(sol[pump.base.P_elec], digits=2), " W")
end

# run the test

# ╔═╡ de9faed9-2f88-44e7-94da-098238b99f09
test_pump_phx()

# ╔═╡ 7960a423-b601-4aec-b0ba-b35eebce82c5
md" # Sizing components"

# ╔═╡ 2d8d7be2-6826-4258-8e19-cc96ddf79a32
low_T_sink = 273.15+20

# ╔═╡ 4f40cdce-ff84-49a6-bd17-56d2595e6866
high_T_source = 273.15+10

# ╔═╡ 83439a01-bafd-4a82-a751-2e25caa0828e
# ╠═╡ disabled = true
#=╠═╡
sol_size = size_ev()
  ╠═╡ =#

# ╔═╡ ab30e6b5-8ac9-4ae0-a855-299eee9e9543
"""
Test function for SinglePhasePlateHX_WithProperties
"""
function size_ev()
    @named hx = Evaporator()

    eqs = [
		
		hx.b ~ 0.002,
		hx.np ~ 25,#*factor,
		hx.W ~ 0.10,
		hx.ϕ ~1.2,
		hx.β ~ π/6,
		hx.th ~ 0.5e-3,
		hx.λ_pl ~ 16,
		hx.λ_pitch ~ 4e-3,
		hx.χ ~ 0.28,
		#hx.L ~ 0.2,
		
        #hx.γṁ_1   ~ log(.4),
        hx.γṁ_2   ~ log(0.025),
        hx.T_in_1 ~ 273.15 + 10,
        hx.T_out_1 ~ hx.T_in_1 - 3,
        #hx.h2i ~ hx.ev.propref.hb,
        hx.h2i ~ 225.0e3,
        hx.T_out_2 ~ hx.T_evap + 7,
		hx.logp_2 ~ log(2e5),
    ]

    @named sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], []; name=:test_plate), [hx])
    simplified = structural_simplify(sys)
	#return simplified
    prob = NonlinearProblem(simplified, [])

	solver = NewtonRaphson(; 
		linesearch = BackTracking(; maxstep=20), 
		autodiff = AutoFiniteDiff()
	)

	solver = RobustMultiNewton(
        autodiff = AutoFiniteDiff())
	sol = NonlinearSolve.solve(prob, solver, abstol = 1e-5)
	#sol = NonlinearSolve.solve(prob, abstol = 1e-6)

    println("=== Plate HX with Properties Test ===")
	
    println("Return code: ", sol.retcode)
    println("Q = ", sol[hx.Q], " W")
    println("Q  ev = ", sol[hx.ev.Q], " W")
    println("Q sh = ", sol[hx.sh.Q], " W")
    println("Tin brine= ", sol[hx.T_in_1]-273.15, " C")
    println("T mid brine= ", sol[hx.T_mid_1]-273.15, " C")
    println("T out brine= ", sol[hx.T_out_1]-273.15, " C")
    println("T evaporation = ", sol[hx.sh.propref.Tsat]-273.15, " C")
    println("T out ev = ", sol[hx.sh.T_out_2]-273.15, " C")
    println("h in = ", sol[hx.h2i], " J/kg")
    println("hmid = ", sol[hx.h2mid], " J/kg")
    println("h out = ", sol[hx.h2o], " J/kg")
    println("=== Length fraction ===")
    println("L = ", sol[hx.L], " m")
    println("EV = ", (1-zap(sol[hx.uz_fr_sh]))*100, " %")
    println("SH = ", zap(sol[hx.uz_fr_sh])*100, " %")
    println("=== HTC ===")
    println("EV h1 = ", sol[hx.ev.h_1], " W/m2/K")
    println("EV h2 = ", sol[hx.ev.h_2], " W/m2/K")
    println("SH h1 = ", sol[hx.sh.h_1], " W/m2/K")
    println("SH h2 = ", sol[hx.sh.h_2], " W/m2/K")
    println("=== Flow ===")
    println("p2 = ", exp(sol[hx.logp_2]), " Pa")
    println("p2 = ", sol[hx.sh.propref.p], " Pa")
	
    println("ṁ 1 = ", exp(sol[hx.γṁ_1]), " kg/s")
    println("ṁ 2 = ", exp(sol[hx.γṁ_2]), " kg/s")
    return sol
end

# ╔═╡ e563b30e-38b6-46db-9a74-184dd0a9117c
md" # Quick analysis: parametric study"

# ╔═╡ efea74fa-e2f8-4e68-944c-60ac7dabd571
# ╠═╡ disabled = true
#=╠═╡
begin
	m_Y = [median(y.particles) for y in YY]
	q_05 = [quantile(y.particles, 0.025) for y in YY]
	q_95 = [quantile(y.particles, 0.975) for y in YY]
end
  ╠═╡ =#

# ╔═╡ ef2eec7f-a94f-4ce1-90fe-2190efc78a5d
#=╠═╡
begin
fig_study = Makie.lines(XX, m_Y)
	band!(XX, q_05, q_95)
	fig_study
end
  ╠═╡ =#

# ╔═╡ 37e3be87-dff9-456d-b921-e0e1ecbb4909

TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
DynamicExpressions = "a40a106e-89c9-4ca8-8020-a735e8728b6b"
DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
MLJ = "add582a8-e3ab-11e8-2d5e-e98b27df1bc7"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
NonlinearSolve = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
QuasiMonteCarlo = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SharedArrays = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
SymbolicRegression = "8254be44-1295-4e6a-a16d-46603ac705cb"

[compat]
AlgebraOfGraphics = "~0.11.3"
CairoMakie = "~0.15.3"
Colors = "~0.13.1"
DataFrames = "~1.7.0"
Distributions = "~0.25.120"
DynamicExpressions = "~1.10.3"
DynamicQuantities = "~1.8.0"
Interpolations = "~0.15.1"
LaTeXStrings = "~1.4.0"
LoopVectorization = "~0.12.172"
MLJ = "~0.20.9"
ModelingToolkit = "~10.8.0"
MonteCarloMeasurements = "~1.4.7"
NonlinearSolve = "~4.9.0"
PlutoDevMacros = "~0.9.0"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.68"
ProgressMeter = "~1.10.4"
QuasiMonteCarlo = "~0.3.3"
Roots = "~2.2.8"
SpecialFunctions = "~2.5.1"
SymbolicRegression = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "fe87b14dbcedb0222b92eb2d8233096787340081"

[[deps.ADTypes]]
git-tree-sha1 = "be7ae030256b8ef14a441726c4c37766b90b93a3"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.15.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.ARFFFiles]]
deps = ["CategoricalArrays", "Dates", "Parsers", "Tables"]
git-tree-sha1 = "678eb18590a8bc6674363da4d5faa4ac09c40a18"
uuid = "da404889-ca92-49ff-9e8b-0aa6b4d38dc8"
version = "1.5.0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AlgebraOfGraphics]]
deps = ["Accessors", "Colors", "DataAPI", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "Isoband", "KernelDensity", "Loess", "Makie", "NaturalSort", "PlotUtils", "PolygonOps", "PooledArrays", "PrecompileTools", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "b24e6d3bd76b6811b195557fe6325c95e010dd4e"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.11.3"
weakdeps = ["DynamicQuantities", "Unitful"]

    [deps.AlgebraOfGraphics.extensions]
    AlgebraOfGraphicsDynamicQuantitiesExt = "DynamicQuantities"
    AlgebraOfGraphicsUnitfulExt = "Unitful"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgCheck]]
git-tree-sha1 = "f9e9a66c9b7be1ad7372bbd9b062d9230c30c5ce"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "9606d7832795cbef89e06a550475be300364a8aa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.19.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "120e392af69350960b1d3b89d41dcc1d66543858"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "b5bb4dc6248fde467be2a863eb8452993e74d402"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "1.1.1"

    [deps.Atomix.extensions]
    AtomixCUDAExt = "CUDA"
    AtomixMetalExt = "Metal"
    AtomixOpenCLExt = "OpenCL"
    AtomixoneAPIExt = "oneAPI"

    [deps.Atomix.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    OpenCL = "08131aa3-fb12-5dee-8b74-c09406e224a2"
    oneAPI = "8f75cd03-7ff8-4ecb-9b8f-daf728133b1b"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra"]
git-tree-sha1 = "26f41e1df02c330c4fa1e98d4aa2168fdafc9b1f"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.4"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "0b3c211ec48050496075523ea7ee045db9341a42"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.1"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "291532989f81db780e435452ccb2a5f902ff665f"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "1.7.0"

    [deps.BlockArrays.extensions]
    BlockArraysAdaptExt = "Adapt"
    BlockArraysBandedMatricesExt = "BandedMatrices"

    [deps.BlockArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a9014924595b7a2c1dd14aac516e38fa10ada656"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "22c63a303e8d6ff1a02a0933d579dbdb6e0f2925"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

    [deps.CategoricalArrays.extensions]
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStructTypesExt = "StructTypes"

    [deps.CategoricalArrays.weakdeps]
    JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SentinelArrays = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
    StructTypes = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"

[[deps.CategoricalDistributions]]
deps = ["CategoricalArrays", "Distributions", "Missings", "OrderedCollections", "Random", "ScientificTypes"]
git-tree-sha1 = "926862f549a82d6c3a7145bc7f1adff2a91a39f0"
uuid = "af321ab8-2d2e-40a6-b165-3d674595d28e"
version = "0.1.15"

    [deps.CategoricalDistributions.extensions]
    UnivariateFiniteDisplayExt = "UnicodePlots"

    [deps.CategoricalDistributions.weakdeps]
    UnicodePlots = "b8865327-cd53-5732-bb35-84acbb429228"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "06ee8d1aa558d2833aa799f6f0b31b30cada405f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.2"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "062c5e1a5bf6ada13db96a4ae4749a4c2234f521"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.9"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "e771a63cc8b539eca78c85b0cabd9233d6c8f06f"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "8010b6bb3388abe68d95743dcbea77650bb2eddf"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.3"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "351d6f4eaf273b753001b2de4dffb8279b100769"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.9.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "3a3dfb30697e96a440e4149c8c51bf32f818c0f3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.17.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "e215ba0e9a9e9377f2ed87cf3eb26840c8990585"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.3"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "5620ff4ee0084a6ab7097a27ba0c19290200b037"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.4"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "a86af9c4c4f33e16a2b2ff43c2113b2f390081fa"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.4.5"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "52af5ee5af4eb6ca03b3782bf65c1e5fc3024c86"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.178.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseMooncakeExt = "Mooncake"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "Functors", "LinearAlgebra", "Markdown", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "80a782f3e65d4900dcf5f2cb71f5e19d9459c04a"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.8.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "516d553f5deee7c55b2945b5edf05b6542837887"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.24.1"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "c8d85ecfcbaef899308706bebdd8b00107f3fb43"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.6.54"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DispatchDoctor]]
deps = ["MacroTools", "Preferences"]
git-tree-sha1 = "fc34127e78323c49984e1a146d577d0f890dd2b4"
uuid = "8d63f2c5-f18a-4cf2-ba9d-b3f60fc568c8"
version = "0.4.26"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.DispatchDoctor.extensions]
    DispatchDoctorChainRulesCoreExt = "ChainRulesCore"
    DispatchDoctorEnzymeCoreExt = "EnzymeCore"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "a0c11553b6c58bfc0283d8279dd3798150b77ac5"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.12"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicDiff]]
deps = ["Compat", "DispatchDoctor", "DynamicExpressions", "ForwardDiff"]
git-tree-sha1 = "fc058cee42edd2880ee1cb1d8eba753924c65363"
uuid = "7317a516-7a03-4707-b902-c6dba1468ba0"
version = "0.2.1"

[[deps.DynamicExpressions]]
deps = ["ChainRulesCore", "DispatchDoctor", "Interfaces", "MacroTools", "PrecompileTools", "Random", "Reexport", "TOML"]
git-tree-sha1 = "9f9a7092218985fc68ea66ff4aa1896b06d97757"
uuid = "a40a106e-89c9-4ca8-8020-a735e8728b6b"
version = "1.10.3"

    [deps.DynamicExpressions.extensions]
    DynamicExpressionsBumperExt = "Bumper"
    DynamicExpressionsLoopVectorizationExt = "LoopVectorization"
    DynamicExpressionsOptimExt = "Optim"
    DynamicExpressionsSymbolicUtilsExt = "SymbolicUtils"
    DynamicExpressionsZygoteExt = "Zygote"

    [deps.DynamicExpressions.weakdeps]
    Bumper = "8ce10254-0962-460f-a3d8-1f77fea1446e"
    LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"
    SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "98c4bb95af37e5d980129261fdd6dab0392c6607"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.2"

[[deps.DynamicQuantities]]
deps = ["DispatchDoctor", "PrecompileTools", "TestItems", "Tricks"]
git-tree-sha1 = "44ec2bcde862031866a9f43ee477eaa1ddb0cccc"
uuid = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
version = "1.8.0"

    [deps.DynamicQuantities.extensions]
    DynamicQuantitiesLinearAlgebraExt = "LinearAlgebra"
    DynamicQuantitiesMeasurementsExt = "Measurements"
    DynamicQuantitiesScientificTypesExt = "ScientificTypes"
    DynamicQuantitiesUnitfulExt = "Unitful"

    [deps.DynamicQuantities.weakdeps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    ScientificTypes = "321657f4-b219-11e9-178b-2701a2544e81"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EarlyStopping]]
deps = ["Dates", "Statistics"]
git-tree-sha1 = "98fdf08b707aaf69f524a6cd0a67858cefe0cfb6"
uuid = "792122b4-ca99-40de-a6bc-6742525f08b6"
version = "0.3.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EnzymeCore]]
git-tree-sha1 = "8272a687bca7b5c601c0c24fc0c71bff10aafdfd"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.12"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "eaa040768ea663ca695d442be1bc97edfe6824f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "6.1.3+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "0a2e5873e9a5f54abb06418d57a8df689336a660"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.2"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastPower]]
git-tree-sha1 = "5f7afd4b1a3969dc34d692da2ed856047325b06e"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.3"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerMooncakeExt = "Mooncake"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FeatureSelection]]
deps = ["MLJModelInterface", "ScientificTypesBase", "Tables"]
git-tree-sha1 = "4c51a6d75ffe689d4ddbf0ff136f9f80abef70ad"
uuid = "33837fe5-dbff-4c9e-8c2f-c5612fe2b8b6"
version = "0.2.3"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FindFirstFunctions]]
git-tree-sha1 = "670e1d9ceaa4a3161d32fe2d2fb2177f8d78b330"
uuid = "64ca27bc-2ba2-4a57-88aa-44e436879224"
version = "1.4.1"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "f089ab1f834470c525562030c8cfde4025d5e915"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.27.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "273bd1cd30768a2fddfa3fd63bbc746ed7249e5f"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "f88e0ba1f6b42121a7c1dfe93a9687d8e164c91b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.5"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "8e233d5167e63d708d41f87597433f59a0f213fe"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.4"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "0f265264b9287a19715dc5d491dbe3aff00c1e71"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.5.0"
weakdeps = ["GeometryBasics", "Makie", "RecipesBase"]

    [deps.GeoInterface.extensions]
    GeoInterfaceMakieExt = ["Makie", "GeometryBasics"]
    GeoInterfaceRecipesBaseExt = "RecipesBase"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"
weakdeps = ["GeoInterface"]

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "35fbd0cefb04a516104b8e183ce0df11b70a3f1a"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.3+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "dc6bed05c15523624909b3953686c5f5ffa10adc"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.ImplicitDiscreteSolve]]
deps = ["DiffEqBase", "OrdinaryDiffEqCore", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "1d49fec0d4d325f7ce56e6bd1fa7cbed24a54064"
uuid = "3263718b-31ed-49cf-8a0f-35a466e8af96"
version = "0.1.3"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
git-tree-sha1 = "8594fac023c5ce1ef78260f24d1ad18b4327b420"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.4"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interfaces]]
git-tree-sha1 = "331ff37738aea1a3cf841ddf085442f31b84324f"
uuid = "85a1e053-f937-4924-92a5-1367d23b7b87"
version = "0.3.2"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Random", "RoundingEmulator"]
git-tree-sha1 = "79342df41c3c24664e5bf29395cfdf2f2a599412"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.36"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IterationControl]]
deps = ["EarlyStopping", "InteractiveUtils"]
git-tree-sha1 = "e663925ebc3d93c1150a7570d114f9ea2f664726"
uuid = "b3c1a2ee-3fec-4384-bf48-272ea71de57c"
version = "0.5.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JuliaFormatter]]
deps = ["CommonMark", "Glob", "JuliaSyntax", "PrecompileTools", "TOML"]
git-tree-sha1 = "56b382cd34b1a80f63211a0b009461915915bf9e"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "2.1.2"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "6ac9e4acc417a5b534ace12690bc6973c25b862f"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.3"

[[deps.JuliaSyntax]]
git-tree-sha1 = "937da4713526b96ac9a178e2035019d3b78ead4a"
uuid = "70703baa-626e-46a2-a12c-08ffd08c73b4"
version = "0.4.10"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "f8da88993c914357031daf0023f18748ff473924"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.16.1"
weakdeps = ["FastBroadcast"]

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "MacroTools", "PrecompileTools", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "38a03910123867c11af988e8718d12c98bf6a234"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.37"
weakdeps = ["EnzymeCore", "LinearAlgebra", "SparseArrays"]

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b94257a1a8737099ca40bc7271a8b374033473ed"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.1"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "4f34eaabe49ecb3fb0d58d6015e32fd31a733199"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.8"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LatinHypercubeSampling]]
deps = ["Random", "StableRNGs", "StatsBase", "Test"]
git-tree-sha1 = "825289d43c753c7f1bf9bed334c253e9913997f8"
uuid = "a5e1c1ea-c99a-51d3-a14d-a9a37257b02d"
version = "1.9.0"

[[deps.LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "76627adb8c542c6b73f68d4bfd0aa71c9893a079"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.2"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LearnAPI]]
git-tree-sha1 = "c276ed13346953f3035a6d35519ec816b820cd44"
uuid = "92ad9a40-7767-427a-9ee6-6e577f1266cb"
version = "1.0.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "22e90c33c5297d6162ee54e2383584849379aa53"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.23.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "f749e7351f120b3566e5923fefdf8e52ba5ec7f9"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.4"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "e5afce7eaf5b5ca0d444bcb4dc4fd78c54cbbac0"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.172"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.LossFunctions]]
deps = ["Markdown", "Requires", "Statistics"]
git-tree-sha1 = "3e91117bd9da5e68f1fbf3d08808226ae71981b7"
uuid = "30fc2ffe-d236-52d8-8643-a9d8f7c094a7"
version = "1.0.1"
weakdeps = ["CategoricalArrays"]

    [deps.LossFunctions.extensions]
    LossFunctionsCategoricalArraysExt = "CategoricalArrays"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MLCore]]
deps = ["DataAPI", "SimpleTraits", "Tables"]
git-tree-sha1 = "73907695f35bc7ffd9f11f6c4f2ee8c1302084be"
uuid = "c2834f40-e789-41da-a90e-33b280584a8c"
version = "1.0.0"

[[deps.MLFlowClient]]
deps = ["Dates", "FilePathsBase", "HTTP", "JSON", "ShowCases", "URIs", "UUIDs"]
git-tree-sha1 = "9abb12b62debc27261c008daa13627255bf79967"
uuid = "64a0f543-368b-4a9a-827a-e71edb2a0b83"
version = "0.5.1"

[[deps.MLJ]]
deps = ["CategoricalArrays", "ComputationalResources", "Distributed", "Distributions", "FeatureSelection", "LinearAlgebra", "MLJBalancing", "MLJBase", "MLJEnsembles", "MLJFlow", "MLJIteration", "MLJModels", "MLJTuning", "OpenML", "Pkg", "ProgressMeter", "Random", "Reexport", "ScientificTypes", "StatisticalMeasures", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "acc37bedb0a1eec793829581448e100707bef3cf"
uuid = "add582a8-e3ab-11e8-2d5e-e98b27df1bc7"
version = "0.20.9"

[[deps.MLJBalancing]]
deps = ["MLJBase", "MLJModelInterface", "MLUtils", "OrderedCollections", "Random", "StatsBase"]
git-tree-sha1 = "f707a01a92d664479522313907c07afa5d81df19"
uuid = "45f359ea-796d-4f51-95a5-deb1a414c586"
version = "0.1.5"

[[deps.MLJBase]]
deps = ["CategoricalArrays", "CategoricalDistributions", "ComputationalResources", "Dates", "DelimitedFiles", "Distributed", "Distributions", "InteractiveUtils", "InvertedIndices", "LearnAPI", "LinearAlgebra", "MLJModelInterface", "Missings", "OrderedCollections", "Parameters", "PrettyTables", "ProgressMeter", "Random", "RecipesBase", "Reexport", "ScientificTypes", "Serialization", "StatisticalMeasuresBase", "StatisticalTraits", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "9fa62ebd16bfa86a3f7821a9ec49d28ad8fb6e74"
uuid = "a7f614a8-145f-11e9-1d2a-a57a1082229d"
version = "1.8.2"
weakdeps = ["StatisticalMeasures"]

    [deps.MLJBase.extensions]
    DefaultMeasuresExt = "StatisticalMeasures"

[[deps.MLJEnsembles]]
deps = ["CategoricalArrays", "CategoricalDistributions", "ComputationalResources", "Distributed", "Distributions", "MLJModelInterface", "ProgressMeter", "Random", "ScientificTypesBase", "StatisticalMeasuresBase", "StatsBase"]
git-tree-sha1 = "84a5be55a364bb6b6dc7780bbd64317ebdd3ad1e"
uuid = "50ed68f4-41fd-4504-931a-ed422449fee0"
version = "0.4.3"

[[deps.MLJFlow]]
deps = ["MLFlowClient", "MLJBase", "MLJModelInterface"]
git-tree-sha1 = "508bff8071d7d1902d6f1b9d1e868d58821f1cfe"
uuid = "7b7b8358-b45c-48ea-a8ef-7ca328ad328f"
version = "0.5.0"

[[deps.MLJIteration]]
deps = ["IterationControl", "MLJBase", "Random", "Serialization"]
git-tree-sha1 = "ad16cfd261e28204847f509d1221a581286448ae"
uuid = "614be32b-d00c-4edb-bd02-1eb411ab5e55"
version = "0.6.3"

[[deps.MLJModelInterface]]
deps = ["REPL", "Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "66626f80d5807921045d539b4f7153b1d47c5f8a"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.11.1"

[[deps.MLJModels]]
deps = ["CategoricalArrays", "CategoricalDistributions", "Combinatorics", "Dates", "Distances", "Distributions", "InteractiveUtils", "LinearAlgebra", "MLJModelInterface", "Markdown", "OrderedCollections", "Parameters", "Pkg", "PrettyPrinting", "REPL", "Random", "RelocatableFolders", "ScientificTypes", "StatisticalTraits", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "09381923be5ed34416ed77badbc26e1adf295492"
uuid = "d491faf4-2d78-11e9-2867-c94bc002c0b7"
version = "0.17.9"

[[deps.MLJTuning]]
deps = ["ComputationalResources", "Distributed", "Distributions", "LatinHypercubeSampling", "MLJBase", "ProgressMeter", "Random", "RecipesBase", "StatisticalMeasuresBase"]
git-tree-sha1 = "38aab60b1274ce7d6da784808e3be69e585dbbf6"
uuid = "03970b2e-30c4-11ea-3135-d1576263f10f"
version = "0.8.8"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "MLCore", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "a772d8d1987433538a5c226f79393324b55f7846"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.8"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "96d73e05b6f3079df0963b66c2844d162263a896"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.4"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "6e64d2321257cc52f47e193407d0659ea1b2b431"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.5"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.ModelingToolkit]]
deps = ["ADTypes", "AbstractTrees", "ArrayInterface", "BlockArrays", "ChainRulesCore", "Combinatorics", "CommonSolve", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "DiffRules", "DifferentiationInterface", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "DynamicQuantities", "EnumX", "ExprTools", "FindFirstFunctions", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "Graphs", "ImplicitDiscreteSolve", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "Latexify", "Libdl", "LinearAlgebra", "MLStyle", "Moshi", "NaNMath", "NonlinearSolve", "OffsetArrays", "OrderedCollections", "OrdinaryDiffEqCore", "PrecompileTools", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SCCNonlinearSolve", "SciMLBase", "SciMLPublic", "SciMLStructures", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "URIs", "UnPack", "Unitful"]
git-tree-sha1 = "b08ce77612d4f32f08a77f7d00cc8e512a0ef153"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "10.8.0"

    [deps.ModelingToolkit.extensions]
    MTKBifurcationKitExt = "BifurcationKit"
    MTKCasADiDynamicOptExt = "CasADi"
    MTKDeepDiffsExt = "DeepDiffs"
    MTKFMIExt = "FMI"
    MTKInfiniteOptExt = "InfiniteOpt"
    MTKLabelledArraysExt = "LabelledArrays"
    MTKPyomoDynamicOptExt = "Pyomo"

    [deps.ModelingToolkit.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    CasADi = "c49709b8-5c63-11e9-2fb2-69db5844192f"
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"
    FMI = "14a09403-18e3-468f-ad8a-74f8dda2d9ac"
    InfiniteOpt = "20393b10-9daf-11e9-18c9-8db751c92c57"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    Pyomo = "0e8e1daf-01b5-4eba-a626-3897743a3816"

[[deps.MonteCarloMeasurements]]
deps = ["Distributions", "GenericSchur", "LinearAlgebra", "MacroTools", "Random", "SLEEFPirates", "StaticArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "de98932a1016cbd6c926ed76b3e740a46b153c70"
uuid = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
version = "1.4.7"

    [deps.MonteCarloMeasurements.extensions]
    ForwardDiffExt = "ForwardDiff"
    MakieExt = "Makie"
    MeasurementsExt = "Measurements"
    RecipesBaseExt = "RecipesBase"
    UnitfulExt = "Unitful"

    [deps.MonteCarloMeasurements.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "fade91fe9bee7b142d332fc6ab3f0deea29f637b"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.9"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Random", "ScopedValues", "Statistics"]
git-tree-sha1 = "4abc63cdd8dd9dd925d8e879cda280bedc8013ca"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.30"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"
    NNlibFFTWExt = "FFTW"
    NNlibForwardDiffExt = "ForwardDiff"
    NNlibSpecialFunctionsExt = "SpecialFunctions"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "aeb6fb02e63b4d4f90337ed90ce54ceb4c0efe77"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.9.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "ee395563ae6ffaecbdf86d430440fddc779253a4"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.13.0"

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseDiffEqBaseExt = "DiffEqBase"
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"

    [deps.NonlinearSolveBase.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LineSearch = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
    LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "65101a20b135616a13625ae6f84b052ef5780363"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.6.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "3e04c917d4e3cd48b2a5091b6f76c720bb3f7362"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.7.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "DiffEqBase", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "3398222199e4b9ca0b5840907fb509f28f1a2fdc"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.2.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenML]]
deps = ["ARFFFiles", "HTTP", "JSON", "Markdown", "Pkg", "Scratch"]
git-tree-sha1 = "63603b2b367107e87dbceda4e33c67aed17e50e0"
uuid = "8b6db2d4-7670-4922-a472-f9537c81ab66"
version = "0.3.2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "87510f7292a2b21aeff97912b0898f9553cc5c2c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "c1f51f704f689f87f28b33836fd460ecf9b34583"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.11.0"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "1bd20b621e8dee5f2d170ae31631bf573ab77eec"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.26.2"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.PlutoDevMacros]]
deps = ["JuliaInterpreter", "Logging", "MacroTools", "Pkg", "TOML"]
git-tree-sha1 = "1cb861c9295d79dc6e23170d4b33bce013f69643"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.9.1"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ec9e63bd098c50e4ad28e7cb95ca7a4860603298"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.68"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "Random"]
git-tree-sha1 = "bb178012780b34046c6d1600a315d8dbee89d83d"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.5"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.PrettyPrinting]]
git-tree-sha1 = "142ee93724a9c5d04d78df7006670a93ed1b244e"
uuid = "54e16d92-306c-5ea0-a30b-337be88ac337"
version = "0.4.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "13c5103482a8ed1536a54c08d0e742ae3dca2d42"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.4"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.QuasiMonteCarlo]]
deps = ["Accessors", "ConcreteStructs", "LatticeRules", "LinearAlgebra", "Primes", "Random", "Requires", "Sobol", "StatsBase"]
git-tree-sha1 = "cc086f8485bce77b6187141e1413c3b55f9a4341"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.3.3"
weakdeps = ["Distributions"]

    [deps.QuasiMonteCarlo.extensions]
    QuasiMonteCarloDistributionsExt = "Distributions"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "dbe5fd0b334694e905cb9fda73cd8554333c46e2"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.1"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "4dd1a95cc16d5abdccc4eac5faf6bc73904be1a2"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.35.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "668e411c0616a70860249b4c96e5d35296631a1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.8"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SCCNonlinearSolve]]
deps = ["CommonSolve", "PrecompileTools", "Reexport", "SciMLBase", "SymbolicIndexingInterface"]
git-tree-sha1 = "5595105cef621942aceb1aa546b883c79ccbfa8f"
uuid = "9dfe8606-65a1-4bb3-9748-cb89d1561431"
version = "1.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "c9dc4c04bcb0146a35dd6af726073c5738b80e3b"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.104.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "3414071e3458f3065de7fa5aed55283b236b4907"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.8"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "3249fe77f322fe539e935ecb388c8290cd38a3fc"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.3.1"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.ScientificTypes]]
deps = ["CategoricalArrays", "ColorTypes", "Dates", "Distributions", "PrettyTables", "Reexport", "ScientificTypesBase", "StatisticalTraits", "Tables"]
git-tree-sha1 = "4d083ffec53dbd5097a6723b0699b175be2b61fb"
uuid = "321657f4-b219-11e9-178b-2701a2544e81"
version = "3.1.0"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "7f44eef6b1d284465fafc66baf4d9bdcc239a15b"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.4.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.ShowCases]]
git-tree-sha1 = "7f534ad62ab2bd48591bdeac81994ea8c445e4a5"
uuid = "605ecd9f-84a6-4c9e-81e2-4798472b76a3"
version = "0.1.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "09d986e27a606f172c5b6cffbd8b8b2f10bf1c75"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.7.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveDiffEqBaseExt = "DiffEqBase"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "9de43e0b9b976f1019bf7a879a686c4514520078"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.21"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "cbea8a6bd7bed51b1619658dec70035e07b8502f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.14"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.StatisticalMeasures]]
deps = ["CategoricalArrays", "CategoricalDistributions", "Distributions", "LearnAPI", "LinearAlgebra", "MacroTools", "OrderedCollections", "PrecompileTools", "ScientificTypesBase", "StatisticalMeasuresBase", "Statistics", "StatsBase"]
git-tree-sha1 = "d15f71761dce0a2833d458824e6b60edbfd4f138"
uuid = "a19d573c-0a75-4610-95b3-7071388c7541"
version = "0.2.1"
weakdeps = ["LossFunctions", "ScientificTypes"]

    [deps.StatisticalMeasures.extensions]
    LossFunctionsExt = "LossFunctions"
    ScientificTypesExt = "ScientificTypes"

[[deps.StatisticalMeasuresBase]]
deps = ["CategoricalArrays", "InteractiveUtils", "MLUtils", "MacroTools", "OrderedCollections", "PrecompileTools", "ScientificTypesBase", "Statistics"]
git-tree-sha1 = "e4f508cf3b3253f3eb357274fe36fb3332ca9896"
uuid = "c062fc1d-0d66-479b-b6ac-8b44719de4cc"
version = "0.1.2"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "89f86d9376acd18a1a4fbef66a56335a3a7633b8"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.5.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "b81c5035922cc89c2d9523afc6c54be512411466"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.5"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "9022bcaa2fc1d484f1326eaa4db8db543ca8c66d"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.4"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f35f6ab602df8413a50c4a25ca14de821e8605fb"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "8ad2e38cbb812e29348719cc63580ec1dfeb9de4"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.1"
weakdeps = ["Adapt", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrettyTables", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "59ca6eddaaa9849e7de9fd1153b6faf0b1db7b80"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.42"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicRegression]]
deps = ["ADTypes", "Compat", "ConstructionBase", "Dates", "DifferentiationInterface", "DispatchDoctor", "Distributed", "DynamicDiff", "DynamicExpressions", "DynamicQuantities", "LineSearches", "Logging", "LossFunctions", "MLJModelInterface", "MacroTools", "Optim", "Pkg", "PrecompileTools", "Printf", "ProgressMeter", "Random", "Reexport", "SpecialFunctions", "StatsBase", "StyledStrings", "TOML"]
git-tree-sha1 = "11582404a6e7ffcf889a9b2699784edc780ec38e"
uuid = "8254be44-1295-4e6a-a16d-46603ac705cb"
version = "1.9.0"

    [deps.SymbolicRegression.extensions]
    SymbolicRegressionEnzymeExt = "Enzyme"
    SymbolicRegressionJSON3Ext = "JSON3"
    SymbolicRegressionSymbolicUtilsExt = "SymbolicUtils"

    [deps.SymbolicRegression.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
    SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper", "WeakValueDicts"]
git-tree-sha1 = "fa63e8f55e99aee528951ba26544403b09645979"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.29.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "8922f50f13eac24583a445968d3cfbdcc0c621ac"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.45.0"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestItems]]
git-tree-sha1 = "42fd9023fef18b9b78c8343a4e2f3813ffbcefcb"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.0.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "02aca429c9885d1109e58f400c333521c13d48a0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.4"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Transducers]]
deps = ["Accessors", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "SplittablesBase", "Tables"]
git-tree-sha1 = "7deeab4ff96b85c5f72c824cae53a1398da3d1cb"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.84"

    [deps.Transducers.extensions]
    TransducersAdaptExt = "Adapt"
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d2282232f8a4d71f79e85dc4dd45e5b12a6297fb"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.23.1"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "b13c4edda90890e5b04ba24e20a310fbe6f249ff"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.3.0"

    [deps.UnsafeAtomics.extensions]
    UnsafeAtomicsLLVM = ["LLVM"]

    [deps.UnsafeAtomics.weakdeps]
    LLVM = "929cbde3-209d-540e-8aea-75f648917ca0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "4ab62a49f1d8d9548a1c8d1a75e5f55cf196f64e"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.71"

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╠═6f32ff1e-5faa-4551-a898-01167021e4a4
# ╠═f5ba6b5c-9f00-43a2-ae87-4d86444c4b3f
# ╠═9539a581-42e3-43dc-996d-4d9316aa9c92
# ╠═8527ed5d-3ff6-435d-a93d-548c824dee78
# ╠═eadd49a1-7bca-40a2-abdb-064a0fb710e6
# ╠═cedbf94b-e85d-4662-a64e-16fba9ef2ed4
# ╠═d47bf527-6960-4e14-b9e5-a422f8888e66
# ╠═910c36c9-4ec1-4b0f-b044-823a229fa536
# ╠═00bcb75f-8da5-46a8-9c20-6b94861342c4
# ╠═1499b12a-8a87-4354-b4b5-d181ffd63f1b
# ╠═fa83b49d-d4ad-470a-88f8-524ae0e3191f
# ╠═3ef17b0f-7ba5-4cef-98fe-2abf3711c907
# ╠═00996e13-d7a2-4f9e-be84-4eec83428486
# ╠═2aa4e5cf-f6d9-4608-8000-ec078c7ae7b9
# ╠═bbaed5c6-445a-43ed-82f2-c7b4b17fbf19
# ╠═fe20b40a-ce2a-408d-89c4-3bd17255ad62
# ╠═f0fa9e95-4369-480d-a6d4-d6b545504686
# ╠═68df8cdc-5652-4fd7-9391-b9de65722f10
# ╠═6e948a61-d8bb-4d93-b6e6-da19256cb24f
# ╠═f3d29aa1-ee85-46a9-8258-f6b9f03d87a2
# ╠═164a20af-4421-4cb6-a4fe-0f8b00ad7059
# ╠═a2cd3b4a-713b-4d5d-8739-139d0a64214c
# ╠═f3070a9a-d98b-49b8-b300-631a4065211e
# ╠═734d74db-7b19-42fd-89ad-78816a939565
# ╠═e31e4221-54e9-4049-b80e-eca118e9bb00
# ╠═eacfb867-b8fe-4166-84df-038848a9d837
# ╠═2ab2c592-3e29-4031-8d19-8856f3cdbc6e
# ╠═004dd778-a212-4147-a9c0-ad60f1eb539b
# ╟─63cf7ccc-e17f-48d7-bbf9-2419a6ebe336
# ╟─6d536fb0-e23f-4d51-8839-8d9dc2422cbd
# ╠═96465cb7-82db-4287-bdb7-9260d39e12ef
# ╠═917d72da-5985-42b2-a079-9dd14b97a791
# ╠═b2fdf520-379d-426a-8271-7170ee7c7ea5
# ╠═1e443538-f5b0-4de7-8ec1-df8737e5ba70
# ╠═e4c92fbe-af30-48e9-b501-1ab71bc0d182
# ╠═39697eef-1b94-46ac-9d75-943997cbb913
# ╠═43489277-7150-4ad1-aabc-b76bd5472b8d
# ╠═18475dd7-225b-48d5-b88a-81f4223b8dd8
# ╠═e9829f7f-cfd0-4a4a-ae4c-00ff405cca42
# ╟─0c1298e5-440c-467b-b334-a4e0929dd83d
# ╠═86df6265-c57f-4ef8-bc0c-522957653fba
# ╠═53fb2d9b-13ea-44c5-9cff-03d9562dbfb1
# ╠═2f2fa303-81a3-4bc4-b8a5-caa21851bf23
# ╠═6a83b5bc-4e25-4781-afa8-8b378d6c838e
# ╠═87178d51-5926-4ae4-aea4-cbfc211d7200
# ╠═010a0e27-2e68-427b-9fef-569a32590d87
# ╠═11df4133-c9f6-436a-9b2b-ada02f872bf0
# ╠═5818a465-6685-40da-b0eb-a7e035c0170c
# ╟─d70ebc3d-7082-4043-925d-220f0ac753f7
# ╟─65e682c7-af7f-43ec-aa70-e12553026888
# ╠═62cea410-f0af-4219-a7df-92a0e574d74d
# ╠═2ed52c9e-a7ec-4793-a4c4-4a1a1611636f
# ╠═a1eec254-6b86-440d-a2a4-2151d588dc23
# ╠═e51c75dc-07f3-4682-9ce2-ed2e72977ea5
# ╠═d31ad2dc-528d-46f6-a19f-90658a1615d0
# ╠═19a9664b-cfa8-4a1a-b2bb-ed3f8c357100
# ╠═70ca56b4-dcde-44a2-aa6c-cd5b84c28c16
# ╠═ce7c4e4d-e0ef-4b49-8ae8-0f63acfe4dc4
# ╠═6be93b5a-be7d-47ad-81ef-72e3d7350fd6
# ╠═e438b966-15eb-42e4-a610-0fc426385bc7
# ╠═96966964-d666-4dc2-a554-8c359e3db528
# ╠═b41466ba-fcff-4def-957b-2b96032643aa
# ╠═9b117e8b-7a38-49ba-9586-f3aacaafa0f3
# ╠═143ccb02-56b0-4b50-a1e5-0d06b53f1b46
# ╠═98c20212-3447-481c-b0d4-326a1703997c
# ╠═f7eeabcd-23e6-4733-81aa-a6bbde82837e
# ╠═9d502c24-7e38-4c3f-96ed-0e18425d82d9
# ╟─0907ee2b-5ce8-4789-936a-ec14a4c3c059
# ╠═039f80ce-620b-4863-aa98-1f0f2ed49a7e
# ╠═eb386e80-6a11-4b38-b55a-50bbad670474
# ╠═07f01adf-8cc7-4a40-aa14-a90cef491b95
# ╠═11823168-8288-49db-88e8-319df4535273
# ╠═409aa2af-c6b7-451c-8557-ced24eccc9c1
# ╠═51de15da-535e-4908-b346-d69d34e2ce01
# ╠═9b9c3b3b-6e9e-49c0-8498-99217be15b2e
# ╠═fcfeffb2-35aa-4209-a495-ce79949d08e0
# ╠═93114c1c-393c-4569-bcad-936ad8938e44
# ╠═86eeb20e-7da2-489d-84af-317f31f63ef4
# ╠═7d1b4554-571a-4113-afdd-243a24260b7c
# ╠═d01db8eb-17dc-405e-b566-63fa0e86682b
# ╠═9418798c-4603-498a-ab7a-a45d509aa67f
# ╠═42bda7b4-ea1c-4717-a514-a6d97aa68f93
# ╠═1bd62dd3-e84e-4e0b-8731-890cc33286b8
# ╠═2d05d7cf-fa80-4552-a2ca-a5ceac3939f1
# ╠═6e831997-cf76-4ea4-815f-8f5ffd48d3e1
# ╠═e8abe6ea-9eda-4f83-abdc-da2bd91dd26a
# ╠═6e3b84d1-afe1-4c56-98ba-50b2f13cfc46
# ╠═c81e49af-a294-45d2-b2fd-c605dad70ef3
# ╠═1adba409-b8cb-430d-9ab6-c29aad55035f
# ╠═fddf235e-b540-4bc9-99d6-b75a7c49bdf6
# ╠═6466aac2-8c05-4bea-a0bc-d2c4723529e0
# ╠═192d94f7-1dba-4ccf-8a14-7f1e090e84f5
# ╠═e4f83718-32dd-4170-bfae-36c5edd28997
# ╠═e5614c61-6260-414f-b100-3a3683de0524
# ╠═832139e7-83e3-496e-8e84-3499c4b9c04b
# ╠═e09caba7-c987-4c7b-af32-beec786c2504
# ╠═fbfcb718-1d8d-462f-8c9d-c7e340152d50
# ╠═3d841b73-85ed-40dd-b78a-4c85dc471164
# ╠═93a62e8e-200c-401f-bcbf-31d4fa7c9275
# ╠═7ac5bbd6-d4a9-4b17-bd0e-b473d97c8a57
# ╠═b12af1c4-53f5-49ed-94f9-dd9f8434bcf0
# ╠═5eae6b7f-ee19-4f7e-8943-e8828ef027fb
# ╠═f66a54da-5889-4dfb-a093-5a8eb8964b42
# ╠═d2b14ef0-9b53-4846-a169-27de59cb9f0f
# ╠═d3cc6b23-4a6d-451c-8a6e-81bf1b71abfb
# ╠═d930fddf-92b8-468c-a556-f0126a93eb43
# ╠═b5d46480-3583-4eff-aeec-861e2b1bd1cb
# ╠═757853b5-e844-432a-9538-b7afab10bb72
# ╠═b9a9a9a0-dca2-45ee-8250-24bb680d839e
# ╠═e1c53c82-3ede-40a7-ae3f-e848fc763372
# ╠═f60e1452-e13c-472d-85bb-fbb64c3e1311
# ╠═364eaab7-4812-4467-8319-777f0b12d1c7
# ╟─819264a6-e4b7-406c-9d97-8ac0f36c3ea8
# ╟─80885b7b-9ab2-4010-a954-ed4790c277f0
# ╠═3ad72686-9a35-4486-a71c-180cba7d46a2
# ╟─241265f4-f66f-4158-855f-89f002f23d4e
# ╟─52e4fb9f-6e5e-4c7e-b34e-9f0bace4db63
# ╠═9aa6b70a-9464-4f48-a340-61030ca9310b
# ╠═5f7f03fb-20e2-4bd2-8bb7-f09a8c2d94e6
# ╠═2a4c04e6-875e-453a-9be5-6194a5530c06
# ╟─22e21ab4-f5a9-4ccd-8fdf-280279e12096
# ╠═8d3a070d-682b-47f4-bfd6-9feb53ce9ccc
# ╠═1ac95f91-934a-4dc7-97a8-f6cc13ecb864
# ╠═df4eb258-955d-4f9e-bfdd-6d202c5b8d90
# ╠═79b5e151-be4a-4d38-ae46-c975b4bf557e
# ╠═04ec6d0d-d812-464e-b244-8458207283da
# ╠═a0e018cb-a75a-4844-8716-98492b29b6d9
# ╠═b9a1f8e0-7ca4-4435-b2fa-741b0a25da82
# ╠═fbc3e326-5684-4029-ba49-6c0f065d151b
# ╠═7ef1715e-2fd3-4deb-853d-233496735474
# ╠═88624925-6092-44c8-9539-28527aa43a7c
# ╟─526d0dad-939f-48af-a09f-c1cae30d6770
# ╟─f6aed83e-02a5-4d62-874e-55f333266b42
# ╟─34d568ce-25d5-4a0b-b1db-0ca3b8e8607a
# ╟─76cc998b-9a7f-4b1c-8403-c16bccbda74e
# ╟─5bf0b015-e0f5-4e49-baa8-c34ee5d75a73
# ╠═d7090f5b-50f6-4198-9237-8447ae442789
# ╠═fdddfe24-f560-4a52-9756-ee1d1d6946d3
# ╠═104492fc-1c63-4571-80e4-2df35b6abdcf
# ╠═0ab5ba9b-ee1d-4e7b-b16d-29ab59563ff2
# ╠═a4b78e57-3f4a-4bd1-a4fd-344aed782596
# ╠═58a5d982-5d07-4a63-950a-ca5639cff43f
# ╠═89dc461f-71b0-4b91-ab47-a70d2929fec2
# ╠═ceb2b429-90ae-4bc5-b346-e65e5727d441
# ╠═80898f01-7945-48d9-abd5-15de3c76d465
# ╠═b3bb8c8b-3825-4a21-94f7-a2d0d00fbb09
# ╠═a496bd64-5a9a-4b3f-bb8f-1a0567943d1c
# ╠═6ac1ecbf-4833-449f-9933-c782a8566fb0
# ╠═8b6e0e2e-ecc5-475e-bd31-1f23893dcb38
# ╠═892ef213-3229-4d91-9ba7-742127f366c8
# ╠═0ddc1362-ee46-452b-a5d7-6d03df36e31d
# ╠═30d6e371-c921-4fc2-8d8f-a6c5d7c34385
# ╠═19556fe1-b85d-4da1-ac62-a56b5a72618b
# ╠═af94b9d0-44ca-42e4-9e8c-96f5864bc80d
# ╠═8a4bc20a-195f-49f5-8379-ca6be9bdb1a1
# ╠═a208968e-8f31-4a82-b8f3-9db4a8d27e8b
# ╠═4b47374a-f1a3-402c-a4ba-bae8edcedc01
# ╠═b48efffc-d19c-4b66-936b-6fd9c9de5a88
# ╠═37d9352d-2357-4f36-b651-5dc73246fd42
# ╠═bfe4471c-7033-4ca5-beb1-530a8be5b757
# ╠═670cd198-239a-4240-be5d-2def55383a36
# ╠═0f6b83f7-167a-4e3f-9fd5-d4f4695fbaea
# ╠═798e29fe-16e5-4f68-baad-9edc88509245
# ╠═05bbe96e-023c-413b-a347-d37e74fad7fe
# ╠═fecd7f3d-bf0b-464e-a1ad-238e508b0f38
# ╟─0fe18091-e8f6-42d4-9961-5783452ddefd
# ╟─7d25aec4-c657-410a-8ae4-03195cefbee7
# ╠═5575055e-9136-4b4e-925a-ddf82ca1156c
# ╠═4412cbfc-d470-4003-a2bd-606df1be7022
# ╠═f100b275-2eb2-4d10-93e2-8954f314ae9e
# ╠═f8b2fee6-eb35-493c-ac30-331176d9a795
# ╠═058c533c-104d-4ae9-80e3-f356787d0d13
# ╠═4899159f-1c86-478c-b874-ef8db1808bc5
# ╠═c8044a03-837b-411c-8919-7d47f83111c3
# ╠═b220d228-e3cb-4dd2-843e-03db76f46fe9
# ╠═d7384891-c5c0-4614-99cf-0f1cf58f759c
# ╠═9b2f80c6-18d0-4b10-827a-b9a665f7f1c4
# ╠═86429d81-3516-488f-be3d-df2d83bdbae8
# ╠═4982df2c-016a-4136-a074-b39c1befffb7
# ╠═c18c39fe-6896-433e-90d4-4a8149ae9b5b
# ╠═ce0ce754-c26a-4b8d-81d9-8234a165a0c6
# ╠═bed552fa-1a51-492b-8552-81186c96c6a1
# ╟─473f4d2b-92a3-495c-975d-38738f714fb9
# ╠═c2755bcf-548b-47ff-be58-b6cbfda319f7
# ╠═afba0093-bdb2-4ee9-95da-0f4d1933d22c
# ╠═bdbf26be-8a8b-4c4d-be59-a67ea68226ef
# ╠═85978ad5-c80a-4eed-b04f-72b57afcc13e
# ╠═2ce9d158-4156-47db-97bc-efa0e18618b3
# ╠═aa8498ba-14aa-421d-abb8-67ac1286bd41
# ╠═0761ab79-e4d9-43c2-b772-fe6df8f333b7
# ╠═b92cc8df-5884-4129-99fe-d22038b4eb0c
# ╠═1e116b39-5019-4e87-80eb-13b1a9d9ae58
# ╠═fb29ede8-6f49-4e14-a6d9-d8caee90e5ba
# ╠═870d81bf-2254-427b-8471-3c9d23b13a72
# ╠═02497be0-7d2f-43d1-b365-d4a8d72f045e
# ╠═fc5900c1-9704-4384-9a29-be545faf74d6
# ╠═8262db7e-b994-4afa-b823-661a0bb4949b
# ╠═eb8cd9d6-9237-4723-84d3-3f293d74d1e9
# ╠═0b7c9c92-e263-4a91-bbf9-d32cb0ac7024
# ╠═6eab6576-10af-432d-abf6-42f424541e5b
# ╠═8dea2e5f-82b5-408e-843d-c25b7d265d48
# ╠═b153b909-a8b2-4e5d-abe7-a3c1bad81358
# ╠═2be85177-c03c-48de-8d4f-c2053fcf2a13
# ╠═24b94de6-1d03-4225-8fdd-95cac5beb33c
# ╠═a24d0373-b198-42f6-a6bb-c2ede7259244
# ╠═74455461-5277-4fee-a53b-edf34e9e26f4
# ╠═dc49f979-5f43-4234-864d-acfef268f5f9
# ╠═50af6ba3-363d-4ca6-977a-7eb35919b27b
# ╠═88cbfed8-9de9-4421-8dd0-8f42b82234bd
# ╠═9035962a-af30-4870-8ccb-c340892e81c3
# ╠═b6a370b1-bfe7-4324-b5a6-6018ddf55e12
# ╠═d7b463f3-4a9d-4572-9ae6-eb5d4ea47531
# ╠═60d5e58e-237b-4f13-ba00-c50a5315ffd7
# ╠═b9879207-d8b7-4465-b899-e7656e48e5fc
# ╠═09b7509b-a308-49a4-966d-33d5ef7f9cb1
# ╠═1945625f-482a-4eda-985d-7112b4cf432d
# ╠═7386981b-7839-4588-aa68-5c4b3d461b06
# ╠═977b3df3-0fa4-452f-884b-94b24ec54cf3
# ╠═1320fac5-2cf2-484d-a01e-98ecd559588d
# ╠═a3fe0b30-007d-415d-af04-f5712e3d0489
# ╠═c21e1161-99f5-4d15-adb8-7f7c669eadb0
# ╠═298f0a4b-d6e7-4772-9862-54895ba7aab8
# ╟─ea33f723-8f17-4000-b483-eb8f1f9354f3
# ╠═daed3493-da2b-489c-b7ae-5bcd66bfc03b
# ╠═b72f0888-4716-48ec-8a99-fc068a988abf
# ╠═e086552e-9ce0-4353-a16b-395d2986d1a9
# ╠═dfc85639-a09b-44fa-99af-61706a192549
# ╠═dc7866d3-68e8-40de-b82d-26717fff53d9
# ╠═5e034021-1116-4d1d-8034-d582c8a9643f
# ╠═e6a2deea-3658-42b3-8ef9-83e31be5a300
# ╠═99e84662-a903-424d-b270-d665a52f479f
# ╠═4217415e-772c-44d2-8c3d-2e85cf01fb1e
# ╠═21a3e541-8ad9-41de-995a-e2b43221d15a
# ╠═1ce3b2c5-71ba-4a1f-ab48-961585e57943
# ╠═5d5f3c09-164b-4a55-9bf1-b2bb90221465
# ╠═4175c8ba-346c-4afa-9267-ed7e689cff35
# ╠═629aaec1-eca1-4b91-910d-5579b49c8df2
# ╠═16928040-8730-4f0c-a67c-50f9d1ae852b
# ╟─712bf240-09c1-4afb-be0d-5d76f9e5bbb1
# ╠═05d571f0-c990-430e-971c-965ba527e814
# ╠═e174ba39-dd8c-4af9-88c1-caa504af2c7b
# ╠═d9a2618d-e923-4370-a20d-b8faf44bf06c
# ╠═387e3115-866a-4831-bb59-a18edf749375
# ╠═6ae2052e-907c-49b9-bb51-e2e3b0c641bc
# ╟─5ce8337c-aae0-476f-914b-854a054b26db
# ╠═9fa61b9f-1d33-444f-9c1e-d7ed70933000
# ╠═a0332504-e854-4b16-bb9d-e5f4ca277f40
# ╠═3b86a38d-9b62-4db0-b943-a1b973eb9b3c
# ╠═a966eea4-f3e4-46f9-8230-7ab19eb7f905
# ╠═a7e7ede3-1caa-4483-91f3-c26393f2b45a
# ╠═e2504f4a-6606-4051-b91e-189bbd218edf
# ╠═fdfb7456-6d8e-44c8-a883-1d1faeb96ede
# ╟─2f8b1672-934b-4b93-9e05-6169ba730cc8
# ╟─9beedddb-8b5a-4811-98cb-9a26ecd77b52
# ╠═5eed01c3-561f-490b-80d5-23db0846bfe3
# ╠═6ebdb8db-d5cc-4d02-a471-1e29a29ab42b
# ╠═1a9a9e58-41ac-45fb-9891-2e2e531621e5
# ╠═97b66208-6b79-40dd-bfd9-b710807a4e3f
# ╠═8bc61318-43a1-447d-81a5-d2fed06799b2
# ╠═0bc8ed9a-6d4d-4983-a204-9b54ef7a1fa9
# ╠═bba28e41-70b6-409a-a4ec-1a1246261e82
# ╟─8c67c3c7-1452-42cf-a7f6-0a40b5bdd071
# ╠═960c408c-3040-4a00-b285-c01080cb1c02
# ╠═644901ea-7778-4e7d-a5cf-44e34df77827
# ╠═57e985e5-023c-4389-8e08-2185bed0751a
# ╠═eac9124f-ffa0-4a9a-948d-3c6e5ec0c7a2
# ╠═bb90f897-17d3-4a9d-a462-fd6d3f4d98e7
# ╠═6acbadd7-94f6-4399-ac5f-3a2724ce641c
# ╠═29e4a45f-7a0b-425b-90d5-b2fc84f43d10
# ╠═4e1204d9-eab8-4890-a957-92341c8edfc4
# ╟─987290e0-8b32-481a-85a0-442b075463cb
# ╠═6a624b9b-1e84-4e98-b1a8-60fd8e95cb13
# ╠═1f4b84eb-7ca1-478e-8e4f-376ac069ea1e
# ╠═8ddb2096-a2b0-4369-a188-5d972af91428
# ╠═195ab453-b1f2-4632-9420-091f7f702ed7
# ╠═517540bb-548e-4b97-8e96-bd674f8a9d18
# ╟─44e95eb8-ef96-4832-9168-df5c504b0569
# ╠═8a1007ab-a556-4670-a2a6-18292d721c89
# ╠═cf7f3e02-3ac6-495d-a6a8-81711ae0fe57
# ╠═6011e1fa-1058-4bce-bb65-fcdf25951bfa
# ╠═0d69490f-83ee-46fc-b034-8370e8a66361
# ╠═6d9c073c-e05b-45ec-a4b4-1b47c43cc5ba
# ╠═f28598ab-9f11-4bd5-b5e4-6e4ae2aea350
# ╠═941dcb04-88be-43ae-bc21-ae3ef9dd0db0
# ╠═6b0ca393-a200-4ee0-b71e-ea3f3b68c3f4
# ╠═880e6614-2556-4d3a-9408-f59b45ba94fe
# ╠═592a2cfb-9d8a-481e-a39d-3147dc76b643
# ╠═c2a86042-1d4b-498c-8370-0d38f73a2bfd
# ╟─c5875432-f3f0-4495-a0f0-853d87d73353
# ╠═acbfba23-bbb5-4d9e-93b7-1120677332bd
# ╠═b0b5c52c-a795-40a0-9c8d-f678010c454a
# ╠═08d88b40-980e-4d17-91d1-2b1a06180bdd
# ╠═f0e43670-bbd9-487a-b8b6-8b5d6804039c
# ╠═d7d4f78c-a845-42ae-86ac-1fbca7c3cfcb
# ╠═09952ea8-ef5d-414f-ac15-478ad8ea6343
# ╠═16868428-1547-4b8d-a900-5e5120f38088
# ╠═62b6a403-53ac-42c1-9fff-055a9b7a478b
# ╠═849d55d9-3e76-4ccb-8635-c8b70376d0b1
# ╠═34ed86a4-20e9-4b2b-9344-7ec6dfaab815
# ╠═41fe5c48-3b28-4d41-8969-e3d59a84cbfa
# ╠═54809c43-f7c3-493b-b06a-a13bf0a44bc6
# ╠═bf3a0c4a-5fea-465d-a3ca-e742b74dc78f
# ╠═c8c35ed3-bd6c-40ac-8049-b603b003f603
# ╠═7ceca66f-f75a-4ba5-aa44-1f2b6939da77
# ╠═72ec6c2e-498a-47c1-a15f-31d1e5cd14c7
# ╠═33c1996a-aade-400a-9b2f-259db12c0835
# ╠═79b17009-14ee-422e-8eb6-9e81c270b602
# ╟─94138b2d-0fd6-494e-9fd1-ef5aecadaabf
# ╟─2097d037-7f79-4fe2-b683-b74b77d98752
# ╠═599e27a4-1aba-4de4-a76d-fb09511de993
# ╠═4b5d8e48-6775-4ec6-ab28-3e96b9954af3
# ╠═c9aa37be-53b4-4f23-9898-b5179f7117f2
# ╠═bb0f32e5-f72f-41d1-af47-ed1cbb571bb1
# ╠═3ec7080c-2972-484b-805e-3d9fe0a7d9d1
# ╠═7148e421-b7e8-4f3c-a0a9-da3d14dd32e5
# ╠═23f9f68a-9e3c-4bbe-8abe-1b851206162f
# ╠═7520ab74-810f-431a-b24e-724e24595b13
# ╠═1d8f15e6-0462-4f86-9a97-3364e3e73bc0
# ╠═3afecaa8-1c77-4b26-b534-2c419ddfb027
# ╠═2164a2d3-59a0-4573-bd12-b62431f46c04
# ╠═318a7ca3-1d6a-4a32-9fc7-c2a0072c92d3
# ╠═c1333a2c-28a8-4a31-a049-8479e43f6c0e
# ╠═1177eea2-1d56-4090-8b70-a636f404b4fb
# ╠═1d782694-8f7d-4d1f-a1bc-b7646ac01b23
# ╠═a35e92ea-bee1-4117-9ccf-f2016be29dc3
# ╠═17aa8588-7089-4247-aaad-0ca54e1d888a
# ╠═286269bd-f4db-4690-81f4-b4c45edeb24d
# ╠═34522b07-e28b-4fa6-b2bf-2c7c73895aa2
# ╠═9d0ffbe1-deff-4de9-a0d7-9e7c8df5bfaf
# ╠═bf82a17d-f524-4aa7-ab17-eadcc045d6cb
# ╟─b570d453-43f2-4a12-a9d6-889752dd5aeb
# ╠═21cd24f0-4d8b-4298-89e1-7ba7817c2eef
# ╟─fff91fa6-8ffb-458a-946d-244b5b13c3cd
# ╠═698de8c1-b7a4-477c-ba42-88e78cabddf2
# ╠═dfdcded3-026e-4925-ac02-ba801488b267
# ╠═1abb42cf-af03-49a9-837b-902999a89744
# ╠═51d3286a-d89d-41f5-88d5-bb361e76ff09
# ╠═5e9e8e04-4aca-416d-a4b4-410b243fd722
# ╟─47ab6638-1135-452c-ae79-bb5682094a51
# ╠═4ef8beb3-621c-41f4-97cf-5a634ec87769
# ╠═011f212e-63b1-4df3-a646-37bcf23717f8
# ╠═12ca0c64-931b-4989-9594-6bd09ed6fea9
# ╠═0aa504e7-2687-4fae-a9a0-3614455e6db4
# ╠═0ef59446-cc78-47f8-8196-32fdfe199980
# ╠═f6b71da4-db72-4896-9e34-b7afc2e5677c
# ╠═312118cb-b321-4d8a-b1d9-62411134a980
# ╠═23899c17-5414-405e-9bb8-7f505b45e2f6
# ╟─e319da9f-acb5-43e4-830f-01822b9a19ab
# ╟─31aa8e9b-9122-497f-866f-ef359f7604ad
# ╠═6ccf5d2b-1e32-4aa6-bea2-f1830e795280
# ╠═0e0ec7a4-31a2-4583-a8e5-0a44f48d2c59
# ╟─794b5b4d-62d7-4dda-a379-a7f1d9af7573
# ╠═36f7dc30-ec0b-493b-bdd0-092c393342fa
# ╟─e78e5fac-9201-4208-a643-3583eecae761
# ╠═a8a0c668-be27-4e64-8a03-039c3bf07df8
# ╠═92b07487-ebec-4026-989a-bf1bcdc9d1c9
# ╟─3812a920-c1f8-4ecd-ae82-8d5ab03b3522
# ╠═76482fe3-b1c0-4a40-a41a-a95362d181ee
# ╟─dd40bbe0-6b5c-4bcc-a371-f6a956b7759b
# ╠═1a0f0c34-0fe2-44ef-b67b-8ecb51732f7f
# ╠═65f2bffb-80eb-42f6-9685-bd719b0f61ed
# ╠═e7f2b4bd-7b17-41d8-9c60-a41d88181456
# ╠═3a4edf4f-d73b-4dd0-a706-da7f14c4801b
# ╠═01ee0a4e-6d90-44ea-8119-bfd48892bff1
# ╠═19f8c59d-f243-45db-837c-f664193b24a4
# ╟─e2a509cf-f4eb-4db2-9cf9-b390e3c36cd8
# ╠═ce47fdb8-b5ef-4a17-b1a7-a47240576545
# ╟─c6469779-0f69-4845-9183-719ecd3b7c26
# ╟─60a39d0c-f909-48b7-922f-eb88bc11b088
# ╟─ae544fa4-9b4a-48a4-ade9-01881147fa8f
# ╟─c5dfcd1d-29de-48ba-b292-2933f3d0f0e6
# ╟─580b4506-f37d-43c6-a75d-56086e03f9b9
# ╠═77e0847b-d179-4f91-8c2b-ca1e407895ba
# ╠═50b5bb13-afd8-4ad8-8ec5-8fd2bbd3b050
# ╟─39f35083-8fbd-4c1b-8e3b-796752dc8850
# ╠═c4355205-1929-48ec-9582-5b3d2b9255d3
# ╠═106e37c3-8c25-4fae-b906-00b86831e58d
# ╠═d6386d9d-fc6b-495c-aed3-9030073281a4
# ╠═c0287358-04ff-49b5-926d-75c10a97da44
# ╠═a2456a78-9ef7-4834-a535-a03c47ea735c
# ╟─d56e6679-1745-4efe-8205-aab080a2c0b8
# ╠═de9faed9-2f88-44e7-94da-098238b99f09
# ╠═7960a423-b601-4aec-b0ba-b35eebce82c5
# ╠═2d8d7be2-6826-4258-8e19-cc96ddf79a32
# ╠═4f40cdce-ff84-49a6-bd17-56d2595e6866
# ╠═83439a01-bafd-4a82-a751-2e25caa0828e
# ╠═ab30e6b5-8ac9-4ae0-a855-299eee9e9543
# ╠═e563b30e-38b6-46db-9a74-184dd0a9117c
# ╠═efea74fa-e2f8-4e68-944c-60ac7dabd571
# ╠═ef2eec7f-a94f-4ce1-90fe-2190efc78a5d
# ╠═e9a4105f-abec-47ca-9170-ca6aef7a9c49
# ╠═4bbe463f-7e79-44f8-8abe-5f104b346e45
# ╠═5870e251-b3da-42b7-b72f-211e74c669fc
# ╠═16436692-04a1-4915-8e89-b8c0805b6bd5
# ╠═37e3be87-dff9-456d-b921-e0e1ecbb4909
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
