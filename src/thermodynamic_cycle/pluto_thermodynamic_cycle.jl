### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

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

# ╔═╡ 734d74db-7b19-42fd-89ad-78816a939565
@fromparent begin
	using ..Properties
end

# ╔═╡ fa83b49d-d4ad-470a-88f8-524ae0e3191f
using ModelingToolkit

# ╔═╡ 3ef17b0f-7ba5-4cef-98fe-2abf3711c907
using NonlinearSolve

# ╔═╡ e9a4105f-abec-47ca-9170-ca6aef7a9c49
using MonteCarloMeasurements

# ╔═╡ 4bbe463f-7e79-44f8-8abe-5f104b346e45
using Distributions

# ╔═╡ 5870e251-b3da-42b7-b72f-211e74c669fc
using LinearAlgebra

# ╔═╡ 229b5dd1-b47b-4b67-8bc5-1694cdec228f


# ╔═╡ 43349e67-c35f-4416-a8d8-76284c25042b


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

# ╔═╡ 18f3eb6f-d15e-4992-a071-eb335930fe38
"""
    amalfi_convective_boiling_htc(D_h, λ, β, ṁ, A_ch,
                                  μₗ, kₗ, Prₗ, h_fg, q″, ρₗ, ρᵥ, x)

Compute the convective-boiling HTC h_cb [W/m²·K] for fully evaporating flow:

h_cb = Gₑ₁ * (D_h * kₗ) * Reₑq * Gₑ₂ * Prₗ^0.4 * Boₑq^0.3

where
  Gₑ₁ = 2.81 * (D_h/λ)^(-0.041) * (2π - β)^(-2.83)
  Gₑ₂ = 0.746 * (D_h/λ)^(-0.082) * (2π - β)^(0.61)

Reₑq = Gₑq * D_h / μₗ  
Boₑq = Gₑq * h_fg / q″  

Gₑq = ṁ/A_ch * [ (1 - x) + x * sqrt(ρᵥ/ρₗ) ]  

# Arguments
- D_h    Hydraulic diameter [m]  
- λ      Corrugation pitch [m]  
- β      Chevron angle [rad]  
- ṁ      Mass flow rate [kg/s]  
- A_ch   Channel flow area [m²]  
- μₗ     Liquid viscosity [Pa·s]  
- kₗ     Liquid thermal conductivity [W/m·K]  
- Prₗ    Liquid Prandtl number  
- h_fg   Latent heat [J/kg]  
- q″     Heat flux [W/m²]  
- ρₗ, ρᵥ Liquid/vapor density [kg/m³]  
- x      Vapor quality (mass fraction)

# Returns
- h_cb   Convective-boiling HTC [W/m²·K]

# References
- Amalfi, R. L., Vakili-Farahani, F., & Thome, J. R. (2016). Flow Boiling and Frictional Pressure Gradients in Plate Heat Exchangers.
https://ht.readthedocs.io/en/release/ht.boiling_plate.html
"""
function amalfi_convective_boiling_htc(D_h,
                                       λ,
                                       β,
                                       ṁ,
                                       A_ch,
                                       μₗ,
                                       kₗ,
                                       Prₗ,
                                       h_fg,
                                       q″,
                                       ρₗ,
                                       ρᵥ,
                                       x)
    # Geometry factors
    Gₑ₁ = 2.81 * (D_h/λ)^(-0.041) * (2π - β)^(-2.83)
    Gₑ₂ = 0.746 * (D_h/λ)^(-0.082) * (2π - β)^(0.61)
    # Equivalent mass flux
    Gₑq = (ṁ / A_ch) * ((1 - x) + x * sqrt(ρᵥ/ρₗ))
    # Dimensional numbers
    Reₑq = Gₑq * D_h / μₗ
    Boₑq = Gₑq * h_fg / q″
    # HTC
    return Gₑ₁ * (D_h * kₗ) * Reₑq * Gₑ₂ * Prₗ^0.4 * Boₑq^0.3
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

# ╔═╡ e31c8c25-3ca9-43cd-a865-5c1500677c24


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


# ╔═╡ 270071c2-f672-4f31-974c-ac304471b517


# ╔═╡ 1ad977dd-ebd6-464d-98c8-ab2148576012


# ╔═╡ 54d7923d-fd7b-482f-b48e-c1e6641f2572


# ╔═╡ b7ee8628-d8d7-4a93-864b-64fdd68d72b0


# ╔═╡ b2fdf520-379d-426a-8271-7170ee7c7ea5
function ExpansionValve(; name,
	fluid= "R134A",
)
	
	ps = @parameters begin
	end
	sts = @variables begin
		p_i = 3e5
		p_o = 1e5
	end
	
	@named i = Properties.RefrigerantProps(fluid = fluid)
	@named o = Properties.RefrigerantProps(fluid = fluid)
	
	eqs = [
		p_i ~ i.p,
		p_o ~ o.p,
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
		exp.i.p ~ 5e5,
		exp.o.p ~ 1e5,
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
	V = 0.02/50,
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
	
		ṁ ~ V̇ * i.ρ,
		
		V̇ ~ f*V,
		
		W ~ ṁ*(o.h-i.h),

		# Enforce the isentropic condition:
		o_is.p ~ o.p,		
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
		compressor.ṁ ~ 0.05,
		compressor.f ~ 50,
		compressor.i.h ~ 400e3,
		#compressor.i.p ~ 1e5,
		compressor.o.p ~ 5e5,
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
test_compressor()

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
        hx.C_1 ~ 1250,     # Heat capacity rate side 1 (W/K)
        hx.C_2 ~ 1500.0,     # Heat capacity rate side 2 (W/K)
        
        # Inlet temperatures
        hx.T_in_1 ~ 280.0,   # Inlet temperature side 1 (K) (hot side)
        hx.T_in_2 ~ 290.0    # Inlet temperature side 2 (K) (cold side)
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
    
    return sol
end

# ╔═╡ 5818a465-6685-40da-b0eb-a7e035c0170c
test_single_phase_hx_temps()

# ╔═╡ 65e682c7-af7f-43ec-aa70-e12553026888
md"""
## Same with properties
"""

# ╔═╡ 62cea410-f0af-4219-a7df-92a0e574d74d
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

# ╔═╡ 2ed52c9e-a7ec-4793-a4c4-4a1a1611636f
"""
    h_plate(ṁ, ρ, cp, μ, λ, b, np, W, ϕ, β)

Compute heat‐transfer coefficient in a plate channel.
- ṁ: mass‐flow rate (kg/s)
- ρ, cp, μ, λ: fluid properties
- b: half‐channel depth (m)
- np, W, ϕ: plate geometry
- β: chevron angle (rad)
"""
function h_plate(ṁ, ρ, cp, μ, λ, b, np, W, ϕ, β)
    # 1) Flow area: (np-1) channels, each W × 2b, times surface enlargement
    A_ch = (np - 1) * W * (2b) * ϕ  # [m²]

    # 2) Volumetric velocity
    v = abs(ṁ) / (ρ * A_ch)         # [m/s]

    # 3) Dimensionless numbers
    dh = 2b                         # hydraulic diameter [m]
    Re = v * dh / μ
    Pr = cp * μ / λ

    # 4) Nusselt number (your chevron correlation)
    Nu = nusselt_plate(Re, Pr, β)

    # 5) Heat‐transfer coefficient
    return Nu * λ / dh              # [W/m²·K]
end


# ╔═╡ e226a6e5-d8ba-4f33-8b6b-dfde3d8d404d


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
function SinglePhasePlateHX_WithProperties(; name, b=0.001, np=10, W=0.2, ϕ=1.1, β=π/6, th=0.0005, λ_Cu=400.0)
    # ——— Parameters ———
    pars = @parameters begin
		b = b
		np = np
		W = W
		ϕ = ϕ
		β = β
		th = th
		λ_Cu = λ_Cu
	end

    # ——— State variables ———
    vars = @variables begin
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
		NTU = 2
    end

    # ——— Subsystems ———
    @named base    = SinglePhaseHX()
    @named prop1i  = Properties.IncompressibleProps(fluid="MEA", χ=0.28)
    @named prop1o  = Properties.IncompressibleProps(fluid="MEA", χ=0.28)
    @named prop2i  = Properties.IncompressibleProps(fluid="MEA", χ=0.28)
    @named prop2o  = Properties.IncompressibleProps(fluid="MEA", χ=0.28)

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
        h_1         ~ h_plate(ṁ_1, 0.5*(prop1i.ρ+prop1o.ρ), 0.5*(prop1i.cp+prop1o.cp), 0.5*(prop1i.μ+prop1o.μ), 0.5*(prop1i.λ+prop1o.λ), b, np, W, ϕ, β),
        h_2         ~ h_plate(ṁ_2, 0.5*(prop2i.ρ+prop2o.ρ), 0.5*(prop2i.cp+prop2o.cp), 0.5*(prop2i.μ+prop2o.μ), 0.5*(prop2i.λ+prop2o.λ), b, np, W, ϕ, β),
        # overall coefficient
        U           ~ 1/(1/h_1 + th/λ_Cu + 1/h_2),
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
    @named hx = SinglePhasePlateHX_WithProperties(b=0.001, np=12, W=0.05, ϕ=1.0, β=π/6, th=0.0005, λ_Cu=400.0)

    eqs = [
        # geometry
        hx.L ~ 0.2,
        # mass flows
        hx.γṁ_1 ~ log(0.5),    # 0.1 kg/s water
        hx.γṁ_2 ~ log(0.5),   # 0.15 kg/s MEA
        # inlet temps
        hx.T_in_1 ~ 320.0,
        hx.T_in_2 ~ 290.0
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name=:test), [hx])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])
    sol = NonlinearSolve.solve(prob, NewtonRaphson())

    println("=== Plate HX with Properties Test ===")
    println("Q = ", sol[hx.Q], " W")
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


# ╔═╡ b41466ba-fcff-4def-957b-2b96032643aa
test_single_phase_plate_hx()

# ╔═╡ c1235068-16e0-401a-b8c8-a5da5cc99a6c


# ╔═╡ 7ac463ee-57d1-4d79-a9f6-dae84f35b426


# ╔═╡ 9b117e8b-7a38-49ba-9586-f3aacaafa0f3
function PHX_Superheater(; name,
                         b=0.001, np=10, W=0.2, ϕ=1.1, β=π/6,
                         th=0.0005, λ_Cu=400.0,
                         λ_pitch=0.002,      # new corrugation pitch [m]
                         Re_low=3000,        # blender limits
                         Re_high=8000)
    # ——— Parameters ———
    pars = @parameters begin
        b = b
        np = np
        W = W
        ϕ = ϕ
        β = β
        th = th
        λ_Cu = λ_Cu
        λ_pitch = λ_pitch
        Re_low = Re_low
        Re_high = Re_high
    end

    # ——— State variables ———
    vars = @variables begin
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
        NTU = 2
        # refrigerant inlet/outlet enthalpies
        h2i = 300e3
        h2o = 300e3
        pref = 1e5
    end

    # ——— Subsystems ———
    @named base    = SinglePhaseHX()
    @named prop1i  = Properties.IncompressibleProps(fluid="MEA", χ=0.28)
    @named prop1o  = Properties.IncompressibleProps(fluid="MEA", χ=0.28)
    @named propref = Properties.RefrigerantProps()

    # ——— Derived expressions ———
    ṁ_1 = exp(γṁ_1)
    ṁ_2 = exp(γṁ_2)
    dh   = 2*b

    # refrigerant quality and two‐phase inputs
    x_ref  = (propref.h - propref.hb)/(propref.hd - propref.hb)
    h_fg   = propref.hd - propref.hb
    ΔT_ref = T_out_2 - propref.Tsat
    q_flux = Q / A
	Prb = propref.cpb*propref.ηb/propref.λb

    # ——— Equations ———
    eqs = [
        # geometry mapping
        A           ~ L*np*W*ϕ,
        A           ~ base.A,

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
            b, np, W, ϕ, β
        ),

        h_2 ~ h_plate(
            ṁ_2,
            propref.ρ,
            propref.cp,
            propref.η,
            propref.λ,
            b, np, W, ϕ, β
        ),

        # overall U‐value
        U ~ 1/(1/h_1 + th/λ_Cu + 1/h_2),
        U ~ base.U,

        # heat‐transfer balances
        Q       ~ base.Q,
        Q       ~ ṁ_2*(h2o - h2i),
        propref.p ~ pref,
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
        b      = 0.001,
        np     = 2,
        W      = 0.05,
        ϕ      = 1.0,
        β      = π/6,
        th     = 0.0005,
        λ_Cu   = 400.0,
        λ_pitch= 0.002,   # corrugation pitch
        Re_low = 3000,
        Re_high= 8000
    )

    # — build a tiny set of boundary conditions to fully specify the superheater
    eqs = [
        # geometry
        hx.L         ~ 0.2,

        # mass flows (0.5 kg/s each side)
        hx.γṁ_1      ~ log(0.5),
        hx.γṁ_2      ~ log(0.5),

        # hot‐side temps (MEA side)
        hx.T_in_1    ~ 293.15,
        hx.T_out_1   ~ 313.15,

        # refrigerant side: inlet/outlet temperatures
        hx.T_in_2    ~ 280.0,
        hx.T_out_2   ~ 290.0,

        # refrigerant pressure and enthalpies
        hx.pref      ~ 2e5,       # 2 bar
        hx.h2i       ~ 200e3,     # inlet enthalpy [J/kg]
        hx.h2o       ~ 250e3,     # outlet enthalpy [J/kg]
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
    println("NTU   = ", sol[hx.base.NTU])
    println("ε     = ", sol[hx.base.ε])
    println("T₁,in = ", sol[hx.T_in_1]-273.15, " °C")
    println("T₁,out= ", sol[hx.T_out_1]-273.15, " °C")
    println("T₂,in = ", sol[hx.T_in_2]-273.15, " °C")
    println("T₂,out= ", sol[hx.T_out_2]-273.15, " °C")
    println("h2i   = ", sol[hx.h2i],       " J/kg")
    println("h2o   = ", sol[hx.h2o],       " J/kg")
    println("p_ref = ", sol[hx.pref]/1e5,  " bar")
    println("T_sat = ", sol[hx.propref.T_b]-273.15, " °C")
    println("ρ_b   = ", sol[hx.propref.ρ_b], " kg/m³")
    println("ρ_d   = ", sol[hx.propref.ρ_d], " kg/m³")
    println("h_fg  = ", sol[hx.propref.h_d - hx.propref.h_b], " J/kg")
    return sol
end


# ╔═╡ d24d1392-d47c-43c8-aee2-7001e66f4540
sol_sh = test_PHX_Superheater()

# ╔═╡ c9ceff67-86ff-4cff-8f4d-5c32b11fe2db
sol_sh.retcode

# ╔═╡ f8d39170-abda-4664-b17e-4c71e90e462a
sol_sh.resid

# ╔═╡ cafa5fd8-7d67-4cfd-82a2-bd2f6c94b737


# ╔═╡ f7eeabcd-23e6-4733-81aa-a6bbde82837e


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
        Q        ~ C_1 * (T_out_1 - T_in_1),
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
        hx.T_2  ~ 360.0,
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
    println("T_out_pc = ", sol[hx.T_2], " K (== T_in_pc)")
    println("T_out_o  = ", sol[hx.T_in_1], " K")
    println("T_out_o  = ", sol[hx.T_out_1], " K")
    return sol
end

# run test

# ╔═╡ 039f80ce-620b-4863-aa98-1f0f2ed49a7e
test_phase_change_hx()

# ╔═╡ 9fc18265-9496-4033-aed5-2717bcfd3971


# ╔═╡ 5cf6a034-9868-448e-9a6b-0c607d49194e


# ╔═╡ b2d62ab4-212a-4ddc-a8e2-45cb144715b3
"""
    X_tt(x, μℓ, μg)

Lockhart–Martinelli parameter for evaporation.
- x: vapor quality
- μℓ, μg: liquid, vapor viscosities
"""
X_tt(x, μℓ, μg) = ((1 - x)/x)^0.9 * sqrt(μg/μℓ)

# ╔═╡ d9f8989e-99c5-428c-9913-f690b313e44d
"""
    Phi2_lam(Xtt)

Two‐phase multiplier in “laminar” regime (Re_eq ≤ 2000).
"""
Phi2_lam(Xtt) = 1 + 1500/Xtt + 1/Xtt^2

# ╔═╡ ba435e92-b5cc-4f05-8b41-4608f1176c23
"""
    Phi2_turb(Xtt)

Two‐phase multiplier in “turbulent” regime (Re_eq > 2000).
"""
Phi2_turb(Xtt) = 1 + 200/sqrt(Xtt)

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

# ╔═╡ 7f4bf9ca-4444-41b9-b71c-94aa626e9ca0
"""
    h_evap_plate(ṁ, A, D_h, μℓ, μg, θ, hℓ, x)

Compute boiling HTC for a fully‐evaporating region:
1. Re_eq = (ṁ/A)*D_h/μℓ  
2. Xtt = X_tt(x, μℓ, μg)  
3. choose Φ² = Phi2_lam or Phi2_turb  
4. return htp = Φ² * hℓ
- ṁ: mass flow [kg/s]
- A: flow area [m²]
- D_h: hydraulic diameter [m]
- μℓ, μg: viscosities [Pa·s]
- θ: corrugation parameter (°)
- hℓ: single-phase HTC from Wanniarachchi [W/m²K]
- x: vapor quality
"""
function h_evap_plate_2(ṁ, A, D_h, μℓ, μg, Re_low, Re_high, hℓ, x)
	G    = ṁ / A
	Re   = G * D_h / μℓ
    Xtt  = X_tt(x, μℓ, μg)
    w    = blend_weight(Re, Re_low, Re_high)
    Φ_2   = (1 - w)*Phi2_lam(Xtt) + w*Phi2_turb(Xtt)
    return Φ_2 * hℓ
end

# ╔═╡ 11823168-8288-49db-88e8-319df4535273
function PhaseChangePlateHX_WithProperties(; name,
        b=0.001, np=10, W=0.2, ϕ=1.1, β=π/6, th=0.0005,
        λ_Cu=400.0, fluid_inc="MEA", χ_inc=0.28, refrigerant = "R134A",
        λ_pitch = 0.002,       # ➤ corrugation pitch [m]
        Re_low   = 2000.0,     # ➤ blender lower Re limit
        Re_high  = 2100.0      # ➤ blender upper Re limit
    )
    # ——— Parameters ———
    pars = @parameters begin
        b = b
        np = np
        W = W
        ϕ = ϕ
        β = β
        th = th
        λ_Cu = λ_Cu
        χ_inc = χ_inc
        λ_pitch = λ_pitch
        Re_low   = Re_low
        Re_high  = Re_high
    end

    # ——— State variables ———
    vars = @variables begin
        L = 1.0
        γṁ_2 = 0.0          # log mass flow for refrigerant
        γṁ_1 = 0.0          # log mass flow for incompressible fluid
        h_2 = 600.0         # film coefficients (init)
        h_1 = 600.0
        T_in_1 = 280.0
        T_out_1 = 280.0
        T_2    = 300.0      # refrigerant bulk T
        p_2    = 1e5        # refrigerant pressure
        h2i    = 300e3
        h2o    = 300e3
        Q      = 1000.0
    end

    # ——— Subsystems ———
    @named base   = PhaseChangeHX()
    @named pref   = Properties.RefrigerantProps(fluid = refrigerant)
    @named prop1i = Properties.IncompressibleProps(fluid = fluid_inc, χ = χ_inc)
    @named prop1o = Properties.IncompressibleProps(fluid = fluid_inc, χ = χ_inc)

    # ——— Derived expressions ———
    ṁ_1 = exp(γṁ_1)
    ṁ_2 = exp(γṁ_2)
    A    = L * np * W * ϕ       # total heat-transfer area
	D_h = 2b/ϕ

    # ——— Equations ———
    eqs = [
        # geometry
        base.A       ~ A,

        # incompressible side mapping
        T_in_1       ~ base.T_in_1,
        prop1i.T     ~ T_in_1,
        T_out_1      ~ base.T_out_1,
        prop1o.T     ~ T_out_1,

        # refrigerant side mapping
        base.T_2     ~ pref.Tsat,     # refrigerant bulk = saturation T
        T_2          ~ base.T_2,
        pref.h       ~ h2i,           # set to inlet enthalpy for sat props

        # capacity rates (incompressible only)
        base.C_1     ~ ṁ_1 * 0.5*(prop1i.cp + prop1o.cp),

		
        h_2 ~ h_evap_plate_2(
		    ṁ_2,        # refrigerant mass flow
		    W*b*np,          # flow area
		    D_h,        # hydraulic diameter (add as a parameter)
		    pref.ηb,    # μₗ at bubble
		    pref.ηd,    # μ_g at dew
		    Re_low,     # blender lower Re limit
		    Re_high,    # blender upper Re limit
		    pref.hb,         # single-phase HTC via wanniarachchi_phe
		    pref.X       # vapor quality
		),
		
        # film coeff: incompressible side unchanged
        h_1 ~ h_plate(
            ṁ_1,
            0.5*(prop1i.ρ + prop1o.ρ),
            0.5*(prop1i.cp + prop1o.cp),
            0.5*(prop1i.μ + prop1o.μ),
            0.5*(prop1i.λ + prop1o.λ),
            b, np, W, ϕ, β
        ),

        # overall U linking
        base.U       ~ 1/(1/h_1 + th/λ_Cu + 1/h_2),

        # heat-transfer balances
        Q            ~ base.Q,
        Q            ~ ṁ_2*(h2o - h2i),
    ]

    sys = NonlinearSystem(eqs, vars, pars; name = name)
    return ModelingToolkit.compose(sys, [base, pref, prop1i, prop1o]; name = name)
end


# ╔═╡ d01db8eb-17dc-405e-b566-63fa0e86682b
#=╠═╡
function test_pcpwp_display()
    @named hx = PhaseChangePlateHX_WithProperties(
        b=0.001, np=4, W=0.05, ϕ=1.0, β=π/6,
        th=0.0005, λ_Cu=400.0
    )

    eqs = [
        hx.L        ~ 0.2,
        hx.γṁ_1     ~ log(0.5),
        hx.γṁ_2     ~ log(0.5),
        hx.T_in_1   ~ 320.0,
        hx.h2i      ~ 300.0e3,
        hx.pref.p   ~ 2e5,
    ]

    @named sys   = ModelingToolkit.compose(
        NonlinearSystem(eqs, [], []; name=:test_plate),
        [hx]
    )
    simpl = structural_simplify(sys)
    prob  = NonlinearProblem(simpl, [])
    solver = NewtonRaphson(
        linesearch = BackTracking(maxstep=20),
        autodiff    = AutoFiniteDiff()
    )
    sol = NonlinearSolve.solve(prob, solver; abstol = 1e-6)

    # Extract everything
    b       = sol.ps[hx.b]
    np      = sol.ps[hx.np]
    W       = sol.ps[hx.W]
    ϕ       = sol.ps[hx.ϕ]
    β       = sol.ps[hx.β]
    λ_pitch = sol.ps[hx.λ_pitch]
    Re_low  = sol.ps[hx.Re_low]
    Re_high = sol.ps[hx.Re_high]

    ṁ1      = exp(sol[hx.γṁ_1])
    ṁ2      = exp(sol[hx.γṁ_2])
    A       = sol[hx.base.A]
    Q       = sol[hx.base.Q]

    Tsat    = sol[hx.pref.Tsat]
    ρb      = sol[hx.pref.ρb]
    ρd      = sol[hx.pref.ρd]
    h_fg    = sol[hx.pref.hd] - sol[hx.pref.hb]
    ηb      = sol[hx.pref.ηb]
    ηd      = sol[hx.pref.ηd]
    Prb     = sol[hx.pref.Prb]
    x_ref   = (sol[hx.pref.h] - sol[hx.pref.hb])/(sol[hx.pref.hd] - sol[hx.pref.hb])

    h1      = sol[hx.h_1]
    h2      = sol[hx.h_2]
    # f       = sol[hx.f]    # if you’ve also extracted the single-phase friction factor

    # Organize into a DataFrame for display
    df = DataFrame(
        Category        = [
            "Geometry", "Geometry", "Geometry", "Geometry", "Geometry", 
            "Solver Limits", "Solver Limits",
            "Flows & Area","Flows & Area","Flows & Area","Flows & Area",
            "Refrig. Props","Refrig. Props","Refrig. Props","Refrig. Props","Refrig. Props","Refrig. Props","Refrig. Props",
            "HTC"
        ],
        Name            = [
            "b", "np", "W", "ϕ", "β (deg)",
            "Re_low", "Re_high",
            "ṁ₁","ṁ₂","A","Q",
            "T_sat","ρ_b","ρ_d","h_fg","η_b","η_d","Pr_b",
            "x_ref"
        ],
        Value           = [
            b, np, W, ϕ, β*180/π,
            Re_low, Re_high,
            ṁ1, ṁ2, A, Q,
            Tsat, ρb, ρd, h_fg, ηb, ηd, Prb,
            x_ref
        ]
    )

    return sol, df
end

  ╠═╡ =#

# ╔═╡ 86eeb20e-7da2-489d-84af-317f31f63ef4
#=╠═╡
sol_pcpwp, df_sol = test_pcpwp_display()
  ╠═╡ =#

# ╔═╡ 9ff7dcd3-b011-43f2-9280-440c58cf8016
#=╠═╡
df_sol
  ╠═╡ =#

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
        b=0.001, np=10, W=0.2, ϕ=1.1, β=π/6, th=0.0005,
        λ_Cu=400.0, fluid_inc="MEA", χ_inc=0.28, refrigerant = "R134A")
    # ——— Parameters ———
    pars = @parameters begin
        b = b
        np = np
        W = W
        ϕ = ϕ
        β = β
        th = th
        λ_Cu = λ_Cu
        χ_inc = χ_inc
    end

    # ——— State variables ———
    vars = @variables begin
        L = 1.0
		
        γṁ_2 = 0.0          # log mass flow for ref
        γṁ_1 = 0.0          # log mass flow for incompressible fluid

		uz_fr_sh = 0.0 # fraction of superheater, unbounded variable R -> 0 1 interval

        #Q = 1e4
        h2i = 600.0
        h2mid = 600.0
        h2o = 600.0
        # PhaseChangeHX variables
        T_in_1 = 280.0
        T_mid_1 = 280.0
        T_out_1 = 280.0
        
		p_2 = 300.0

		Q = 3000
    end

    # ——— Derived expressions ———
	fr_sh = zap(uz_fr_sh)
	fr_ev = 1-zap(uz_fr_sh)
    ṁ_1 = exp(γṁ_1)
    ṁ_2 = exp(γṁ_2)
    dh     = 2*b

	

    @named ev = PhaseChangePlateHX_WithProperties(b=0.001, np=4, W=0.05,
                ϕ=1.0, β=π/6, th=0.0005, λ_Cu=400.0)

    @named sh = PHX_Superheater(b=0.001, np=2, W=0.05, ϕ=1.0, β=π/6, th=0.0005, λ_Cu=400.0)

    eqs = [
        # geometry
        sh.L ~ L * fr_sh,
        ev.L ~ L * fr_ev,
        # mass flows
        sh.γṁ_1 ~ γṁ_1,    # 0.1 kg/s water
        sh.γṁ_2 ~ γṁ_2,
    	ev.γṁ_1 ~ γṁ_1,    # 0.1 kg/s water
        ev.γṁ_2 ~ γṁ_2,    # 0.1 kg/s water
        # inlet temps
        sh.T_in_1 ~ T_in_1,
        sh.T_out_1 ~ T_mid_1,
        ev.T_in_1 ~ T_mid_1,
        ev.T_out_1 ~ T_out_1,
        
		sh.pref ~ p_2,
		ev.pref.p ~ p_2,
		
		sh.h2o ~ h2o,
		sh.h2i ~ h2mid,
		ev.h2o ~ h2mid,
		ev.h2i ~ h2i,
		
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
		
		#hx.uz_fr_sh ~ 0.0,
        hx.L      ~ 0.2,
        hx.γṁ_1   ~ log(0.25),
        hx.γṁ_2   ~ log(0.25),
        #hx.γṁ_2   ~ log(0.25),
        hx.T_in_1 ~ 320.0,
        hx.T_out_1 ~ hx.T_in_1 - 3,
        hx.h2i ~ 350.0e3,
		#hx.p_2 ~ 2e5,
    ]

    @named sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], []; name=:test_plate), [hx])
    simpl = structural_simplify(sys)
    prob = NonlinearProblem(simpl, [])

	solver = NewtonRaphson(; 
		linesearch = BackTracking(; maxstep=20), 
		autodiff = AutoFiniteDiff()
	)

	sol = NonlinearSolve.solve(prob, solver)

    println("=== Plate HX with Properties Test ===")
    println("Q = ", sol[hx.Q], " W")
    println("Q  ev = ", sol[hx.ev.Q], " W")
    println("Q sh = ", sol[hx.sh.Q], " W")
    println("Tin brine= ", sol[hx.T_in_1]-273.15, " C")
    println("T mid brine= ", sol[hx.T_mid_1]-273.15, " C")
    println("T out brine= ", sol[hx.T_out_1]-273.15, " C")
    println("T mid = ", sol[hx.sh.propref.Tsat]-273.15, " C")
    println("T out ev = ", sol[hx.sh.T_out_2]-273.15, " C")
    println("h in = ", sol[hx.h2i], " J/kg")
    println("hmid = ", sol[hx.h2mid], " J/kg")
    println("hmid = ", sol[hx.sh.propref.hd], " J/kg")
    println("h out = ", sol[hx.h2o], " J/kg")
    println("=== Length fraction ===")
    println("EV = ", (1-zap(sol[hx.uz_fr_sh]))*100, " %")
    println("SH = ", zap(sol[hx.uz_fr_sh])*100, " %")
    println("=== HTC ===")
    println("EV h1 = ", sol[hx.ev.h_1], " W/m2/K")
    println("EV h2 = ", sol[hx.ev.h_2], " W/m2/K")
    println("SH h1 = ", sol[hx.sh.h_1], " W/m2/K")
    println("SH h2 = ", sol[hx.sh.h_2], " W/m2/K")
    return sol
end

# ╔═╡ 6e831997-cf76-4ea4-815f-8f5ffd48d3e1
sol_evap = test_ev()

# ╔═╡ 7f6ebcea-f84c-4cae-bdb7-b9f8117f29f1
sol_evap |> get_variables

# ╔═╡ e8abe6ea-9eda-4f83-abdc-da2bd91dd26a
sol_evap.retcode

# ╔═╡ 6e3b84d1-afe1-4c56-98ba-50b2f13cfc46
sol_evap.resid

# ╔═╡ c81e49af-a294-45d2-b2fd-c605dad70ef3
sol_evap.stats

# ╔═╡ 9d816dac-a346-46d0-8394-bd87c39b6f7b
md"""Next step is to link it to new "real" prop elements for the inlet and outlet, this way we'll be able to constrain the enthalpy and delta T? Not sure it's actually necessary"""

# ╔═╡ e9738240-eb5a-49bd-8c37-46f8a8c84f7e
md"""Need to figure out a better htc correlation"""

# ╔═╡ 9b628eb9-cc22-496e-b0f9-0327e52fb0d1
"""

"""
function ThermodynamicCycle(; name,
                     refrigerant = "R134A",
                     fluid       = "MEA",
                     χ           = 0.28,
                     L           = 0.2,
                     Per         = 2π * 0.05 * 1000)
    # ——— Parameters ———
    ps = @parameters begin
    end

    # ——— State variables ———
    sts = @variables begin
		T_i_b = 9+273.154e5
		T_cond = 45
		p_cond = 2e5
    end

    # ——— Property components ———
    @named evaporator  = Evaporator()
    @named compressor = Compressor(fluid = refrigerant)
	@named inlet_exp  = Properties.RefrigerantProps(fluid = refrigerant)
	@named inlet_cond  = Properties.RefrigerantProps(fluid = refrigerant)

    # ——— Governing equations ———
	
    eqs = [
		
		#hx.uz_fr_sh ~ 0.0,
        evaporator.L      ~ 0.2,
        evaporator.γṁ_1   ~ log(0.25),
        evaporator.γṁ_2   ~ log(0.25),
        #hx.γṁ_2   ~ log(0.25),
        evaporator.T_in_1 ~ 320.0,
        evaporator.T_out_1 ~ evaporator.T_in_1 - 3,
        evaporator.h2i ~ 350.0e3,#
		#hx.p_2 ~ 2e5,
		evaporator.h2mid ~ evaporator.sh.propref.hd,
		#hx.uz_fr_sh ~ 0.2,
    ]
    eqs = [
		compressor.i.h ~ evaporator.o_ref.h,
		compressor.o.h ~ inlet_cond.h,
		
		inlet_exp.h ~ inlet_exp.hd,
		
		inlet_exp.h ~ evaporator.i_ref.h,
		
		p_l ~ evaporator.p_ref,
		p_l ~ compressor.i.p,
		
		p_h ~ compressor.o.p,
		p_h ~ inlet_cond.p,
		p_h ~ inlet_exp.p,
		p_h ~4e5,
		
       	evaporator.U_SH    ~ 1000,
        evaporator.U_EV    ~ 3000,

        exp(evaporator.γṁ_ref)  ~ compressor.ṁ,
		
    ]

    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name = name)
    return ModelingToolkit.compose(sys,
        [evaporator, compressor, inlet_exp, inlet_cond], name = name)
end

# ╔═╡ da834208-8f98-4527-924b-28ceb7587fdb
"""
Test routine for the Evaporator model
"""
function test_thermodynamic_cycle()
    
	@named cycle = ThermodynamicCycle();
    eqs = [
		
		#cycle.compressor.f ~ 25,
        # Brine inlet/outlet temps and flows
        cycle.evaporator.i_sf.T  ~ 273.15 + 8,
        cycle.evaporator.o_sf.T  ~ cycle.evaporator.mid_sf.T - 3,
        cycle.evaporator.γṁ_ref  ~ log(0.05),
        # Pseudo-angle to enforce length relation
        cycle.evaporator.αL      ~ 3.0,
    ]

    @named sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [cycle])
    simpl       = structural_simplify(sys)
    prob        = NonlinearProblem(simpl, [])
    #sol         = NonlinearSolve.solve(prob, abstol = 1e-6)
	solver = NewtonRaphson(; 
            linesearch = BackTracking(; maxstep=20), 
            autodiff = AutoFiniteDiff()
        )
    sol         = NonlinearSolve.solve(prob, solver; abstol = 1e-6)


    return sol
end


# ╔═╡ ba5d7f7e-d8ab-4852-a7e3-ad7e9c0e1329
test_thermodynamic_cycle()

# ╔═╡ 69c6ad7d-5b14-4422-b3c7-5ca71add8891
"""
Test routine for the Evaporator model
"""
function test_cycle()
    
	@named cycle = ThermodynamicCycle();
    eqs = [
		
		#cycle.compressor.f ~ 25,
        # Brine inlet/outlet temps and flows
        cycle.evaporator.i_sf.T  ~ 273.15 + 8,
        cycle.evaporator.o_sf.T  ~ cycle.evaporator.mid_sf.T - 3,
        cycle.evaporator.γṁ_ref  ~ log(0.05),
        # Pseudo-angle to enforce length relation
        cycle.evaporator.αL      ~ 3.0,
    ]

    @named sys   = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [cycle])
    simpl       = structural_simplify(sys)
    prob        = NonlinearProblem(simpl, [])
    #sol         = NonlinearSolve.solve(prob, abstol = 1e-6)
	solver = NewtonRaphson(; 
            linesearch = BackTracking(; maxstep=20), 
            autodiff = AutoFiniteDiff()
        )
    sol         = NonlinearSolve.solve(prob, solver; abstol = 1e-6)


    return sol
end


# ╔═╡ d0fea23b-5993-485a-ba84-9f601b4e96b5
"""

"""
function Cycle(; name,
                     refrigerant = "R134A",
                     fluid       = "MEA",
                     χ           = 0.28,
                     L           = 0.2,
                     Per         = 2π * 0.05 * 1000)
    # ——— Parameters ———
    ps = @parameters begin
    end

    # ——— State variables ———
    sts = @variables begin
		p_h = 4e5
		p_l = 2e5
    end


    # ——— Property components ———
    @named evaporator  = Evaporator()
    @named compressor = Compressor(fluid = refrigerant)
	@named inlet_exp  = Properties.RefrigerantProps(fluid = refrigerant)
	@named inlet_cond  = Properties.RefrigerantProps(fluid = refrigerant)

    # ——— Governing equations ———
    eqs = [
		compressor.i.h ~ evaporator.o_ref.h,
		compressor.o.h ~ inlet_cond.h,
		
		inlet_exp.h ~ inlet_exp.hd,
		
		inlet_exp.h ~ evaporator.i_ref.h,
		
		p_l ~ evaporator.p_ref,
		p_l ~ compressor.i.p,
		
		p_h ~ compressor.o.p,
		p_h ~ inlet_cond.p,
		p_h ~ inlet_exp.p,
		p_h ~4e5,
		
       	evaporator.U_SH    ~ 1000,
        evaporator.U_EV    ~ 3000,

        exp(evaporator.γṁ_ref)  ~ compressor.ṁ,
		
    ]

    # ——— Assemble system ———
    sys = NonlinearSystem(eqs, sts, ps, name = name)
    return ModelingToolkit.compose(sys,
        [evaporator, compressor, inlet_exp, inlet_cond], name = name)
end

# ╔═╡ 4437611d-700b-499e-91e0-cbc01afe3adb
md" # Combined problem"

# ╔═╡ 26da88e5-8a54-4b85-aaf3-b02993279aa4
md" # BHE + Evap + Compressor"

# ╔═╡ 8a327f7b-23e9-4af6-b58a-046a32d464a3
md" # Internal forced convection"

# ╔═╡ ef749cbd-01ce-45fc-9707-7b47417ffe58
function h_tube(;name)
	ps = @parameters begin
		cp
		μ
		λ
		ρ
		Di
	end
	sts = @variables begin
		h
		ṁ
	end
	csts = @constants begin
		Re_low = 2300
		Re_high = 5000
	end
	
	V̇ = ṁ / ρ
	S = π * Di ^ 2 / 4
	v = V̇ / S

	Re = ρ * v * Di / μ
	Pr = cp * μ / λ

	Nu_lam = 3.66
	f = 0.3164*Re^(-1//4)
	Nu_tur = (f/8*(Re-1000)*Pr)/(1+12.7*sqrt(f/8)*(Pr^(2//3)-1))
	turb_index = 0.5 * (1 + SpecialFunctions.erf((2 * SpecialFunctions.erfinv(0.98) / (Re_high - Re_low)) * (Re - (Re_low + Re_high)/2)))

	Nu = (1-turb_index) * Nu_lam + turb_index * Nu_tur

    eqs = [h ~ Nu * λ / Di]
	
	return NonlinearSystem(eqs, sts, ps; name = name);
end

# ╔═╡ c73d8d96-4ac7-4748-a844-a8f2bfbe0458
#TDL: make the htc into other components
function Evaporator1(;name, refrigerant::String = "R134A")
	ps = @parameters begin
		W = 0.2, [description = "Plate Heat eXchanger width"]
		np = 20.0, [description = "Total number of plate"]
		ϕ = 1.118, [description = "Surface enlargment factor"]
		b = 2e-3, [description = "Pressing depth "]
		β = 30, [description = "Chevron angle (°)"]
		χ = 0.38, [description = "Secondary fluid concentration"]
		L_PHX = 0.3, [description = "Plate Heat eXchanger length"]
		th = 1e-3, [description = "Plate thickness"]
		λ_Cu = 300.0, [description = "Wall (copper) thermal conductivity"]
	end
	sts = @variables begin
		ṁ_r = 1.0
		ṁ_b = 1.0
		T_p_b = 1.0
		T_p_r = 1.0
		T_evap = 1.0
		T_b_in = 1.0
		αL = 1.0
		T_b_interface = 1.0
		T_r_interface = 1.0
		T_r_out = 1.0
		T_b_out = 1.0
		h_r_interface = 1.0
		h_r_in = 1.0
		K_SH_r_ = 1.0
		K_EV_r_ = 1.0
		K_SH_b_ = 1.0

		h1 = 1.0
		h2 = 1.0
	end

	h_dew = h_dew_(T_evap)
	h_bub = h_bub_(T_evap)
	p_evap = p_PC(T_evap)

	@named gphase = Prop_REF_gas()
	@named lphase = Prop_REF_liq()
	@named brine = Prop_SF(χ0 = χ)

	#aliases
	
	cp_r_gas = gphase.cp
	μ_r_gas = gphase.μ
	λ_r_gas = gphase.λ
	ρ_r_gas = gphase.ρ
	
	cp_r_liq = lphase.cp
	μ_r_liq = lphase.μ
	λ_r_liq = lphase.λ
	ρ_r_liq = lphase.ρ
	
	cp_b = brine.cp
	μ_b = brine.μ
	λ_b = brine.λ
	ρ_b = brine.ρ
	
	dh = 2*b                       # hydraulic diameter 
	S = (b*np)/2 #assuming half volume/section in each channel
	V̇_b = abs(ṁ_b)/ρ_b
	
	v_b = V̇_b/S                  # speed of the secondary fluid 
	
	V̇_r_gas = abs(ṁ_r)/ρ_r_gas
	v_r_gas = V̇_r_gas/S

	Re_raw_b = ρ_b * v_b * dh / μ_b           # Reynolds
	Re_b = sqrt(Re_raw_b^2 + 1e-12) # < maybe solves some NaN problems?
	
	Pr_b = cp_b * μ_b / λ_b
	Nul_b = 3.65/(β^0.455 * Re_b^0.339)
	Nut_b = 12.6/(β^1.142) * Re_b^(0.642 + 0.00111*β)
	Nu_b  = (Nul_b^3 + Nut_b^3)^(1/3) 
	#very wrong correlation ⚠
	Nudh_b = Nu_b * Pr_b^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
	h_b = Nudh_b*λ_b/dh
	
	Re_r_gas = ρ_r_gas * v_r_gas * dh / μ_r_gas           # Reynolds
	
	Pr_r_gas = cp_r_gas * μ_r_gas / λ_r_gas
	Nul_r_gas = 3.65/(β^0.455 * Re_r_gas^0.339)
	Nut_r_gas = 12.6/(β^1.142) * Re_r_gas^(0.642 + 0.00111*β)
	Nu_r_gas  = (Nul_r_gas^3 + Nut_r_gas^3)^(1/3) 
	Nudh_r_gas = Nu_r_gas * Pr_r_gas^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
	h_r_gas = Nudh_r_gas*λ_r_gas/dh


	U_SH = 1/(1/h_b + th/λ_Cu + 1/h_r_gas)
	U_EV = 1/(1/h_b + th/λ_Cu)

	Per_PHX = (W*ϕ)*np #only perimeter in contact with the other channels
	
	L_EV = sin(αL)^2 * L_PHX
	L_SH = cos(αL)^2 * L_PHX
	#L_SH = L_PHX

	UA_SH = U_SH * L_SH * Per_PHX
	UA_SH = L_SH * Per_PHX
	UA_EV = U_EV * L_EV * Per_PHX
	UA_EV = L_EV * Per_PHX

	# "right" side of the HX; superheating.
	# Term regarding the superheated gas
	K_SH_r_ = 2e3 / ṁ_r / cp_r_gas
	# Term regarding the brine (single phase)
	K_SH_b_ = 2e3 / ṁ_b / cp_b

	# "left" side of the HX: evaporation.
	# Term regarding the brine (single phase)
	K_EV_b_ = UA_EV / ṁ_b / cp_b
	
	# --- SH
	
	#variables that will not be present in the end system
	@variables x_SH
	T_b_SH = T_evap
		+ (K_SH_r_/(K_SH_b_+K_SH_r_))*(T_b_in - T_evap)
		+ (K_SH_b_/(K_SH_b_+K_SH_r_))*(T_b_in - T_evap)*exp(-(K_SH_b_+K_SH_r_)*(L_SH-x_SH))
	#T_b_SH = T_evap + (T_b_in - T_evap)*exp(-K_SH_b_*x)
	T_b_interface_ = substitute(T_b_SH, x_SH => 0.0)
	#T_b_interface_ = substitute(T_b_SH, x => L_SH)

	#T_r_SH = T_evap + (α_r_*(T_b_in - T_evap)*exp(-α_SH_*L_SH))/(α_r_ + α_SH_) * (exp((α_r_ + α_SH_)*x) - 1)/exp(α_r_*x)
	#T_r_SH = T_evap + (K_SH_r_*(T_b_in - T_evap))/(K_SH_r_ + K_SH_b_)*(1 - exp(-(K_SH_r_ + K_SH_b_)*x))
	T_r_SH = T_evap + (K_SH_r_/(K_SH_b_+K_SH_r_))*(T_b_in - T_evap)*(1 - exp(-(K_SH_b_+K_SH_r_)*x_SH))

	
	# --- EV

	#variables that will not be present in the end system
	@variables x_EV
	
	T_b_EV = T_evap + (T_b_interface_ - T_evap)*exp(-K_EV_b_*(L_EV-x_EV))
	T_r_EV = T_evap
	h_r_EV = h_dew - ṁ_b*cp_b/ṁ_r*(T_b_interface_-T_b_EV)

	eqs = [
		h1 ~ h_b
		h2 ~ h_r_gas
		K_SH_r_ ~ K_SH_r_
		K_SH_b_ ~ K_SH_b_
		K_EV_r_ ~ K_EV_r_
		p_evap ~ lphase.p
		p_evap ~ gphase.p
		
		T_p_r ~ lphase.T
		T_p_r ~ gphase.T
		
		T_p_b ~ brine.T
		
		T_b_interface ~ T_b_interface_
		T_r_interface ~ substitute(T_r_SH, x_SH => 0.0)
		T_r_out ~ substitute(T_r_SH, x_SH => L_SH)
		T_b_out ~ substitute(T_b_EV, x_EV => 0.0)
		h_r_interface ~ h_dew
		h_r_in ~ substitute(h_r_EV, x_EV=> 0.0)
	]
	return ModelingToolkit.compose(NonlinearSystem(eqs, sts, ps; name = name), [gphase, lphase, brine]);
end

# ╔═╡ 2d93742a-cd89-47d9-91f6-39ef3ca294fd
begin
@named Evap1 = Evaporator1();
	
	eqs_bc = [
		Evap1.ṁ_r ~ 1.0
		Evap1.ṁ_b ~ 1.0
		
		Evap1.T_evap ~ -10+273.15
		Evap1.T_p_r ~ Evap1.T_evap + 5
		Evap1.T_p_b ~ Evap1.T_p_r + 5
		
		Evap1.T_b_in ~ -2+273.15
		Evap1.αL ~ pi/2*0.5
	]
	@named bcevap = ModelingToolkit.compose(NonlinearSystem(eqs_bc, [], []; name = :CJUBUIOname), [Evap1]);
	sbce = structural_simplify(bcevap)
	probev = NonlinearProblem(sbce, [])
	sol_ev = NonlinearSolve.solve(probev, Klement())
end;

# ╔═╡ 78f31f40-509d-472c-9ada-a738c6137dca
sol_ev.retcode

# ╔═╡ d118db4b-a0e1-43b6-85e0-34730a25cd5a
sol_ev.resid

# ╔═╡ d77a1d30-41f4-49fb-9080-e3539d562ef4
sol_ev[Evap1.T_b_interface] .- sol_ev[Evap1.T_b_in]

# ╔═╡ e4a528e3-09b6-4611-ad1c-87e034e9b5e1
sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface]

# ╔═╡ f5df3fe2-06d2-4866-9527-b54e1cfe0e9e
(sol_ev[Evap1.T_b_interface] .- sol_ev[Evap1.T_b_in])/(sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface])

# ╔═╡ acc65d45-264e-4ff7-abab-d8a5eaf473e7
((sol_ev[Evap1.T_b_in] .- sol_ev[Evap1.T_b_interface])/(sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface]))/(sol_ev[Evap1.gphase.cp] * sol_ev[Evap1.ṁ_r] / sol_ev[Evap1.brine.cp] / sol_ev[Evap1.ṁ_b])

# ╔═╡ b1b07235-fc31-4a86-ba6f-d47ae8a52427
sol_ev[Evap1.gphase.cp] / sol_ev[Evap1.brine.cp]

# ╔═╡ 617aac0d-417c-4261-84c9-98c557148bc4
md" # 1"

# ╔═╡ adedee77-c7ee-4353-bcd5-b3ab89b66e61
(sol_ev[Evap1.T_b_in] .- sol_ev[Evap1.T_b_interface]) * sol_ev[Evap1.brine.cp] * sol_ev[Evap1.ṁ_b]

# ╔═╡ 4ce72a1e-5d92-48eb-b868-4941c8135fc6
(sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface]) * sol_ev[Evap1.gphase.cp] * sol_ev[Evap1.ṁ_r]

# ╔═╡ 5a4a7a0b-337e-43ed-b094-c0492089f25d
md" # 2"

# ╔═╡ 1743b928-9e5a-466f-a5a4-260dab7e3234
(sol_ev[Evap1.T_b_interface] .- sol_ev[Evap1.T_b_out]) * sol_ev[Evap1.brine.cp] * sol_ev[Evap1.ṁ_b]

# ╔═╡ 69e6eb61-2d3f-4eb3-943f-a5eaf619a341
(sol_ev[Evap1.h_r_interface] .- sol_ev[Evap1.h_r_in]) * sol_ev[Evap1.ṁ_r]

# ╔═╡ a1216c06-5404-40f3-9c12-4cfa57356d4c
sol_ev[Evap1.T_r_interface]

# ╔═╡ c4e9bd10-4946-4227-b2ec-1083c627c843
sol_ev[Evap1.T_evap]

# ╔═╡ 0255698a-8018-4719-afc9-8f108abd1dd5
sol_ev[Evap1.gphase.cp]

# ╔═╡ c801b8e4-cbb4-426e-bc84-8b50114ae62b
sol_ev[Evap1.h1]

# ╔═╡ 5a137256-ee74-41f7-8758-f9a1e2408fb1
sol_ev[Evap1.h2]

# ╔═╡ f010d662-9cfd-4436-9760-4507371d3fba
md" ### Test props"

# ╔═╡ e563b30e-38b6-46db-9a74-184dd0a9117c
md" # Quick analysis: parametric study"

# ╔═╡ efea74fa-e2f8-4e68-944c-60ac7dabd571
begin
	m_Y = [median(y.particles) for y in YY]
	q_05 = [quantile(y.particles, 0.025) for y in YY]
	q_95 = [quantile(y.particles, 0.975) for y in YY]
end

# ╔═╡ ef2eec7f-a94f-4ce1-90fe-2190efc78a5d
begin
fig_study = Makie.lines(XX, m_Y)
	band!(XX, q_05, q_95)
	fig_study
end

# ╔═╡ f10ce6fe-b32c-44a6-860f-77a10405b764
md" # Temperature profiles"

# ╔═╡ dbf11406-e5be-441d-baed-ccd0a0e94d6b
md"# Utilities"

# ╔═╡ b6e75d0b-ba51-4513-acb1-26b61a5d66c7
md" # Other"

# ╔═╡ 4aa25799-ffec-425a-ba91-d577b7f8b2f7
candidate_solvers = [
    NewtonRaphson(),
    NonlinearSolve.TrustRegion(),
    LevenbergMarquardt(),
    PseudoTransient(),
    RobustMultiNewton(),
    FastShortcutNonlinearPolyalg(),
    Broyden(),
    Klement(),
    LimitedMemoryBroyden(),
    SimpleNewtonRaphson(),
    SimpleBroyden(),
    SimpleLimitedMemoryBroyden(),
    SimpleKlement(),
    SimpleTrustRegion(),
    SimpleDFSane(),
    SimpleHalley(),
    DynamicSS(),
    SSRootfind(),#
    NLsolveJL(),#
    CMINPACK(),
    KINSOL(),
    SIAMFANLEquationsJL(),
    #NLSolversJL(),
    #PETScSNES(),
	];

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
NonlinearSolve = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
DataFrames = "~1.7.0"
Distributions = "~0.25.120"
ModelingToolkit = "~10.2.0"
MonteCarloMeasurements = "~1.4.6"
NonlinearSolve = "~4.9.0"
PlutoDevMacros = "~0.9.0"
PlutoTest = "~0.2.2"
SpecialFunctions = "~2.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "695e9d6af990a7265d5bbcf7b5c6d5053284a26c"

[[deps.ADTypes]]
git-tree-sha1 = "e2478490447631aedba0823d4d7a80b2cc8cdb32"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.14.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

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

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

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
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "4e25216b8fea1908a0ce0f5d87368587899f75be"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "a8c0f363186263d75e97a41878d10dd842797561"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "1.6.3"

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

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
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
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
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

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

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

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "2d87d7bd165c1ca0d11923a9fabe90a9d71e88a6"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.176.0"

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
git-tree-sha1 = "76292e889472e810d40a844b714743c0ffb1c53b"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.6.0"

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
git-tree-sha1 = "f311fe66bfe4e38b2f1c8d1081f06920092b57aa"
uuid = "8d63f2c5-f18a-4cf2-ba9d-b3f60fc568c8"
version = "0.4.19"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.DispatchDoctor.extensions]
    DispatchDoctorChainRulesCoreExt = "ChainRulesCore"
    DispatchDoctorEnzymeCoreExt = "EnzymeCore"

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
git-tree-sha1 = "a7e9f13f33652c533d49868a534bfb2050d1365f"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.15"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

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

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

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
git-tree-sha1 = "df32f07f373f06260cd6af5371385b5ef85dd762"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.2"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

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

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

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

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImplicitDiscreteSolve]]
deps = ["DiffEqBase", "OrdinaryDiffEqCore", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "abb8f57a2768bc152d7acafefcdda01ee374cd11"
uuid = "3263718b-31ed-49cf-8a0f-35a466e8af96"
version = "0.1.2"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

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
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JuliaFormatter]]
deps = ["CommonMark", "Glob", "JuliaSyntax", "PrecompileTools", "TOML"]
git-tree-sha1 = "56b382cd34b1a80f63211a0b009461915915bf9e"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "2.1.2"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "c47892541d03e5dc63467f8964c9f2b415dfe718"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.46"

[[deps.JuliaSyntax]]
git-tree-sha1 = "937da4713526b96ac9a178e2035019d3b78ead4a"
uuid = "70703baa-626e-46a2-a12c-08ffd08c73b4"
version = "0.4.10"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "fb7fd516de38db80f50fe15e57d44da2836365e7"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.16.0"
weakdeps = ["FastBroadcast"]

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b94257a1a8737099ca40bc7271a8b374033473ed"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.1"

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

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "866ce84b15e54d758c11946aacd4e5df0e60b7a3"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.1"

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
git-tree-sha1 = "c0d1a91a50af6778863d320761f807f641f74935"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.17.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
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
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"

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

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

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
git-tree-sha1 = "2c61ddac64df1c7a629b8b222c144e9cbb064453"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "10.2.0"

    [deps.ModelingToolkit.extensions]
    MTKBifurcationKitExt = "BifurcationKit"
    MTKCasADiDynamicOptExt = "CasADi"
    MTKDeepDiffsExt = "DeepDiffs"
    MTKFMIExt = "FMI"
    MTKInfiniteOptExt = "InfiniteOpt"
    MTKLabelledArraysExt = "LabelledArrays"

    [deps.ModelingToolkit.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    CasADi = "c49709b8-5c63-11e9-2fb2-69db5844192f"
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"
    FMI = "14a09403-18e3-468f-ad8a-74f8dda2d9ac"
    InfiniteOpt = "20393b10-9daf-11e9-18c9-8db751c92c57"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"

[[deps.MonteCarloMeasurements]]
deps = ["Distributions", "GenericSchur", "LinearAlgebra", "MacroTools", "Random", "SLEEFPirates", "StaticArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "e9f946f013632de7ba4f399bbdc339bbb28cc4ec"
uuid = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
version = "1.4.6"

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

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "453de0fc2be3d11b9b93ca4d0fddd91196dcf1ed"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.5"

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

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

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
git-tree-sha1 = "404d71dd057759f4d590191a643113485c4a482a"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.12.0"

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
git-tree-sha1 = "9c8cd0a986518ba317af263549b48e34ac8f776d"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.5.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "e3888bdbab6e0bfadbc3164ef4595e40e7b7e954"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.6.0"
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

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "08dac9c6672a4548439048089bac293759a897fd"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.26.1"
weakdeps = ["EnzymeCore"]

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoDevMacros]]
deps = ["JuliaInterpreter", "Logging", "MacroTools", "Pkg", "TOML"]
git-tree-sha1 = "72f65885168722413c7b9a9debc504c7e7df7709"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.9.0"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

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

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

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

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "2e154f7d7e38db1af0a14ec751aba33360c3bef9"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.33.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
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

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SCCNonlinearSolve]]
deps = ["CommonSolve", "PrecompileTools", "Reexport", "SciMLBase", "SymbolicIndexingInterface"]
git-tree-sha1 = "5e67b56af8615d1b47720fdc04e34482279c4559"
uuid = "9dfe8606-65a1-4bb3-9748-cb89d1561431"
version = "1.2.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

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
git-tree-sha1 = "04bbcdc8d1f7d6f667f75fbcc68728231e21fabe"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.101.0"

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
git-tree-sha1 = "7da1216346ad79499d08d7e2a3dbf297dc80c829"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.6"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "7aaa5fe4617271b64fce0466d187f2a72edbd81a"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.5.0"

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
git-tree-sha1 = "ab958b4fec46d1f1d057bb8e2a99bfdb90744646"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.20"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
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
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "PrettyTables", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "b6a641e38efa01355aa721246dd246e10c7dcd4d"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.40"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

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
git-tree-sha1 = "3d9551301d9ecdb8c193aac2ed0a3efc303494ca"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.41.0"

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
git-tree-sha1 = "d155450e6dff2a8bc2fcb81dcb194bd98b0aeb46"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.2"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

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

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "4ab62a49f1d8d9548a1c8d1a75e5f55cf196f64e"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.71"

[[deps.WeakValueDicts]]
git-tree-sha1 = "98528c2610a5479f091d470967a25becfd83edd0"
uuid = "897b6980-f191-5a31-bcb0-bf3c4585e0c1"
version = "0.1.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

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
"""

# ╔═╡ Cell order:
# ╠═d47bf527-6960-4e14-b9e5-a422f8888e66
# ╠═910c36c9-4ec1-4b0f-b044-823a229fa536
# ╠═00bcb75f-8da5-46a8-9c20-6b94861342c4
# ╠═1499b12a-8a87-4354-b4b5-d181ffd63f1b
# ╠═734d74db-7b19-42fd-89ad-78816a939565
# ╠═fa83b49d-d4ad-470a-88f8-524ae0e3191f
# ╠═3ef17b0f-7ba5-4cef-98fe-2abf3711c907
# ╠═229b5dd1-b47b-4b67-8bc5-1694cdec228f
# ╠═43349e67-c35f-4416-a8d8-76284c25042b
# ╟─004dd778-a212-4147-a9c0-ad60f1eb539b
# ╟─18f3eb6f-d15e-4992-a071-eb335930fe38
# ╟─63cf7ccc-e17f-48d7-bbf9-2419a6ebe336
# ╟─6d536fb0-e23f-4d51-8839-8d9dc2422cbd
# ╠═96465cb7-82db-4287-bdb7-9260d39e12ef
# ╠═e31c8c25-3ca9-43cd-a865-5c1500677c24
# ╠═917d72da-5985-42b2-a079-9dd14b97a791
# ╠═270071c2-f672-4f31-974c-ac304471b517
# ╠═1ad977dd-ebd6-464d-98c8-ab2148576012
# ╠═54d7923d-fd7b-482f-b48e-c1e6641f2572
# ╠═b7ee8628-d8d7-4a93-864b-64fdd68d72b0
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
# ╠═5818a465-6685-40da-b0eb-a7e035c0170c
# ╠═d70ebc3d-7082-4043-925d-220f0ac753f7
# ╟─65e682c7-af7f-43ec-aa70-e12553026888
# ╠═62cea410-f0af-4219-a7df-92a0e574d74d
# ╠═2ed52c9e-a7ec-4793-a4c4-4a1a1611636f
# ╠═e226a6e5-d8ba-4f33-8b6b-dfde3d8d404d
# ╠═e438b966-15eb-42e4-a610-0fc426385bc7
# ╠═96966964-d666-4dc2-a554-8c359e3db528
# ╠═b41466ba-fcff-4def-957b-2b96032643aa
# ╠═c1235068-16e0-401a-b8c8-a5da5cc99a6c
# ╠═7ac463ee-57d1-4d79-a9f6-dae84f35b426
# ╠═9b117e8b-7a38-49ba-9586-f3aacaafa0f3
# ╠═d24d1392-d47c-43c8-aee2-7001e66f4540
# ╠═c9ceff67-86ff-4cff-8f4d-5c32b11fe2db
# ╠═f8d39170-abda-4664-b17e-4c71e90e462a
# ╠═98c20212-3447-481c-b0d4-326a1703997c
# ╠═cafa5fd8-7d67-4cfd-82a2-bd2f6c94b737
# ╠═f7eeabcd-23e6-4733-81aa-a6bbde82837e
# ╠═9d502c24-7e38-4c3f-96ed-0e18425d82d9
# ╠═0907ee2b-5ce8-4789-936a-ec14a4c3c059
# ╠═039f80ce-620b-4863-aa98-1f0f2ed49a7e
# ╠═9fc18265-9496-4033-aed5-2717bcfd3971
# ╠═5cf6a034-9868-448e-9a6b-0c607d49194e
# ╠═b2d62ab4-212a-4ddc-a8e2-45cb144715b3
# ╠═d9f8989e-99c5-428c-9913-f690b313e44d
# ╠═ba435e92-b5cc-4f05-8b41-4608f1176c23
# ╠═eb386e80-6a11-4b38-b55a-50bbad670474
# ╠═7f4bf9ca-4444-41b9-b71c-94aa626e9ca0
# ╠═11823168-8288-49db-88e8-319df4535273
# ╠═d01db8eb-17dc-405e-b566-63fa0e86682b
# ╠═86eeb20e-7da2-489d-84af-317f31f63ef4
# ╠═9ff7dcd3-b011-43f2-9280-440c58cf8016
# ╠═9418798c-4603-498a-ab7a-a45d509aa67f
# ╠═42bda7b4-ea1c-4717-a514-a6d97aa68f93
# ╠═1bd62dd3-e84e-4e0b-8731-890cc33286b8
# ╠═6e831997-cf76-4ea4-815f-8f5ffd48d3e1
# ╠═7f6ebcea-f84c-4cae-bdb7-b9f8117f29f1
# ╠═e8abe6ea-9eda-4f83-abdc-da2bd91dd26a
# ╠═6e3b84d1-afe1-4c56-98ba-50b2f13cfc46
# ╠═c81e49af-a294-45d2-b2fd-c605dad70ef3
# ╠═2d05d7cf-fa80-4552-a2ca-a5ceac3939f1
# ╠═9d816dac-a346-46d0-8394-bd87c39b6f7b
# ╠═e9738240-eb5a-49bd-8c37-46f8a8c84f7e
# ╠═9b628eb9-cc22-496e-b0f9-0327e52fb0d1
# ╠═ba5d7f7e-d8ab-4852-a7e3-ad7e9c0e1329
# ╠═da834208-8f98-4527-924b-28ceb7587fdb
# ╠═69c6ad7d-5b14-4422-b3c7-5ca71add8891
# ╠═d0fea23b-5993-485a-ba84-9f601b4e96b5
# ╠═4437611d-700b-499e-91e0-cbc01afe3adb
# ╠═26da88e5-8a54-4b85-aaf3-b02993279aa4
# ╠═8a327f7b-23e9-4af6-b58a-046a32d464a3
# ╠═ef749cbd-01ce-45fc-9707-7b47417ffe58
# ╠═c73d8d96-4ac7-4748-a844-a8f2bfbe0458
# ╠═2d93742a-cd89-47d9-91f6-39ef3ca294fd
# ╠═78f31f40-509d-472c-9ada-a738c6137dca
# ╠═d118db4b-a0e1-43b6-85e0-34730a25cd5a
# ╠═d77a1d30-41f4-49fb-9080-e3539d562ef4
# ╠═e4a528e3-09b6-4611-ad1c-87e034e9b5e1
# ╠═f5df3fe2-06d2-4866-9527-b54e1cfe0e9e
# ╠═acc65d45-264e-4ff7-abab-d8a5eaf473e7
# ╠═b1b07235-fc31-4a86-ba6f-d47ae8a52427
# ╠═617aac0d-417c-4261-84c9-98c557148bc4
# ╠═adedee77-c7ee-4353-bcd5-b3ab89b66e61
# ╠═4ce72a1e-5d92-48eb-b868-4941c8135fc6
# ╠═5a4a7a0b-337e-43ed-b094-c0492089f25d
# ╠═1743b928-9e5a-466f-a5a4-260dab7e3234
# ╠═69e6eb61-2d3f-4eb3-943f-a5eaf619a341
# ╠═a1216c06-5404-40f3-9c12-4cfa57356d4c
# ╠═c4e9bd10-4946-4227-b2ec-1083c627c843
# ╠═0255698a-8018-4719-afc9-8f108abd1dd5
# ╠═c801b8e4-cbb4-426e-bc84-8b50114ae62b
# ╠═5a137256-ee74-41f7-8758-f9a1e2408fb1
# ╠═f010d662-9cfd-4436-9760-4507371d3fba
# ╠═e563b30e-38b6-46db-9a74-184dd0a9117c
# ╠═efea74fa-e2f8-4e68-944c-60ac7dabd571
# ╠═ef2eec7f-a94f-4ce1-90fe-2190efc78a5d
# ╠═e9a4105f-abec-47ca-9170-ca6aef7a9c49
# ╠═4bbe463f-7e79-44f8-8abe-5f104b346e45
# ╠═5870e251-b3da-42b7-b72f-211e74c669fc
# ╠═f10ce6fe-b32c-44a6-860f-77a10405b764
# ╠═dbf11406-e5be-441d-baed-ccd0a0e94d6b
# ╠═b6e75d0b-ba51-4513-acb1-26b61a5d66c7
# ╠═4aa25799-ffec-425a-ba91-d577b7f8b2f7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
