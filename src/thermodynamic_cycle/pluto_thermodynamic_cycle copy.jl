### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ f0dae363-b54b-4f43-8f75-ef1f1e525e7d
using SpecialFunctions

# ╔═╡ 0e6790ed-f7fd-49ac-939a-72ff4f782b98
using Interpolations

# ╔═╡ 42268bb5-c95c-45ab-bdc6-125897bb0fcf
using Clapeyron

# ╔═╡ 1d44ba55-15e5-4691-86e1-bc6429e34938
using SteadyStateDiffEq

# ╔═╡ 42976c82-9438-4f83-85e9-332e63bd9ce4
using NLsolve

# ╔═╡ 31371da9-384d-4ec4-95da-260f22b40bcd
using Sundials

# ╔═╡ 879d3fb5-b766-4bde-80f9-75b676e743c0
using MINPACK

# ╔═╡ b4bd9251-e114-411d-8fde-d528e6176c0c
using SIAMFANLEquations

# ╔═╡ 4fd26428-13ca-48ad-9dae-faa408eb3899
using NLSolvers

# ╔═╡ 43d00f8a-9f71-47d1-b68d-d51f3e4b6bfc
using PETSc

# ╔═╡ bc6d9b90-734e-4d9d-be0f-104a3b99aea2
# ╠═╡ skip_as_script = true
#=╠═╡
using PlutoUI
  ╠═╡ =#

# ╔═╡ 8857ee4a-9620-4071-ae78-6e63ea0b7f0e
# ╠═╡ disabled = true
#=╠═╡
using CoolProp
  ╠═╡ =#

# ╔═╡ a9eddcc7-ac94-483e-8622-e6c9dca8351b
# ╠═╡ disabled = true
#=╠═╡
using Unitful
  ╠═╡ =#

# ╔═╡ 0405ca05-0ca9-4b6c-ac93-91e5f50cab91
# ╠═╡ skip_as_script = true
#=╠═╡
using Markdown
  ╠═╡ =#

# ╔═╡ c691da35-ff83-4dc6-af9f-2be03d95187a
# ╠═╡ skip_as_script = true
#=╠═╡
using InteractiveUtils
  ╠═╡ =#

# ╔═╡ 79bc4f64-38b6-460d-b244-b3d25de13009
# ╠═╡ disabled = true
#=╠═╡
using DataFrames
  ╠═╡ =#

# ╔═╡ 4dd5e489-a155-48c3-b7da-0d459f6c7588
# ╠═╡ skip_as_script = true
#=╠═╡
using Logging
  ╠═╡ =#

# ╔═╡ f6b9843f-bb5e-4275-bea2-05e9d594207f
# ╠═╡ skip_as_script = true
#=╠═╡
begin
    using WGLMakie
    set_theme!(theme_dark())
end
  ╠═╡ =#

# ╔═╡ 1c61fab6-f56b-460b-870f-8e3a83337312
using ModelingToolkit

# ╔═╡ fa00c042-0b93-4635-ac75-4513e61671ec
using Nemo

# ╔═╡ e5f88f64-0777-4a12-8fff-ab0c2b468a7c
# ╠═╡ skip_as_script = true
#=╠═╡
using PlutoExtras
  ╠═╡ =#

# ╔═╡ 3efb8ab9-1d4b-4408-b224-a85611388348
using PlutoDevMacros

# ╔═╡ bfdb6ceb-a0e0-4f87-b650-913a3579085e
# ╠═╡ skip_as_script = true
#=╠═╡
using Printf
  ╠═╡ =#

# ╔═╡ 69504032-712c-4d3c-8231-be9f15d683de
# ╠═╡ skip_as_script = true
#=╠═╡
using Parquet
  ╠═╡ =#

# ╔═╡ e01e0c08-459e-4636-8fef-b0f6790ee964
# ╠═╡ disabled = true
#=╠═╡
using FileIO
  ╠═╡ =#

# ╔═╡ a943d92d-2672-4f52-9b9e-bfa51f042bcf
# ╠═╡ skip_as_script = true
#=╠═╡
using DuckDB
  ╠═╡ =#

# ╔═╡ 93d5f510-e03b-4ced-9ae7-ef1837b2fb34
# ╠═╡ skip_as_script = true
#=╠═╡
using PlutoTest
  ╠═╡ =#

# ╔═╡ e0a8dace-aa48-4ed7-958e-0d2595746609
using NonlinearSolve

# ╔═╡ badf578f-6cc7-48c2-8693-370856d659dd
@fromparent begin
	using ..Properties: *#Tf, Vf, Sf, Uf, Gf, Af, ρf, ρmolf, Hf, pf, Cpf, Cvf, γf, κTf, κSf, cf, αf, μJTf, Zf, λf, ηf, Tsatf, Vbf, Vdf, RefrigerantProps
end

# ╔═╡ 45b4ca38-88af-4057-bad9-c3c69e218ee2
using FiniteDiff

# Create the solver with the desired autodiff setting

# ╔═╡ 57a059fe-3929-4cc5-9591-49fa7cdb18cc
# ╠═╡ skip_as_script = true
#=╠═╡
using MonteCarloMeasurements
  ╠═╡ =#

# ╔═╡ 086bc03f-081e-428b-b3e5-62d90a632bec
# ╠═╡ skip_as_script = true
#=╠═╡
using Distributions
  ╠═╡ =#

# ╔═╡ 679aae37-451e-4bf2-b7ef-c8ba6bb9e542
md" # Imports"

# ╔═╡ 47acf615-6bce-4b01-bf71-6b9a52316942
# ╠═╡ disabled = true
#=╠═╡
# Removed custom assemble function - using ModelingToolkit.compose directly

  ╠═╡ =#

# ╔═╡ 028775a6-0c88-4824-bf6d-50fd183cb775
test_ref_properties(200e3, 10e5)

# ╔═╡ ba2f13d7-0cb8-4556-933c-bb7b5d86db96
md" # Sub Model: Compressor"

# ╔═╡ 9538bfde-6ac2-4557-aaba-6b5e55f5b774
function Compressor(; name,
	fluid= "R134A",
	η_is = 0.7,
	V = 0.2/50,
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
		#f = 50.0
		#W = 0.0
		#V̇ = 1.0
		#ṁ = 1.0
	end
	
	@named i = RefrigerantProps(fluid = fluid)
	@named o_is = RefrigerantProps(fluid = fluid)
	@named o = RefrigerantProps(fluid = fluid)
	
	eqs = [
		#prop_ref.h ~ h_val,
		#prop_ref.p ~ p_val,
		# pulling up variables
		T_i ~ i.T,
		T_o ~ o.T,
		T_o_is ~ o_is.T,
		# Enforce the isentropic condition:
		#i.s ~ o_is.s,
		#o.p ~ o_is.p,
		
		# Efficiency definition: ratio of the isentropic enthalpy change to the actual enthalpy change.
		#η_is ~ (o_is.h - i.h) / (o.h - i.h),
	
		#ṁ ~ V̇ * i.ρ,
		
		#V̇ ~ f*V,
		
		#W ~ ṁ*(o.h-i.h),

		
		
	    ];
	
	
    # Create the parent system with connected components
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [i, o_is, o], name=name)
end

# ╔═╡ 9d89f040-e0ae-4807-876a-2e089363049e
function test_compressor()
    # Create the borehole model
	@named compressor = Compressor()

	eqs=[
		compressor.i.h ~ 200e3,
		compressor.i.p ~ 1e5,
		compressor.o.h ~ 200e3,
		compressor.o.p ~ 1e5,
		compressor.o_is.h ~ 200e3,
		compressor.o_is.p ~ 1e5,
    ]

    
    # Create the parent system with connected components
    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], [], name = :test), [compressor])

    # Simplify and solve
    simpl_sys = structural_simplify(sys)
    problem = NonlinearProblem(simpl_sys, [])
	return problem
    sol = NonlinearSolve.solve(problem)
	sol.retcode
end

# ╔═╡ 11d90544-aa10-463b-ae1d-d1fa70141e42
prob = test_compressor()

# ╔═╡ 10f4aaf8-4741-4daf-8e5c-6f8765d83cef
solver = NewtonRaphson(; autodiff=AutoFiniteDiff())

# Then solve the problem using this solver

# ╔═╡ c89c0195-7f93-4a49-993a-7efa861db428
md" ## Sub model: Evaporator"

# ╔═╡ 1b08d63a-db57-44c7-b1ec-10acbcaa7951
md""" # Solving the system"""

# ╔═╡ ebf3e61a-c624-407f-9af2-c0e03a39a902


# ╔═╡ 5e2946e4-436b-4232-93f6-5549b0a1377d
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

# ╔═╡ d4aba512-dd8c-49a5-a7b1-eccf7ca24eab
begin

		desc = """
	`A` total heat exchange surface area of the plate heat exchanger\\
	`L` plate length\\
	`W` plate width\\
	`H` heat exchanger height\\
	`np` number of plates\\
	`ϕ` surface enlargement factor \\
	`b` pressing depth\\
	`β` chevron angle\\
	"""
	@parameters ṁ_r ṁ_b
	@parameters A, W, H, np, ϕ, b, β, T_p_b, T_p_r
	@parameters T_evap  T_b_in T_p_evap χ
	#@parameters μ_b λ_b ρ_b cp_b
	#@parameters μ_r_gas λ_r_gas ρ_r_gas cp_r_gas

	h_dew = h_dew_(T_evap)
	h_bub = h_bub_(T_evap)
	P_evap = p_PC(T_evap)

	cp_r_gas = cp_r(T_p_r, P_evap)
	μ_r_gas = μ_g(T_p_r, P_evap)
	λ_r_gas = λ_g(T_p_r, P_evap)
	ρ_r_gas = 1/Clapeyron.volume(modelref, P_evap, T_p_r, phase=:gas)
	
	cp_r_liq = cp_r(T_p_r, P_evap)
	μ_r_liq = μ_l(T_p_r, P_evap)
	λ_r_liq = λ_l(T_p_r, P_evap)
	ρ_r_liq = 1/Clapeyron.volume(modelref, P_evap, T_p_r, phase=:liquid)
	
	cp_b = cp_b_(T_p_b, χ)
	μ_b = μ_b_(T_p_b, χ)
	λ_b = λ_b_(T_p_b, χ)
	ρ_b = ρ_b_(T_p_b, χ)
	
	dh = 2*b                       # hydraulic diameter 
	S = (b*np)/2 #assuming half volume/section in each channel
	V̇_b = abs(ṁ_b)/ρ_b
	
	v_b = V̇_b/S                  # speed of the secondary fluid 
	
	V̇_r_gas = abs(ṁ_r)/ρ_r_gas
	v_r_gas = V̇_r_gas/S         # hydraulic diameter 

	Re_b = ρ_b * v_b * dh / μ_b           # Reynolds
	#very wrong correlation ⚠
	Pr_b = cp_b * μ_b / λ_b
	Nul_b = 3.65/(β^0.455 * Re_b^0.339)
	Nut_b = 12.6/(β^1.142) * Re_b^(0.642 + 0.00111*β)
	Nu_b  = (Nul_b^3 + Nut_b^3)^(1/3) 
	Nudh_b = Nu_b * Pr_b^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
	h_b = Nudh_b*λ_b/dh
	
	Re_r_gas = ρ_r_gas * v_r_gas * dh / μ_r_gas           # Reynolds
	
	Pr_r_gas = cp_r_gas * μ_r_gas / λ_r_gas
	Nul_r_gas = 3.65/(β^0.455 * Re_r_gas^0.339)
	Nut_r_gas = 12.6/(β^1.142) * Re_r_gas^(0.642 + 0.00111*β)
	Nu_r_gas  = (Nul_r_gas^3 + Nut_r_gas^3)^(1/3) 
	Nudh_r_gas = Nu_r_gas * Pr_r_gas^(1/3) # * (μ/μw)^0.17 # I approximate this to 1 for simplicity
	h_r_gas = Nudh_r_gas*λ_r_gas/dh


	@parameters th, λ_Cu#, h_r_gas# h_b
	
	U_SH = 1/(1/h_b + th/λ_Cu + 1/h_r_gas)
	U_EV = 1/(1/h_b + th/λ_Cu)

	@parameters L_PHX
	Per_PHX = (W*ϕ)*np #only perimeter in contact with the other channels
	
	@variables αL x
	L_EV = sin(αL)^2 * L_PHX
	L_SH = cos(αL)^2 * L_PHX

	UA_SH = U_SH * L_SH * Per_PHX
	P_r = UA_SH / ṁ_r / cp_r_gas
	α_SH = UA_SH /ṁ_b / cp_b
	
	UA_EV = U_EV * L_EV * Per_PHX
	α_EV = UA_EV / ṁ_b / cp_b
	
	# ---
	
	T_b_SH = T_evap + (T_b_in - T_evap)*exp(-α_SH*(L_SH-x))
	T_b_interface = substitute(T_b_SH, x => 0.0)
	T_b_EV = T_evap + (T_b_interface - T_evap)*exp(-α_EV*(L_EV-x))
	T_r_SH = T_evap + (P_r*(T_b_in - T_evap)*exp(-α_SH*L_SH))/(P_r + α_SH) * (exp((P_r + α_SH)*x) - 1)/exp(P_r*x)
	T_r_out = substitute(T_r_SH, x => L_SH)
	T_b_out = substitute(T_b_EV, x => 0.0)
	T_r_EV = T_evap
	h_r_interface = h_dew
	h_r_EV = h_dew - ṁ_b*cp_b/ṁ_r*(T_b_interface-T_b_EV)
	h_r_in = substitute(h_r_EV, x=> 0.0) |> simplify;
end;

# ╔═╡ 25d99fb7-1c60-4ca8-a929-126ff6140af6
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	eq_superheating = T_r_out ~ T_evap+ΔT_SH_
	eq_inlet = h_r_in ~ h_bub_(p_PC(T_evap))
	
	eqs_tz = [eq_inlet, eq_superheating]
	#eqs_tz = [eq_superheating, eq_prop_r]
	#eqs_tz = [eq_superheating, eq_prop_r, eq_prop_b]
	unknowns = [αL, ṁ_b]#, T_p_r]
	#unknowns = [αL, T_p_r, T_p_b]
	ps = [
		b,
		W,
		H,
		np,
		ϕ,
		β,
		L_PHX,
		
		ṁ_r,
		
		T_evap,
		
		T_b_in,
		T_p_b,
		T_p_r,
		
		λ_Cu,
		th,

		χ,
	]
	@named sys = NonlinearSystem(eqs_tz, unknowns, ps)#, observed = [eq_superheating])
	#sys1_ = structural_simplify(sys)
	sys_ = complete(simplify(sys))
	
	p_ = [
		χ => χ_,
		
		λ_Cu=> 300.0,
		th => 1e-3,

		#PHX
		b => 2e-3,
		W => W_,
		np => np_,
		H => np*b*2, ### meh
		ϕ => 1.189,
		β => 30.0,

		# IV:
		ṁ_r => ṁ_r_,
		ṁ_b => ṁ_r_*4,
		αL => pi/4,

		# Parameters
		L_PHX => L_PHX_,
		T_evap => 273.15 + T_evap_,
		T_b_in => T_evap+ΔT_SH_+ΔT_brine_,
		T_p_r => T_evap,
		T_p_b => T_b_in,
	]
	prob_ = NonlinearProblem(sys_, [], p_, symbolic_u0 = true);
	sol_ = NonlinearSolve.solve(prob_)
	sol_.retcode
	
end
  ╠═╡ =#

# ╔═╡ 989f50ba-1014-42ba-a454-b91bbe076f6b
#=╠═╡
sol_.resid
  ╠═╡ =#

# ╔═╡ 08c9f46e-e683-4429-a5d6-56682232c1e7
#=╠═╡
begin

	# Create a dictionary mapping the unknown to its solution value
	sol_dict = Dict(unknowns .=> map(s->sol_[s], unknowns))  # or use sol_[αL] if defined that way
	swaps = merge(Dict(p_), sol_dict)
	
	# Now substitute into your dimensional temperature expression
	T_b_SH_dim = substitute(T_b_SH, swaps)
	T_b_SH_dim = substitute(T_b_SH_dim, swaps)
	T_b_SH_dim = substitute(T_b_SH_dim, swaps)
	T_b_SH_dim = substitute(T_b_SH_dim, swaps)
	T_b_SH_func = build_function(T_b_SH_dim, x, expression=Val{false});
	
	T_b_EV_dim = substitute(T_b_EV, swaps)
	T_b_EV_dim = substitute(T_b_EV_dim, swaps)
	T_b_EV_dim = substitute(T_b_EV_dim, swaps)
	T_b_EV_dim = substitute(T_b_EV_dim, swaps)
	T_b_EV_func = build_function(T_b_EV_dim, x, expression=Val{false});
	
	T_r_SH_dim = substitute(T_r_SH, swaps)
	T_r_SH_dim = substitute(T_r_SH_dim, swaps)
	T_r_SH_dim = substitute(T_r_SH_dim, swaps)
	T_r_SH_dim = substitute(T_r_SH_dim, swaps)
	T_r_SH_func = build_function(T_r_SH_dim, x, expression=Val{false});
	
	T_r_EV_dim = substitute(T_r_EV, swaps)
	T_r_EV_dim = substitute(T_r_EV_dim, swaps)
	T_r_EV_dim = substitute(T_r_EV_dim, swaps)
	T_r_EV_dim = substitute(T_r_EV_dim, swaps)
	T_r_EV_func = build_function(T_r_EV_dim, x, expression=Val{false});
	
	h_r_EV_dim = substitute(h_r_EV, swaps)
	h_r_EV_dim = substitute(h_r_EV_dim, swaps)
	h_r_EV_dim = substitute(h_r_EV_dim, swaps)
	h_r_EV_dim = substitute(h_r_EV_dim, swaps)
	h_r_EV_func = build_function(h_r_EV_dim, x, expression=Val{false});
end;
  ╠═╡ =#

# ╔═╡ 8903d352-a21d-43bd-9030-e316af7de22d
#=╠═╡
begin

	# Now substitute into your dimensional temperature expression
	L_EV_sol = substitute(L_EV, swaps)
	L_EV_func = build_function(L_EV_sol, expression=Val{false});
	
	L_SH_sol = substitute(L_SH, swaps)
	L_SH_func = build_function(L_SH_sol, expression=Val{false});
	
	T_b_in_sol = substitute(T_b_in, swaps)
	T_b_in_sol = substitute(T_b_in_sol, swaps)
	T_b_in_func = build_function(T_b_in_sol, expression=Val{false});
end;
  ╠═╡ =#

# ╔═╡ ad4c38f4-b78d-42c9-b868-102b1a17386d
#=╠═╡
sol_.ps[T_b_in]
  ╠═╡ =#

# ╔═╡ 352f52d6-70ec-47d7-8b46-892f6f1f6184
#=╠═╡
sin(sol_[αL])^2
  ╠═╡ =#

# ╔═╡ 91478fb6-9ffc-4bc2-94ff-415e11c00c6e
md" # Combined problem"

# ╔═╡ 66def5d1-7150-44fa-96be-0540aab0e097
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	eqs_coupled = [
		T_r_out ~ T_evap+ΔT_SH_,
		h_r_in ~ h_bub_(p_PC(T_evap)),
		T_b_out ~ Tbh_in,
		Tbh_out ~ T_b_in,
		UA_12 ~ 2.0,
		]

	unknowns_coupled = [Tbh_in, T_b_in, αL, ṁ_b, UA_12]

	ps_coupled = [
		b,
		W,
		np,
		ϕ,
		β,
		L_PHX,
		
		ṁ_r,
		#ṁ_b,
		
		T_evap,
		T_wall,
		
		#T_b_in,
		T_p_b,
		T_p_r,
		
		λ_Cu,
		th,

		χ,

		#αL,

		L_bh,

		#UA_12,
		UA_1w,
		UA_2w,

		T_p_1,
		T_p_2,
	]
	@named sys_coupled = NonlinearSystem(eqs_coupled, unknowns_coupled, ps_coupled)

	sys_coupled_ = complete(simplify(sys_coupled))
	
	p_coupled_ = [
		L_bh => 300.0,
		UA_12 => 2,
		UA_1w => 2,
		UA_2w => 2,
		
		χ => χ_,
		
		λ_Cu=> 300.0,
		th => 1e-3,

		#PHX
		b => 2e-3,
		W => W_,
		np => np_,
		ϕ => 1.189,
		β => 30.0,

		# IV:
		ṁ_r => ṁ_r_,
		ṁ_b => ṁ_r_,
		αL => pi/2,

		# Parameters
		L_PHX => L_PHX_,
		
		T_evap => 273.15 + T_evap_,
		T_p_r => T_evap,
		T_wall => T_evap+20,

		
		T_b_in => 0.4T_evap+0.6T_wall,
		Tbh_in => 0.6T_evap+0.4T_wall,

		T_p_b => T_b_in,
		T_p_1 => Tbh_in,
		T_p_2 => Tbh_in,
	]
	prob_coupled_ = NonlinearProblem(sys_coupled_, [], p_coupled_, symbolic_u0 = true);
	sol_coupled_ = NonlinearSolve.solve(prob_coupled_, reltol = 1e-5, abstol = 1e-5)	
end
  ╠═╡ =#

# ╔═╡ bdda6bf9-673d-4127-8de4-4d92576a87b8
#=╠═╡
prob_coupled_
  ╠═╡ =#

# ╔═╡ 51248ccf-dcff-416f-9500-9e0099062c55
#=╠═╡
sol_coupled_.retcode
  ╠═╡ =#

# ╔═╡ 1d1d643d-0c9f-4b89-bf31-61797a59d66f
#=╠═╡
sol_coupled_.resid
  ╠═╡ =#

# ╔═╡ 74af929f-2a7d-4634-9cf6-5532f65861d1
#=╠═╡
sol_coupled_[ṁ_b]
  ╠═╡ =#

# ╔═╡ b61e1f63-465f-4086-866f-6d1cd74fb355
#=╠═╡
sin(sol_coupled_[αL])^2
  ╠═╡ =#

# ╔═╡ 5191ed40-a917-4d63-95e0-3efa29501720
#=╠═╡
sol_coupled_[Tbh_in]-273.15
  ╠═╡ =#

# ╔═╡ be518053-81ce-44eb-9fb7-b4e2c998755f
#=╠═╡
sol_coupled_[T_b_in]-273.15

  ╠═╡ =#

# ╔═╡ 582d6677-b554-4c95-aeac-cfff776a70ad
md" # BHE + Evap + Compressor"

# ╔═╡ bdd82c91-2c36-4b5f-8f0f-e4716bcf61de
# ╠═╡ skip_as_script = true
#=╠═╡
UA_12_cd_tot |> get_variables
#	UA_1w_cd_tot = 1/R_1w_cd_tot
#	UA_2w_cd_tot = 1/R_2w_cd_tot;
#r_p,
#r_b,
#UP_12_Upipe
  ╠═╡ =#

# ╔═╡ 04cf2b12-35c2-4b4b-abe5-4bce6b2ce1cf
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	@variables D_p T_p_bh
	UP_12_Upipe = substitute(UA_12_cd_tot, [r_p => D_p/2, r_b => D_b/2, D_1 => D_p-2th, D_2 => D_p-2th, th_1 => th, th_2 => th, ṁ_1 => ṁ_b, ṁ_2 => ṁ_b, T_p_1 => T_p_bh, T_p_2 => T_p_bh])
	UP_12_Upipe |> get_variables
end
  ╠═╡ =#

# ╔═╡ c900c406-5e25-467b-8771-feedfe933685
begin
	@variables ΔT_SH V_comp f_comp s_comp_in s_comp_out_is h_comp_in h_comp_out_is h_comp_out η_comp ρ_comp_in P_comp_in T_comp_in P_comp_out T_comp_out_is T_comp_out V̇_comp W_comp
	eqs_all = [
		# Define the inlet and outlet isentropic entropies:
		s_comp_in ~ Clapeyron.entropy(modelref, P_comp_in, T_comp_in, phase="unknown"),
		s_comp_out_is ~ Clapeyron.entropy(modelref, P_comp_out, T_comp_out_is, phase="unknown"),
		
		# Define the inlet and outlet enthalpies:
		h_comp_in ~ Clapeyron.enthalpy(modelref, P_comp_in, T_comp_in, phase="unknown"),
		h_comp_out_is ~ Clapeyron.enthalpy(modelref, P_comp_out, T_comp_out_is, phase="unknown"),
		h_comp_out ~ Clapeyron.enthalpy(modelref, P_comp_out, T_comp_out, phase="unknown"),
		
		# Enforce the isentropic condition:
		s_comp_in ~ s_comp_out_is,
		
		# Efficiency definition: ratio of the isentropic enthalpy change to the actual enthalpy change.
		η_comp ~ (h_comp_out_is - h_comp_in) / (h_comp_out - h_comp_in),
	
		ρ_comp_in ~ Clapeyron.mass_density(modelref, P_comp_in, T_comp_in, phase="unknown"),
		ṁ_r ~ V̇_comp * ρ_comp_in,
		V̇_comp ~ f_comp * V_comp,
		W_comp ~ ṁ_r*(h_comp_out-h_comp_in),

		T_r_out ~ T_evap+ΔT_SH,
		h_r_in ~ h_bub_(p_PC(T_evap)),
		T_b_out ~ Tbh_in,
		Tbh_out ~ T_b_in,

		
		T_p_b ~ 0.5*T_b_in+0.5*Tbh_in,
		T_p_1 ~ 0.5*T_b_in+0.5*Tbh_in,
		T_p_2 ~ 0.5*T_b_in+0.5*Tbh_in,

		UA_12 ~ 2.0#UP_12_Upipe,
	    ];

	
	unknowns_all= [T_comp_out, T_comp_out_is, s_comp_in, s_comp_out_is, h_comp_in, h_comp_out_is, h_comp_out, W_comp, ṁ_r, ρ_comp_in,V̇_comp,
		Tbh_in, T_b_in, αL, ṁ_b,
	T_p_b, T_p_1, T_p_2,
	UA_12]

	ps_all = [T_comp_in, P_comp_in, P_comp_out, η_comp, f_comp, V_comp,
		b,
		W,
		np,
		ϕ,
		β,
		L_PHX,
		
		#ṁ_r,
		#ṁ_b,
		
		T_evap,
		T_wall,
		
		#T_b_in,
		#T_p_b,
		T_p_r,
		
		λ_Cu,
		th,

		χ,

		#αL,

		L_bh,

		UA_12,
		UA_1w,
		UA_2w,

		#T_p_1,
		#T_p_2,

		ΔT_SH,

		# BHE
		#ṁ_b
		D_p,
		D_b,
		α,
		th,
		T_p_bh,
		#T_p_2,
		Re_low,
		Re_high,
		λ_ground,
		λ_grout,
		λ_PE,
		#χ,
	]
	
	guesses_all = [
		# BH parameters
		D_p => 40e-3,
		D_b => 115e-3,
		α => 0.5, #adimensional spacing
		th => 2e-3, #some check with consistency to be done across the board
		Re_low => 2300,
		Re_high => 5000,
		λ_ground => 3.5,
		λ_grout => 2.0,
		λ_PE => 0.4,
		
		L_bh => 300.0,
		UA_1w => 2,
		UA_2w => 2,
		
		χ => χ_,

		# Evaporator
		
		λ_Cu=> 300.0,
		th => 1e-3,
		
		b => 2e-3,
		W => W_,
		np => np_,
		ϕ => 1.189,
		β => 30.0,
		L_PHX => L_PHX_,

		#Compressor
		η_comp => 0.8,
		V_comp => 0.1/50,

		
		# Operating conditions
		f_comp => 50,

		
		T_evap => 273.15 + T_evap_,
		T_p_r => T_evap,
		T_wall => T_evap+20,

		
		T_b_in => 0.4T_evap+0.6T_wall,
		Tbh_in => 0.6T_evap+0.4T_wall,

		T_p_b => T_b_in,
		T_p_bh => T_b_in,

		T_p_1 => Tbh_in,
		T_p_2 => Tbh_in,
		
		ΔT_SH => 5.0,
		#---
		
		# Initial guesses:
		ṁ_r => ṁ_r_,
		ṁ_b => ṁ_r_,
		αL => pi/2,

		
		UA_12 => 2,
		
		h_comp_in => Clapeyron.enthalpy(modelref, P_comp_in, T_comp_in, phase="unknown"),
		h_comp_out_is => Clapeyron.enthalpy(modelref, P_comp_out, T_comp_out_is, phase="unknown"),
		h_comp_out => Clapeyron.enthalpy(modelref, P_comp_out, T_comp_out_is, phase="unknown"),
		s_comp_in => s_comp_out_is,
		s_comp_out_is => Clapeyron.entropy(modelref, P_comp_out, T_comp_out_is, phase="unknown"),
		T_comp_out => T_comp_in + 10.0,
		T_comp_out_is => T_comp_in + 10.0,
		T_comp_in => 293.0,
		P_comp_in => 101325.0,
		ρ_comp_in => 1.0,
		P_comp_out => 300000.0,
		V̇_comp => V_comp*f_comp,
		W_comp => 0.0,
	
	];
	
	@named all_coupled = NonlinearSystem(eqs_all, unknowns_all, ps_all);
	all_coupled_ = complete(simplify(all_coupled));
	
	prob_all = NonlinearProblem(all_coupled_, [], guesses_all, symbolic_u0 = true)
end

# ╔═╡ 31f4a311-4b38-4868-803b-e4f263468e55
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	
	sol_all = NonlinearSolve.solve(prob_all, reltol = 1e-5, abstol = 1e-5)
	sol_all.retcode
	
end
  ╠═╡ =#

# ╔═╡ 9ad5a652-0db4-42db-be19-b20407459252
md" # Internal forced convection"

# ╔═╡ beeb113e-3476-41e5-8faa-0123bae41e67
# ╠═╡ skip_as_script = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 2da011ed-131f-48d7-baac-285b2b7d147e
# ╠═╡ skip_as_script = true
#=╠═╡
@named h1 = h_tube()
  ╠═╡ =#

# ╔═╡ 075e8fb2-89ff-4936-8c67-ca8e90159a5a
# ╠═╡ skip_as_script = true
#=╠═╡
@named h2 = h_tube()
  ╠═╡ =#

# ╔═╡ 92c361eb-bb54-4565-8ac6-d6cbad4141d9
# ╠═╡ skip_as_script = true
#=╠═╡
function Prop_SF(;name, fluid::String = "ethanol", χ0 = 0.285)
	ps = @parameters begin
		χ = χ0
	end
	sts = @variables begin
		T = 300.0
		cp = 1.0
		μ = 1.0
		λ = 1.0
		ρ = 1.0
		β = 1.0
	end

    eqs = [
        cp ~ cp_b_(T, χ),
        λ ~ λ_b_(T, χ),
        μ ~ μ_b_(T, χ),
        ρ ~ ρ_b_(T, χ),
		β ~ β_b_(T, χ)
	]

	return NonlinearSystem(eqs, sts, ps; name = name);
end
  ╠═╡ =#

# ╔═╡ 3cd0e1aa-f35e-4455-ae74-5973bbedd54d
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	
		ρ_g_(p, T) = 1/Clapeyron.volume(modelref, p, T, phase=:gas)
		ρ_l_(p, T) = 1/Clapeyron.volume(modelref, p, T, phase=:liquid)
		@register_symbolic ρ_l_(p, T)
		@register_symbolic ρ_g_(p, T)
end
  ╠═╡ =#

# ╔═╡ 569d53eb-45d2-47d0-93cb-866018abeb4f
# ╠═╡ skip_as_script = true
#=╠═╡
function Prop_REF_gas(;name, refrigerant::String = "R134A")
	ps = @parameters begin
	end
	sts = @variables begin
		p = 1e5
		T = 300.0
		cp = 1.0
		μ = 1.0
		λ = 1.0
		ρ = 1.0
	end

    eqs = [
        cp ~ cp_r(T, p),
        λ ~ λ_g(T, p),
        μ ~ μ_g(T, p),
        ρ ~ ρ_g_(p, T)
	]

	return NonlinearSystem(eqs, sts, ps; name = name);
end
  ╠═╡ =#

# ╔═╡ ad17bd58-b5c6-4560-9be4-ff2eadf594d9
# ╠═╡ skip_as_script = true
#=╠═╡
function Prop_REF_liq(;name, refrigerant::String = "R134A")
	ps = @parameters begin
	end
	sts = @variables begin
		p = 1e5
		T = 300.0
		cp = 1.0
		μ = 1.0
		λ = 1.0
		ρ = 1.0
	end
    eqs = [
        cp ~ cp_r(T, p),
        λ ~ λ_l(T, p),
        μ ~ μ_l(T, p),
        ρ ~ ρ_l_(p, T),
	]

	return NonlinearSystem(eqs, sts, ps; name = name);
end
  ╠═╡ =#

# ╔═╡ e83823bb-e069-4e5d-9d94-b22b09135a94
# ╠═╡ skip_as_script = true
#=╠═╡
#TDL: make the htc into other components
function Evaporator(;name, refrigerant::String = "R134A")
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
  ╠═╡ =#

# ╔═╡ 30ed60b4-8ef8-4a0a-995d-1a8e7cb1249a
# ╠═╡ skip_as_script = true
#=╠═╡
begin
@named Evap1 = Evaporator();
	
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
  ╠═╡ =#

# ╔═╡ e8fa9ff4-cf07-4410-9287-ecfedf2e8c32
#=╠═╡
sol_ev.retcode
  ╠═╡ =#

# ╔═╡ 557bea62-35b1-45d6-b1fc-34eecde31d89
#=╠═╡
sol_ev.resid
  ╠═╡ =#

# ╔═╡ 7d05490f-9935-4b28-bf99-d2e523767b30
#=╠═╡
sol_ev[Evap1.T_b_interface] .- sol_ev[Evap1.T_b_in]
  ╠═╡ =#

# ╔═╡ 8d3c89a4-978c-4eda-a6f5-2aab8e86b48c
#=╠═╡
sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface]
  ╠═╡ =#

# ╔═╡ 4f2877d3-f3e2-4926-b56e-bc52115c1dc1
#=╠═╡
(sol_ev[Evap1.T_b_interface] .- sol_ev[Evap1.T_b_in])/(sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface])
  ╠═╡ =#

# ╔═╡ a1cb28de-24b3-419a-9f39-8e9aac89adfa
#=╠═╡
((sol_ev[Evap1.T_b_in] .- sol_ev[Evap1.T_b_interface])/(sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface]))/(sol_ev[Evap1.gphase.cp] * sol_ev[Evap1.ṁ_r] / sol_ev[Evap1.brine.cp] / sol_ev[Evap1.ṁ_b])
  ╠═╡ =#

# ╔═╡ f4df2ef9-62ed-4f7d-808c-0f209339c166
#=╠═╡
sol_ev[Evap1.gphase.cp] / sol_ev[Evap1.brine.cp]
  ╠═╡ =#

# ╔═╡ 54e30cff-1aa0-498b-90d9-bcf3dbe4dd4e
md" # 1"


# ╔═╡ 34bc09e0-e2c9-4c31-8369-784facd0fad3
#=╠═╡
(sol_ev[Evap1.T_b_in] .- sol_ev[Evap1.T_b_interface]) * sol_ev[Evap1.brine.cp] * sol_ev[Evap1.ṁ_b]
  ╠═╡ =#

# ╔═╡ ed25cda0-4975-4909-8c37-1eef16844a4e
#=╠═╡
(sol_ev[Evap1.T_r_out] - sol_ev[Evap1.T_r_interface]) * sol_ev[Evap1.gphase.cp] * sol_ev[Evap1.ṁ_r]
  ╠═╡ =#

# ╔═╡ c96e61c1-92dd-46de-9f84-cfd87348764b
md" # 2"


# ╔═╡ fc2368c4-763c-4322-a3fe-6ed5ac018591
#=╠═╡
(sol_ev[Evap1.T_b_interface] .- sol_ev[Evap1.T_b_out]) * sol_ev[Evap1.brine.cp] * sol_ev[Evap1.ṁ_b]
  ╠═╡ =#

# ╔═╡ 6596bb7e-2282-4f84-8969-2c2cf10428ec
#=╠═╡
(sol_ev[Evap1.h_r_interface] .- sol_ev[Evap1.h_r_in]) * sol_ev[Evap1.ṁ_r]

  ╠═╡ =#

# ╔═╡ 98320e7b-4e49-44df-b3eb-7d7d4800c5a1
#=╠═╡
sol_ev[Evap1.T_r_interface]
  ╠═╡ =#

# ╔═╡ 927cc01d-e543-463a-ae44-9caf39a2394e
#=╠═╡
sol_ev[Evap1.T_evap]
  ╠═╡ =#

# ╔═╡ 9cdd134d-b858-4102-9c81-29851cd0df56
#=╠═╡
sol_ev[Evap1.gphase.cp]
  ╠═╡ =#

# ╔═╡ 1a310cc4-aa70-46c7-9081-55603deeb921
#=╠═╡
sol_ev[Evap1.h1]
  ╠═╡ =#

# ╔═╡ bd3073a0-43d4-42d6-983b-e3e822c1abfb
#=╠═╡
sol_ev[Evap1.h2]

  ╠═╡ =#

# ╔═╡ 516ea25f-a9de-4171-8af7-f7ad5dbe4737
md" ### Test props"

# ╔═╡ 3fde0742-9054-456f-ad90-262e07f13819
md" # Quick analysis: parametric study"

# ╔═╡ 9058e171-6b6f-4478-a492-40cbb7cbf406
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	m_Y = [median(y.particles) for y in YY]
	q_05 = [quantile(y.particles, 0.025) for y in YY]
	q_95 = [quantile(y.particles, 0.975) for y in YY]
end
  ╠═╡ =#

# ╔═╡ 274c189d-8325-4634-927b-f40976190bb2
# ╠═╡ skip_as_script = true
#=╠═╡
begin
fig_study = Makie.lines(XX, m_Y)
	band!(XX, q_05, q_95)
	fig_study
end
  ╠═╡ =#

# ╔═╡ b7c12215-f78a-4b4f-8c4e-46228b5687dc
md" # Temperature profiles"

# ╔═╡ a9587932-391a-46ad-84e0-9c4760e534af
md"# Utilities"

# ╔═╡ 1d57cda4-8948-4b00-be29-4acbc17dd7ae
md" # Other"

# ╔═╡ 8da5d4ea-0951-4e87-8e2b-45000a2cf95a
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

# ╔═╡ f19816e2-4bdc-4467-9504-c5eb0b0d1dd0
#=╠═╡
begin
	#https://docs.sciml.ai/NonlinearSolve/stable/solvers/nonlinear_system_solvers/#NLSolvers.jl
    logger = ConsoleLogger(stderr, Logging.Error)  # Only show errors
	valid_sols = Any[]
    # Get variable names as strings and then as symbols
    names_unk = string.(unknowns)
    var_names_sym = Symbol.(names_unk)
    
    # Global DataFrame columns: Solver and the unknowns
    col_names = [:Solver; var_names_sym...]
    results_df = DataFrame([name => Any[] for name in col_names]...)

    for alg in candidate_solvers
        solver_name = string(typeof(alg))
        try
            sol_ = Logging.with_logger(logger) do
                NonlinearSolve.solve(prob_, alg)
            end
            # Create a dictionary where keys are Symbols
            sol_dict = Dict(var_names_sym .=> sol_.u)
            # Merge with a dictionary for the solver name:
            row_dict = merge(Dict(:Solver => solver_name), sol_dict)
			#row_dict.LA = sol_["LA"]
            # You can push the dictionary directly as a row:
            push!(results_df, row_dict)
            valid_sols= [valid_sols]
        catch e
        end
    end
    results_df;
	
end

  ╠═╡ =#

# ╔═╡ 38384b61-a30c-435c-a13c-228f16e6169e
md"---"

# ╔═╡ f943556f-a2fe-4cb6-a184-7abb799faa31
# ╠═╡ disabled = true
#=╠═╡
sol = NonlinearSolve.solve(prob, solver)
  ╠═╡ =#

# ╔═╡ fbf3c417-8b31-4256-8981-84b0f0ce34b0
#=╠═╡
using LinearAlgebra
  ╠═╡ =#

# ╔═╡ f2310579-6f1d-4b7e-bb72-b99e11d92885
# ╠═╡ disabled = true
#=╠═╡
using LinearAlgebra
  ╠═╡ =#

# ╔═╡ 2341b413-76d0-4b21-8f6f-58b5a5707a65
#=╠═╡
sol = NonlinearSolve.solve(prob)
  ╠═╡ =#

# ╔═╡ 5a71fb17-3bde-4c73-9af9-464ee93a893d
#=╠═╡
sol = NonlinearSolve.solve(prob)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─679aae37-451e-4bf2-b7ef-c8ba6bb9e542
# ╠═f0dae363-b54b-4f43-8f75-ef1f1e525e7d
# ╠═0e6790ed-f7fd-49ac-939a-72ff4f782b98
# ╠═42268bb5-c95c-45ab-bdc6-125897bb0fcf
# ╠═1d44ba55-15e5-4691-86e1-bc6429e34938
# ╠═42976c82-9438-4f83-85e9-332e63bd9ce4
# ╠═31371da9-384d-4ec4-95da-260f22b40bcd
# ╠═879d3fb5-b766-4bde-80f9-75b676e743c0
# ╠═b4bd9251-e114-411d-8fde-d528e6176c0c
# ╠═4fd26428-13ca-48ad-9dae-faa408eb3899
# ╠═43d00f8a-9f71-47d1-b68d-d51f3e4b6bfc
# ╠═bc6d9b90-734e-4d9d-be0f-104a3b99aea2
# ╠═8857ee4a-9620-4071-ae78-6e63ea0b7f0e
# ╠═a9eddcc7-ac94-483e-8622-e6c9dca8351b
# ╠═0405ca05-0ca9-4b6c-ac93-91e5f50cab91
# ╠═c691da35-ff83-4dc6-af9f-2be03d95187a
# ╠═79bc4f64-38b6-460d-b244-b3d25de13009
# ╠═4dd5e489-a155-48c3-b7da-0d459f6c7588
# ╠═f6b9843f-bb5e-4275-bea2-05e9d594207f
# ╠═1c61fab6-f56b-460b-870f-8e3a83337312
# ╠═fbf3c417-8b31-4256-8981-84b0f0ce34b0
# ╠═fa00c042-0b93-4635-ac75-4513e61671ec
# ╠═e5f88f64-0777-4a12-8fff-ab0c2b468a7c
# ╠═3efb8ab9-1d4b-4408-b224-a85611388348
# ╠═bfdb6ceb-a0e0-4f87-b650-913a3579085e
# ╠═69504032-712c-4d3c-8231-be9f15d683de
# ╠═e01e0c08-459e-4636-8fef-b0f6790ee964
# ╠═a943d92d-2672-4f52-9b9e-bfa51f042bcf
# ╠═93d5f510-e03b-4ced-9ae7-ef1837b2fb34
# ╠═e0a8dace-aa48-4ed7-958e-0d2595746609
# ╠═badf578f-6cc7-48c2-8693-370856d659dd
# ╠═47acf615-6bce-4b01-bf71-6b9a52316942
# ╠═028775a6-0c88-4824-bf6d-50fd183cb775
# ╟─ba2f13d7-0cb8-4556-933c-bb7b5d86db96
# ╠═9538bfde-6ac2-4557-aaba-6b5e55f5b774
# ╠═9d89f040-e0ae-4807-876a-2e089363049e
# ╠═2341b413-76d0-4b21-8f6f-58b5a5707a65
# ╠═5a71fb17-3bde-4c73-9af9-464ee93a893d
# ╠═11d90544-aa10-463b-ae1d-d1fa70141e42
# ╠═45b4ca38-88af-4057-bad9-c3c69e218ee2
# ╠═10f4aaf8-4741-4daf-8e5c-6f8765d83cef
# ╠═f943556f-a2fe-4cb6-a184-7abb799faa31
# ╟─c89c0195-7f93-4a49-993a-7efa861db428
# ╟─1b08d63a-db57-44c7-b1ec-10acbcaa7951
# ╠═ebf3e61a-c624-407f-9af2-c0e03a39a902
# ╠═5e2946e4-436b-4232-93f6-5549b0a1377d
# ╠═d4aba512-dd8c-49a5-a7b1-eccf7ca24eab
# ╠═25d99fb7-1c60-4ca8-a929-126ff6140af6
# ╠═989f50ba-1014-42ba-a454-b91bbe076f6b
# ╠═08c9f46e-e683-4429-a5d6-56682232c1e7
# ╠═8903d352-a21d-43bd-9030-e316af7de22d
# ╠═ad4c38f4-b78d-42c9-b868-102b1a17386d
# ╠═352f52d6-70ec-47d7-8b46-892f6f1f6184
# ╠═f19816e2-4bdc-4467-9504-c5eb0b0d1dd0
# ╟─91478fb6-9ffc-4bc2-94ff-415e11c00c6e
# ╠═66def5d1-7150-44fa-96be-0540aab0e097
# ╠═bdda6bf9-673d-4127-8de4-4d92576a87b8
# ╠═51248ccf-dcff-416f-9500-9e0099062c55
# ╠═1d1d643d-0c9f-4b89-bf31-61797a59d66f
# ╠═74af929f-2a7d-4634-9cf6-5532f65861d1
# ╠═b61e1f63-465f-4086-866f-6d1cd74fb355
# ╠═5191ed40-a917-4d63-95e0-3efa29501720
# ╠═be518053-81ce-44eb-9fb7-b4e2c998755f
# ╟─582d6677-b554-4c95-aeac-cfff776a70ad
# ╠═bdd82c91-2c36-4b5f-8f0f-e4716bcf61de
# ╠═04cf2b12-35c2-4b4b-abe5-4bce6b2ce1cf
# ╠═c900c406-5e25-467b-8771-feedfe933685
# ╠═31f4a311-4b38-4868-803b-e4f263468e55
# ╟─9ad5a652-0db4-42db-be19-b20407459252
# ╠═beeb113e-3476-41e5-8faa-0123bae41e67
# ╠═2da011ed-131f-48d7-baac-285b2b7d147e
# ╠═075e8fb2-89ff-4936-8c67-ca8e90159a5a
# ╠═92c361eb-bb54-4565-8ac6-d6cbad4141d9
# ╠═569d53eb-45d2-47d0-93cb-866018abeb4f
# ╠═3cd0e1aa-f35e-4455-ae74-5973bbedd54d
# ╠═ad17bd58-b5c6-4560-9be4-ff2eadf594d9
# ╠═e83823bb-e069-4e5d-9d94-b22b09135a94
# ╠═30ed60b4-8ef8-4a0a-995d-1a8e7cb1249a
# ╠═e8fa9ff4-cf07-4410-9287-ecfedf2e8c32
# ╠═557bea62-35b1-45d6-b1fc-34eecde31d89
# ╠═7d05490f-9935-4b28-bf99-d2e523767b30
# ╠═8d3c89a4-978c-4eda-a6f5-2aab8e86b48c
# ╠═4f2877d3-f3e2-4926-b56e-bc52115c1dc1
# ╠═a1cb28de-24b3-419a-9f39-8e9aac89adfa
# ╠═f4df2ef9-62ed-4f7d-808c-0f209339c166
# ╟─54e30cff-1aa0-498b-90d9-bcf3dbe4dd4e
# ╠═34bc09e0-e2c9-4c31-8369-784facd0fad3
# ╠═ed25cda0-4975-4909-8c37-1eef16844a4e
# ╟─c96e61c1-92dd-46de-9f84-cfd87348764b
# ╠═fc2368c4-763c-4322-a3fe-6ed5ac018591
# ╠═6596bb7e-2282-4f84-8969-2c2cf10428ec
# ╠═98320e7b-4e49-44df-b3eb-7d7d4800c5a1
# ╠═927cc01d-e543-463a-ae44-9caf39a2394e
# ╠═9cdd134d-b858-4102-9c81-29851cd0df56
# ╠═1a310cc4-aa70-46c7-9081-55603deeb921
# ╠═bd3073a0-43d4-42d6-983b-e3e822c1abfb
# ╟─516ea25f-a9de-4171-8af7-f7ad5dbe4737
# ╟─3fde0742-9054-456f-ad90-262e07f13819
# ╠═9058e171-6b6f-4478-a492-40cbb7cbf406
# ╠═274c189d-8325-4634-927b-f40976190bb2
# ╠═57a059fe-3929-4cc5-9591-49fa7cdb18cc
# ╠═086bc03f-081e-428b-b3e5-62d90a632bec
# ╠═f2310579-6f1d-4b7e-bb72-b99e11d92885
# ╟─b7c12215-f78a-4b4f-8c4e-46228b5687dc
# ╟─a9587932-391a-46ad-84e0-9c4760e534af
# ╟─1d57cda4-8948-4b00-be29-4acbc17dd7ae
# ╠═8da5d4ea-0951-4e87-8e2b-45000a2cf95a
# ╟─38384b61-a30c-435c-a13c-228f16e6169e
