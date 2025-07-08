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

# ╔═╡ fa83b49d-d4ad-470a-88f8-524ae0e3191f
using ModelingToolkit

# ╔═╡ 3ef17b0f-7ba5-4cef-98fe-2abf3711c907
using NonlinearSolve

# ╔═╡ 734d74db-7b19-42fd-89ad-78816a939565
@fromparent begin
	using ..Properties
	using ..Thermal  
	using ..Hydraulics
	using ..ThermodynamicCycle
end

# ╔═╡ 00996e13-d7a2-4f9e-be84-4eec83428486
using Printf

# ╔═╡ e9a4105f-abec-47ca-9170-ca6aef7a9c49
using MonteCarloMeasurements

# ╔═╡ 4bbe463f-7e79-44f8-8abe-5f104b346e45
using Distributions

# ╔═╡ 5870e251-b3da-42b7-b72f-211e74c669fc
using LinearAlgebra

# ╔═╡ 4c6f7298-8b1a-4d72-9f87-6b3d9e5a8e12
md"""
# Coupled System Components

This module contains coupled components that compose existing thermal and hydraulic components:
- `CoupledEvaporator`: Couples plate heat exchanger thermal behavior with pressure drop calculations
- `CoupledBHE`: Couples U-tube thermal model with borehole hydraulics
- Parent components enforce variable inheritance and coordinate sub-components
"""

# ╔═╡ 1a2b3c4d-5e6f-7g8h-9i0j-1k2l3m4n5o6p
function CoupledEvaporator(; name, sec_fluid="MEA", χ=0.28, refrigerant="R134A")
    """
    Coupled evaporator that composes existing components:
    - Uses existing Evaporator from ThermodynamicCycle module 
    - Adds hydraulic pressure drop components from Hydraulics module
    - Creates parent component that enforces variable inheritance
    """
    
    # ——— Variables (parent level for inheritance) ———
    vars = @variables begin
        # Mass flow rates 
        ṁ_ref = 0.1       # Refrigerant mass flow rate [kg/s]
        ṁ_sf = 1.0        # Secondary fluid mass flow rate [kg/s]
        
        # Thermal states
        T_ref_in = 270.0   # Refrigerant inlet temperature [K]
        T_ref_out = 280.0  # Refrigerant outlet temperature [K] 
        T_sf_in = 290.0    # Secondary fluid inlet temperature [K]
        T_sf_out = 285.0   # Secondary fluid outlet temperature [K]
        p_ref = 5e5        # Refrigerant pressure [Pa]
        
        # Heat transfer
        Q_total = 5000.0   # Total heat transfer [W]
        
        # Hydraulic results
        Δp_ref = 5000.0    # Refrigerant pressure drop [Pa]
        Δp_sf = 10000.0    # Secondary fluid pressure drop [Pa]
        
        # Geometry (shared between thermal and hydraulic)
        L_hx = 0.2         # Heat exchanger length [m]
        
        # Average properties for hydraulics
        ρ_ref_avg = 10.0   # Average refrigerant density [kg/m³]  
        μ_ref_avg = 1e-5   # Average refrigerant viscosity [Pa·s]
        ρ_sf_avg = 1050.0  # Average secondary fluid density [kg/m³]
        μ_sf_avg = 2e-3    # Average secondary fluid viscosity [Pa·s]
    end
    
    # ——— Create sub-components from existing modules ———
    @named thermal = ThermodynamicCycle.Evaporator(sec_fluid=sec_fluid, refrigerant=refrigerant)
    @named hydraulic_ref = Hydraulics.PlateHeatExchangerPressureDrop()    # Refrigerant side pressure drop
    @named hydraulic_sf = Hydraulics.PlateHeatExchangerPressureDrop()     # Secondary fluid side pressure drop

    eqs = [
        # ——— Enforce variable inheritance (parent -> children) ———
        # Mass flows
        thermal.γṁ_1 ~ log(ṁ_sf),     # Secondary fluid to thermal component
        thermal.γṁ_2 ~ log(ṁ_ref),    # Refrigerant to thermal component
        
        # Temperatures  
        thermal.T_in_1 ~ T_sf_in,
        thermal.T_out_1 ~ T_sf_out,
        thermal.T_out_2 ~ T_ref_out,
        
        # Pressure and heat transfer
        thermal.p_2 ~ p_ref,
        Q_total ~ thermal.Q,
        
        # ——— Connect hydraulic components ———
        # Refrigerant side hydraulics
        hydraulic_ref.ṁ ~ ṁ_ref,
        hydraulic_ref.b ~ thermal.b,        # Share plate spacing
        hydraulic_ref.W ~ thermal.W,        # Share plate width
        hydraulic_ref.L ~ thermal.L,        # Share plate length
        hydraulic_ref.np ~ thermal.np,      # Share number of plates
        hydraulic_ref.ρ ~ ρ_ref_avg,
        hydraulic_ref.μ ~ μ_ref_avg,
        Δp_ref ~ hydraulic_ref.Δp,
        
        # Secondary fluid hydraulics
        hydraulic_sf.ṁ ~ ṁ_sf,
        hydraulic_sf.b ~ thermal.b,         # Share plate spacing
        hydraulic_sf.W ~ thermal.W,         # Share plate width
        hydraulic_sf.L ~ thermal.L,         # Share plate length
        hydraulic_sf.np ~ thermal.np,       # Share number of plates
        hydraulic_sf.ρ ~ ρ_sf_avg,
        hydraulic_sf.μ ~ μ_sf_avg,
        Δp_sf ~ hydraulic_sf.Δp,
        
        # ——— Geometry consistency ———
        L_hx ~ thermal.L,              # Share length with thermal component
    ]

    sys = NonlinearSystem(eqs, vars, []; name=name)
    return ModelingToolkit.compose(sys, [thermal, hydraulic_ref, hydraulic_sf]; name=name)
end

# ╔═╡ 2d3e4f5g-6h7i-8j9k-0l1m-2n3o4p5q6r7s
function CoupledBHE(; name, sec_fluid="MEA", χ=0.28)
    """
    Coupled BHE that composes existing components:
    - Uses existing BoreholeThermal component from Thermal module
    - Adds pressure drop components from Hydraulics module 
    - Creates parent component that enforces variable inheritance
    """
    
    # ——— Variables (parent level for inheritance) ———
    vars = @variables begin
        # Mass flow and fluid properties
        ṁ = 0.5           # Mass flow rate [kg/s]
        
        # Thermal states
        T_in = 290.0      # Inlet temperature [K]
        T_out = 285.0     # Outlet temperature [K] 
        T_ground = 283.15 # Ground temperature [K]
        
        # Heat transfer
        Q_total = 2500.0  # Total heat transfer [W]
        
        # Pressure states
        p_in = 2e5        # Inlet pressure [Pa]
        p_out = 1.8e5     # Outlet pressure [Pa]
        
        # Hydraulic results
        Δp_down = 4000.0  # Pressure drop downward pipe [Pa]
        Δp_up = 4000.0    # Pressure drop upward pipe [Pa]
        Δp_bend = 2000.0  # Pressure drop at U-bend [Pa]
        Δp_total = 10000.0 # Total pressure drop [Pa]
        
        # Geometry (shared)
        L_bh = 100.0      # Borehole length [m]
        D_pipe = 0.026    # Pipe inner diameter [m]
        
        # Average properties for hydraulics
        ρ_avg = 1050.0    # Average fluid density [kg/m³]
        μ_avg = 2e-3      # Average fluid viscosity [Pa·s]
    end
    
    # ——— Create sub-components from existing modules ———
    @named thermal = Thermal.BoreholeThermal(sec_fluid=sec_fluid, χ=χ)
    @named hydraulic = Hydraulics.BoreholeHydraulics()  # Complete borehole hydraulics

    eqs = [
        # ——— Enforce variable inheritance (parent -> children) ———
        # Flow conditions
        thermal.ṁ ~ ṁ,
        thermal.T_in ~ T_in,
        thermal.T_out ~ T_out,
        thermal.T_ground ~ T_ground,
        
        # Heat transfer
        Q_total ~ thermal.Q,
        
        # Geometry
        thermal.L_bh ~ L_bh,
        
        # ——— Connect hydraulic components ———
        # Borehole hydraulics
        hydraulic.ṁ ~ ṁ,
        hydraulic.D_pipe ~ D_pipe,
        hydraulic.L_bh ~ L_bh,
        hydraulic.ρ ~ ρ_avg,
        hydraulic.μ ~ μ_avg,
        
        # Connect pressure drops
        Δp_down ~ hydraulic.Δp_down,
        Δp_up ~ hydraulic.Δp_up,
        Δp_bend ~ hydraulic.Δp_bend,
        Δp_total ~ hydraulic.Δp_total,
        p_out ~ p_in - Δp_total,
        
        # ——— Geometry consistency ———
        D_pipe ~ thermal.D_pipe,   # Share pipe diameter
    ]

    sys = NonlinearSystem(eqs, vars, []; name=name)
    return ModelingToolkit.compose(sys, [thermal, hydraulic]; name=name)
end

# ╔═╡ 3e4f5g6h-7i8j-9k0l-1m2n-3o4p5q6r7s8t
"""
Test function for CoupledEvaporator
"""
function test_coupled_evaporator()
    @named evap = CoupledEvaporator(sec_fluid="MEA", refrigerant="R134A")

    eqs = [
        # Geometry
        evap.b ~ 0.002,
        evap.np ~ 25,
        evap.W ~ 0.10,
        evap.ϕ ~ 1.2,
        evap.β ~ π/6,
        evap.th ~ 0.5e-3,
        evap.λ_pl ~ 16,
        evap.λ_pitch ~ 4e-3,
        evap.χ ~ 0.28,
        evap.L ~ 0.2,
        
        # Flow conditions
        evap.γṁ_1 ~ log(0.365),    # Secondary fluid mass flow
        evap.γṁ_2 ~ log(0.025),    # Refrigerant mass flow
        evap.T_in_1 ~ 273.15 + 10, # Secondary fluid inlet temp
        evap.T_out_1 ~ evap.T_in_1 - 3, # Secondary fluid outlet temp
        evap.h2i ~ 220.0e3,        # Refrigerant inlet enthalpy
        evap.T_out_2 ~ evap.T_evap + 7, # Refrigerant outlet temp (superheat)
        
        # Fluid properties for pressure drop
        evap.ρ_ref_avg ~ 10.0,     # Low density for vapor
        evap.μ_ref_avg ~ 1e-5,     # Low viscosity for vapor
        evap.ρ_sec_avg ~ 1050.0,   # MEA solution density
        evap.μ_sec_avg ~ 2e-3,     # MEA solution viscosity
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], []; name=:test_coupled_evap), [evap])
    simplified = structural_simplify(sys)
    
    prob = NonlinearProblem(simplified, [])
    solver = NewtonRaphson(; 
        linesearch = BackTracking(; maxstep=20), 
        autodiff = AutoFiniteDiff()
    )
    
    sol = NonlinearSolve.solve(prob, solver, abstol = 1e-6)

    println("=== Coupled Evaporator Test ===")
    println("Total Q = ", sol[evap.Q_total], " W")
    println("Refrigerant pressure drop = ", sol[evap.Δp_ref], " Pa")
    println("Secondary fluid pressure drop = ", sol[evap.Δp_sf], " Pa")
    
    return sol
end

# ╔═╡ 4f5g6h7i-8j9k-0l1m-2n3o-4p5q6r7s8t9u
"""
Test function for CoupledBHE
"""
function test_coupled_bhe()
    @named bhe = CoupledBHE(sec_fluid="MEA")

    eqs = [
        # Geometry
        bhe.L_bh ~ 100.0,
        bhe.D_bh ~ 0.15,
        bhe.D_pipe_o ~ 0.032,
        bhe.D_pipe_i ~ 0.026,
        bhe.k_grout ~ 1.5,
        bhe.k_ground ~ 2.5,
        
        # Thermal conditions
        bhe.T_ground ~ 283.15,  # 10°C
        bhe.R_f1_b ~ 0.1,
        bhe.R_f2_b ~ 0.1,
        bhe.R_12 ~ 0.05,
        
        # Flow conditions
        bhe.γṁ ~ log(0.5),      # 0.5 kg/s mass flow
        bhe.T_in ~ 278.15,      # 5°C inlet
        bhe.T_out ~ 283.15,     # 10°C outlet (heating mode)
        
        # Fluid properties
        bhe.ρ ~ 1050.0,         # MEA solution density
        bhe.μ ~ 2e-3,           # MEA solution viscosity
        bhe.ε ~ 1e-6,           # Smooth pipes
        
        # Evaluation position
        bhe.x_eval ~ 50.0,      # Middle of borehole
    ]

    @named sys = ModelingToolkit.compose(NonlinearSystem(eqs, [], []; name=:test_coupled_bhe), [bhe])
    simplified = structural_simplify(sys)
    
    prob = NonlinearProblem(simplified, [])
    solver = NewtonRaphson(; 
        linesearch = BackTracking(; maxstep=20), 
        autodiff = AutoFiniteDiff()
    )
    
    sol = NonlinearSolve.solve(prob, solver, abstol = 1e-6)

    println("=== Coupled BHE Test ===")
    println("Heat transfer rate Q = ", sol[bhe.Q_total], " W")
    println("Total pressure drop = ", sol[bhe.Δp_total], " Pa")
    println("  - Downward pipe = ", sol[bhe.Δp_down], " Pa")
    println("  - Upward pipe = ", sol[bhe.Δp_up], " Pa") 
    println("  - U-bend = ", sol[bhe.Δp_bend], " Pa")
    
    return sol
end

# ╔═╡ cell_order = [
# 	"d47bf527-6960-4e14-b9e5-a422f8888e66",
# 	"910c36c9-4ec1-4b0f-b044-823a229fa536", 
# 	"00bcb75f-8da5-46a8-9c20-6b94861342c4",
# 	"fa83b49d-d4ad-470a-88f8-524ae0e3191f",
# 	"3ef17b0f-7ba5-4cef-98fe-2abf3711c907",
# 	"734d74db-7b19-42fd-89ad-78816a939565",
# 	"00996e13-d7a2-4f9e-be84-4eec83428486",
# 	"e9a4105f-abec-47ca-9170-ca6aef7a9c49",
# 	"4bbe463f-7e79-44f8-8abe-5f104b346e45",
# 	"4c6f7298-8b1a-4d72-9f87-6b3d9e5a8e12",
# 	"1a2b3c4d-5e6f-7g8h-9i0j-1k2l3m4n5o6p",
# 	"2d3e4f5g-6h7i-8j9k-0l1m-2n3o4p5q6r7s",
# 	"3e4f5g6h-7i8j-9k0l-1m2n-3o4p5q6r7s8t",
# 	"4f5g6h7i-8j9k-0l1m-2n3o-4p5q6r7s8t9u",
# ]
# 	"5870e251-b3da-42b7-b72f-211e74c669fc",
# 	"complete-evaporator",
# 	"complete-bhe",
# 	"test-complete-evaporator",
# 	"test-complete-bhe"
# ]

# ╔═╡ complete-evaporator =#
"""
CompleteEvaporator: Integrated evaporator with both thermal and hydraulic behavior.

Combines:
- Two-phase heat transfer (evaporation + superheating)
- Pressure drop calculations for both fluid sides
- Plate heat exchanger geometry
- Property calculations for both fluids

Arguments:
- name: Component name
- refrigerant: Refrigerant fluid name (default: "R134A")
- sec_fluid: Secondary fluid name (default: "MEA")
- χ: Secondary fluid concentration (default: 0.28)
"""
function CompleteEvaporator(; name,
        refrigerant="R134A", 
        sec_fluid="MEA", 
        χ=0.28)
    
    # ——— Parameters (mostly will become variables for sizing) ———
    ps = @parameters begin
        # Plate heat exchanger geometry
        A = 1.0           # Heat transfer area [m²]
        b = 2e-3          # Channel gap [m] 
        np = 30           # Number of plates [-]
        W = 10e-2         # Plate width [m]
        L = 0.2           # Plate length [m]
        ϕ = 1.138         # Surface enlargement factor [-]
        β = π/6           # Chevron angle [rad]
        th = 5e-4         # Plate thickness [m]
        λ_pl = 16         # Plate thermal conductivity [W/m·K]
        
        # Secondary fluid properties
        χ = χ             # Concentration [-]
        
        # Operational parameters
        roughness = 1e-6  # Surface roughness [m]
    end
    
    # ——— Variables (for optimization/sizing) ———
    sts = @variables begin
        # Mass flow rates
        ṁ_ref = 0.1       # Refrigerant mass flow rate [kg/s]
        ṁ_sf = 1.0        # Secondary fluid mass flow rate [kg/s]
        
        # Temperatures
        T_ref_in = 270.0   # Refrigerant inlet temperature [K]
        T_ref_evap = 275.0 # Refrigerant evaporation temperature [K] 
        T_ref_out = 280.0  # Refrigerant outlet temperature [K]
        T_sf_in = 290.0    # Secondary fluid inlet temperature [K]
        T_sf_out = 285.0   # Secondary fluid outlet temperature [K]
        
        # Pressures
        p_ref = 5e5       # Refrigerant pressure [Pa]
        p_sf_in = 1.5e5   # Secondary fluid inlet pressure [Pa]
        p_sf_out = 1.4e5  # Secondary fluid outlet pressure [Pa]
        
        # Heat transfer
        Q_total = 5000.0  # Total heat transfer [W]
        Q_evap = 4000.0   # Heat transfer in evaporation zone [W] 
        Q_sh = 1000.0     # Heat transfer in superheat zone [W]
        
        # Heat transfer coefficients
        h_ref_evap = 2000.0  # Refrigerant HTC in evaporation [W/m²·K]
        h_ref_sh = 1000.0    # Refrigerant HTC in superheat [W/m²·K]
        h_sf = 1500.0        # Secondary fluid HTC [W/m²·K]
        
        # Hydraulic variables
        Δp_ref = 5000.0   # Refrigerant pressure drop [Pa]
        Δp_sf = 10000.0   # Secondary fluid pressure drop [Pa]
        
        # Moving boundary
        x_evap = 0.7      # Evaporation length fraction [-]
    end
    
    # ——— Property components ———
    @named ref_in = Properties.RefrigerantProps(fluid=refrigerant)
    @named ref_evap = Properties.RefrigerantProps(fluid=refrigerant) 
    @named ref_out = Properties.RefrigerantProps(fluid=refrigerant)
    @named sf_in = Properties.IncompressibleProps(fluid=sec_fluid)
    @named sf_out = Properties.IncompressibleProps(fluid=sec_fluid)
    
    # ——— Derived hydraulic parameters ———
    dh = 2*b  # Hydraulic diameter for plate channels
    v_sf = ṁ_sf / (sf_in.ρ * b * W * np/2)  # Secondary fluid velocity
    v_ref = ṁ_ref / (ref_in.ρ * b * W * np/2)  # Refrigerant velocity
    
    Re_sf = sf_in.ρ * v_sf * dh / sf_in.μ     # Reynolds number secondary fluid
    Re_ref = ref_in.ρ * v_ref * dh / ref_in.μ # Reynolds number refrigerant
    
    # ——— Equations ———
    eqs = [
        # ——— Property connections ———
        # Refrigerant inlet
        ref_in.T ~ T_ref_in,
        ref_in.p ~ p_ref,
        
        # Refrigerant at evaporation
        ref_evap.T ~ T_ref_evap,
        ref_evap.p ~ p_ref,
        ref_evap.Q ~ 1.0,  # Saturated vapor
        
        # Refrigerant outlet (superheated)
        ref_out.T ~ T_ref_out,
        ref_out.p ~ p_ref,
        
        # Secondary fluid
        sf_in.T ~ T_sf_in,
        sf_in.χ ~ χ,
        sf_out.T ~ T_sf_out, 
        sf_out.χ ~ χ,
        
        # ——— Energy balances ———
        Q_evap ~ ṁ_ref * (ref_evap.h - ref_in.h),
        Q_sh ~ ṁ_ref * (ref_out.h - ref_evap.h),
        Q_total ~ Q_evap + Q_sh,
        Q_total ~ ṁ_sf * sf_in.cp * (T_sf_in - T_sf_out),
        
        # ——— Heat transfer correlations ———
        # Secondary fluid HTC (simplified Wanniarachchi correlation)
        h_sf ~ 0.724 * (sf_in.λ/dh) * Re_sf^0.583 * (sf_in.cp*sf_in.μ/sf_in.λ)^0.4,
        
        # Refrigerant HTC in superheat (simplified)
        h_ref_sh ~ 0.023 * (ref_out.λ/dh) * Re_ref^0.8 * (ref_out.cp*ref_out.μ/ref_out.λ)^0.4,
        
        # ——— Pressure drop calculations (using Hydraulics module functions) ———
        # Secondary fluid pressure drop (simplified plate correlation)
        Δp_sf ~ (sf_in.ρ * v_sf^2 / 2) * (4 * 0.045 * Re_sf^(-0.2) * L / dh + 2.5),
        
        # Refrigerant pressure drop (simplified)
        Δp_ref ~ (ref_in.ρ * v_ref^2 / 2) * (4 * 0.045 * Re_ref^(-0.2) * L / dh + 2.5),
        
        # Pressure boundary conditions
        p_sf_out ~ p_sf_in - Δp_sf,
        
        # ——— Geometry constraints ———
        A ~ L * W * np * ϕ,  # Total heat transfer area
    ]
    
    # Create the system
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [ref_in, ref_evap, ref_out, sf_in, sf_out], name=name)
end

# ╔═╡ complete-bhe =#
"""
CompleteBHE: Complete Borehole Heat Exchanger with thermal and hydraulic behavior.

Combines:
- U-tube thermal model with ground interaction
- Pressure drop in pipes and bends
- Natural convection in annulus
- Pipe wall thermal resistance

Arguments:
- name: Component name
- fluid: Working fluid name (default: "MEA")
- χ: Fluid concentration (default: 0.28)
- depth: Borehole depth [m] (default: 100.0)
- D_bore: Borehole diameter [m] (default: 0.15)
"""
function CompleteBHE(; name,
        fluid="MEA",
        χ=0.28,
        depth=100.0,
        D_bore=0.15)
    
    # ——— Parameters ———
    ps = @parameters begin
        # Borehole geometry
        depth = depth         # Borehole depth [m]
        D_bore = D_bore       # Borehole diameter [m]
        
        # Pipe geometry
        D_pipe_o = 0.032      # Pipe outer diameter [m]
        D_pipe_i = 0.026      # Pipe inner diameter [m]
        pipe_spacing = 0.05   # Distance between pipe centers [m]
        λ_pipe = 0.4          # Pipe thermal conductivity [W/m·K]
        
        # Ground properties
        T_ground = 283.15     # Undisturbed ground temperature [K]
        λ_ground = 2.0        # Ground thermal conductivity [W/m·K]
        
        # Grout properties  
        λ_grout = 1.5         # Grout thermal conductivity [W/m·K]
        
        # Fluid properties
        χ = χ                 # Fluid concentration [-]
        
        # Surface roughness
        roughness = 1e-6      # Pipe internal roughness [m]
    end
    
    # ——— Variables ———
    sts = @variables begin
        # Flow conditions
        ṁ = 0.5               # Mass flow rate [kg/s]
        
        # Temperatures
        T_in = 290.0          # Inlet temperature [K]
        T_out = 285.0         # Outlet temperature [K]
        T_wall = 284.0        # Borehole wall temperature [K]
        
        # Pressures
        p_in = 2e5            # Inlet pressure [Pa]
        p_out = 1.8e5         # Outlet pressure [Pa]
        
        # Heat transfer
        Q_total = 2500.0      # Total heat extraction [W]
        
        # Thermal resistances
        R_pipe = 0.05         # Pipe wall resistance [K·m/W]
        R_grout = 0.15        # Grout resistance to borehole wall [K·m/W]
        R_fluid = 0.03        # Fluid convection resistance [K·m/W]
        
        # Hydraulic variables
        Δp_straight = 8000.0  # Pressure drop in straight sections [Pa]
        Δp_bend = 2000.0      # Pressure drop in U-bend [Pa]
        Δp_total = 10000.0    # Total pressure drop [Pa]
        
        # Heat transfer coefficients
        h_fluid = 2000.0      # Fluid-to-pipe HTC [W/m²·K]
        
        # Velocities and Reynolds numbers
        v_fluid = 1.0         # Fluid velocity [m/s]
        Re = 10000.0          # Reynolds number [-]
        Nu = 50.0             # Nusselt number [-]
    end
    
    # ——— Property components ———
    @named fluid_in = Properties.IncompressibleProps(fluid=fluid)
    @named fluid_out = Properties.IncompressibleProps(fluid=fluid)
    @named fluid_avg = Properties.IncompressibleProps(fluid=fluid)
    
    # ——— Thermal components from base modules ———
    @named thermal_params = Thermal.BoreholeTemperatureParameters(
        L_bh=depth, T_wall=T_wall,
        R_f1_b=R_grout, R_f2_b=R_grout, R_f1_f2=R_grout*2
    )
    
    # ——— Derived parameters ———
    A_pipe_i = π * D_pipe_i^2 / 4           # Internal pipe area
    A_pipe_o = π * D_pipe_o^2 / 4           # External pipe area
    P_pipe_i = π * D_pipe_i                 # Internal pipe perimeter
    
    # ——— Equations ———
    eqs = [
        # ——— Property connections ———
        fluid_in.T ~ T_in,
        fluid_in.χ ~ χ,
        fluid_out.T ~ T_out,
        fluid_out.χ ~ χ,
        fluid_avg.T ~ (T_in + T_out)/2,
        fluid_avg.χ ~ χ,
        
        # ——— Energy balance ———
        Q_total ~ ṁ * fluid_avg.cp * (T_in - T_out),
        
        # ——— Thermal model ———
        # Pipe wall thermal resistance
        R_pipe ~ log(D_pipe_o/D_pipe_i) / (2*π*λ_pipe),
        
        # Fluid convection resistance  
        R_fluid ~ 1 / (h_fluid * P_pipe_i),
        
        # Heat transfer coefficient correlation (Dittus-Boelter for turbulent)
        Nu ~ 0.023 * Re^0.8 * (fluid_avg.cp*fluid_avg.μ/fluid_avg.λ)^0.4,
        h_fluid ~ Nu * fluid_avg.λ / D_pipe_i,
        
        # Overall thermal resistance and temperature difference
        (T_in + T_out)/2 - T_wall ~ Q_total * (R_fluid + R_pipe + R_grout) / depth,
        
        # ——— Hydraulic model ———
        # Velocity and Reynolds number
        v_fluid ~ ṁ / (fluid_avg.ρ * A_pipe_i),
        Re ~ fluid_avg.ρ * v_fluid * D_pipe_i / fluid_avg.μ,
        
        # Pressure drop in straight sections (downward + upward)
        Δp_straight ~ 2 * (fluid_avg.ρ * v_fluid^2 / 2) * 
                      (4 * 0.0791 * Re^(-0.25) * depth / D_pipe_i),  # Simplified friction factor
        
        # Pressure drop in U-bend (180° bend with expansion/contraction)
        Δp_bend ~ (fluid_avg.ρ * v_fluid^2 / 2) * 1.5,  # Simplified bend loss coefficient
        
        # Total pressure drop
        Δp_total ~ Δp_straight + Δp_bend,
        p_out ~ p_in - Δp_total,
        
        # ——— Thermal boundary condition ———
        T_wall ~ T_ground,  # Simplified: wall temperature equals ground temperature
        
        # ——— Connect to thermal parameters ———
        thermal_params.L_bh ~ depth,
        thermal_params.T_wall ~ T_wall,
        thermal_params.T_1_in ~ T_in,
        thermal_params.ṁ ~ ṁ,
        thermal_params.cp_1 ~ fluid_avg.cp,
        thermal_params.cp_2 ~ fluid_avg.cp,
    ]
    
    # Create the system
    sys = NonlinearSystem(eqs, sts, ps, name=name)
    return ModelingToolkit.compose(sys, [fluid_in, fluid_out, fluid_avg, thermal_params], name=name)
end

# ╔═╡ test-complete-evaporator =#
"""
Test function for CompleteEvaporator component
"""
function test_complete_evaporator()
    @named evap = CompleteEvaporator(refrigerant="R134A", sec_fluid="MEA", χ=0.28)
    
    # Test conditions
    test_eqs = [
        evap.T_ref_evapor ~ 275.15,  # 2°C evaporation
        evap.T_sf_in ~ 283.15,     # 10°C secondary fluid inlet
        evap.ṁ_ref ~ 0.1,          # 0.1 kg/s refrigerant
        evap.ṁ_sf ~ 1.0,           # 1.0 kg/s secondary fluid
        evap.p_ref ~ 3.5e5,        # 3.5 bar refrigerant pressure
    ]
    
    sys = NonlinearSystem(test_eqs, [], [], name=:test_evap)
    complete_sys = ModelingToolkit.compose(sys, [evap], name=:test_evap)
    
    # Simplify and check
    simple_sys = structural_simplify(complete_sys)
    
    println("CompleteEvaporator test system created successfully")
    println("Variables: $(length(unknowns(simple_sys)))")
    println("Equations: $(length(equations(simple_sys)))")
    
    return simple_sys
end

# ╔═╡ test-complete-bhe =#
"""
Test function for CompleteBHE component
"""
function test_complete_bhe()
    @named bhe = CompleteBHE(fluid="MEA", χ=0.28, depth=100.0, D_bore=0.15)
    
    # Test conditions
    test_eqs = [
        bhe.T_in ~ 288.15,         # 15°C inlet
        bhe.ṁ ~ 0.5,               # 0.5 kg/s mass flow
        bhe.T_ground ~ 283.15,     # 10°C ground temperature
        bhe.p_in ~ 2e5,            # 2 bar inlet pressure
    ]
    
    sys = NonlinearSystem(test_eqs, [], [], name=:test_bhe)
    complete_sys = ModelingToolkit.compose(sys, [bhe], name=:test_bhe)
    
    # Simplify and check
    simple_sys = structural_simplify(complete_sys)
    
    println("CompleteBHE test system created successfully")
    println("Variables: $(length(unknowns(simple_sys)))")
    println("Equations: $(length(equations(simple_sys)))")
    
    return simple_sys
end

# ╔═╡ test_pcpwp_enhanced =#
"""
Enhanced test function for PhaseChangePlateHX_WithProperties:
prints geometry, flows, HTCs, and ALL refrigerant properties used in h2_raw calculation.
"""
function test_pcpwp_enhanced()
    @named hx = PhaseChangePlateHX_WithProperties(
        Re_low  = 2000.0,
        Re_high = 2100.0
    )

    vars = @variables begin
        # No additional variables needed
    end
    
    # Drive the model with same conditions
    eqs = [
        hx.b ~ 0.0015,
        hx.np ~ 25,
        hx.W ~ 0.10,
        hx.ϕ ~ 1.2,
        hx.β ~ π/6,
        hx.th ~ 0.5e-3,
        hx.λ_pl ~ 16,
        hx.λ_pitch ~ 4e-3,
        hx.χ ~ 0.28,
        hx.L ~ 0.2,              # plate length [m]
        
        hx.γṁ_1 ~ log(0.9),         # brine flow [kg/s]
        hx.T_in_1 ~ 273.15 + 10.0,  # brine inlet [K]
        hx.h2i ~ 250e3,             # refrigerant inlet enthalpy
        hx.h2o ~ hx.h2i + 200e3,    # refrigerant outlet enthalpy
        hx.propref.Tsat ~ hx.T_in_1 - 10,  # refrigerant evaporation temperature
    ]

    @named sys = ModelingToolkit.compose(
        NonlinearSystem(eqs, vars, []; name=:test_plate_enhanced),
        [hx]
    )
    simplified = structural_simplify(sys)

    prob = NonlinearProblem(simplified, [])
    solver = NewtonRaphson(; 
        linesearch = BackTracking(; maxstep=20), 
        autodiff = AutoFiniteDiff()
    )
    sol = NonlinearSolve.solve(prob, solver, abstol = 1e-6)
    
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
    println("  p_evap         = ", round(pev, digits=2), " bar")
    println()
    
    println("TEMPERATURES:")
    println("  T_brine_in     = ", round(T_in_1, digits=2), " °C")
    println("  T_brine_out    = ", round(T_out_1, digits=2), " °C")
    println("  T_evaporation  = ", round(T_sat, digits=2), " °C")
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
