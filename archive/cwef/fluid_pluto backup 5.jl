### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 14abfdab-a03a-4649-8aa6-14e8a4cea7bc
using Interpolations

# ╔═╡ fedaaebd-cf14-48fa-a2c9-a528f0c65033
using CoolProp

# ╔═╡ f177d2ef-ff07-444e-b353-e2484e7e8b00
using Parquet

# ╔═╡ 122c8036-63f1-4deb-8e1d-4dcf471bd298
using ModelingToolkit

# ╔═╡ 924b8d55-a8f9-4f1a-b1e0-31601ff811c9
using FileIO

# Optional DuckDB support

# ╔═╡ ee9f209e-8f89-487d-9755-67d06b42308a
begin
const USE_DUCKDB = Ref(false)
	try
    using DuckDB
    global USE_DUCKDB = Ref(true)
catch
    @warn "DuckDB not available. Long-format queries will use DataFrames filtering instead."
end
end

# ╔═╡ 0335e609-cfa7-4961-b85b-7214aeecf92e
using Makie

# ╔═╡ 5fb03f74-ad6a-4d0b-9c6e-96c9d96132ec
export FluidInterpolants, load_fluid, regenerate_fluid_data

# ╔═╡ 2dc9195f-aebc-4407-9802-35cc85cfef85
export ρ_b_, λ_b_, cp_b_, μ_b_, β_b_

# ╔═╡ e20e2ca5-9ad4-4eda-8242-229496d0c244
export list_available_fluids, query_fluid_property

# Global variable to store loaded fluid interpolants

# ╔═╡ 7af52b81-cdaa-4947-b684-3cd0b658b0a4
const LOADED_FLUIDS = Dict{String, Any}()
# DuckDB connection for fluid data querying

# ╔═╡ 7fdb2dc7-af0b-4035-bd9d-15af34f27fd1
const DB_CONN = Ref{Any}(nothing)

# Pre-declare the symbolic versions of the fluid property functions
# These will be properly defined once data is loaded

# ╔═╡ 8c986457-9d64-4be5-b5b0-d39b9702e05c
"""
    ρ_b_(T, χ)

Calculate the density of the fluid at temperature T and concentration χ.
"""
function ρ_b_ end

# ╔═╡ 3a793d39-1f74-43f6-950a-b46dcb8a3de8
"""
    λ_b_(T, χ)

Calculate the thermal conductivity of the fluid at temperature T and concentration χ.
"""
function λ_b_ end

# ╔═╡ 806ebe0f-7be4-4f35-8827-eb111d3a8891
"""
    cp_b_(T, χ)

Calculate the specific heat of the fluid at temperature T and concentration χ.
"""
function cp_b_ end

# ╔═╡ d370d15e-07a4-49a0-a6b2-b18bb760d22f
"""
    μ_b_(T, χ)

Calculate the viscosity of the fluid at temperature T and concentration χ.
"""
function μ_b_ end

# ╔═╡ 64dcc6ea-d858-4c29-ad36-534afdf3c8d8
"""
    β_b_(T, χ)

Calculate the thermal expansion coefficient of the fluid at temperature T and concentration χ.
"""
function β_b_ end

# ╔═╡ 1a4a329b-9002-471e-9013-653d14044bf1
"""
    β_b_(T)

Calculate the thermal expansion coefficient with default concentration (0.0).
"""
β_b_(T) = β_b_(T, 0.0)

# Register the symbolic functions

# ╔═╡ b49b1c08-a158-4721-bb2e-d2d2fc90201d
ModelingToolkit.@register_symbolic ρ_b_(T, χ)

# ╔═╡ 2009ccd1-5cec-4a3e-8865-876e0a4e50f8
ModelingToolkit.@register_symbolic λ_b_(T, χ)

# ╔═╡ ea8f6eb0-49c2-41a2-b5da-56b59ec858bd
ModelingToolkit.@register_symbolic cp_b_(T, χ)

# ╔═╡ e9540434-f7c5-486a-af94-0eda47ddae0d
ModelingToolkit.@register_symbolic μ_b_(T, χ)

# ╔═╡ 1d1cdfa3-a995-4286-bcf7-8f1f741810df
ModelingToolkit.@register_symbolic β_b_(T, χ)

# ╔═╡ 69f773f7-3602-4215-9caa-019406b19268
ModelingToolkit.@register_symbolic β_b_(T)

# ╔═╡ 739d9847-5a6a-4f2e-8143-167269d44be4
"""
    generate_all_fluid_data(output_dir=nothing; n_chi=100, n_temp=50)

Generate and save data for all supported fluids.
This is intended for package development to create or update all fluid parquet files at once.

# Arguments
- `output_dir`: Directory to save the fluid data (default: package data directory)
- `n_chi`: Number of concentration points to sample (default: 100)
- `n_temp`: Number of temperature points to sample (default: 50)

# Returns
- Dictionary mapping fluid names to their data directories
"""
function generate_all_fluid_data(output_dir=nothing; n_chi=100, n_temp=50)
    if isnothing(output_dir)
        output_dir = joinpath(@__DIR__, "..", "data", "fluids")
    end
    
    # List of currently supported fluids
    supported_fluids = ["water", "mea"]
    
    # Ensure directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
        @info "Created fluid data directory: $output_dir"
    end
    
    results = Dict{String, String}()
    
    # Generate data for each fluid
    for fluid_name in supported_fluids
        @info "Generating data for $fluid_name..."
        try
            # Generate fluid data
            fluid_data = generate_fluid_data(fluid_name, n_chi=n_chi, n_temp=n_temp)
            
            # Save to disk
            fluid_dir = save_fluid_data(fluid_data, output_dir)
            
            # Build interpolants to verify the data
            interpolants = build_interpolants(fluid_data)
            
            # Store the fluid directory in results
            results[fluid_name] = fluid_dir
            
            @info "Successfully generated data for $fluid_name at $fluid_dir"
        catch e
            @error "Failed to generate data for $fluid_name: $e"
        end
    end
    
    @info "Completed generating fluid data for $(length(results)) out of $(length(supported_fluids)) fluids"
    return results
end


# ╔═╡ 44ab047d-d1c0-4165-9525-b2c80572d796
test = load_fluid("mea", data_dir=joinpath(@__DIR__(), "data", "fluids"))

# ╔═╡ b9689240-5e8b-4e9e-8b9c-fef6c20ddfe8
test.itp_Tfreeze

# ╔═╡ bc04af84-5e7c-46c1-ab08-0f234d17e440
lines(test.itp_Tfreeze, 0, 1.0)

# ╔═╡ 20433cc4-97b7-44d0-958a-aff03fdd12ba
pwd()

# ╔═╡ 1022217d-2e51-4a88-ae93-f7bf47f0959d
 USE_DUCKDB[]

# ╔═╡ 72aef11f-3a07-43ae-821b-87f48dec2874
init_database()

# ╔═╡ 4d01cde3-4840-4e9c-ac3a-8c36c8a14f7d
data_dir

# ╔═╡ 432e6cac-14f9-4d09-9864-6500b306d5d0
generate_fluid_data("mea")

# ╔═╡ 426d5e79-f04c-4518-8faf-9832bb896d4f
"""
    FluidInterpolants

Stores the interpolation functions for a fluid's properties.
"""
struct FluidInterpolants
    name::String
    chi_valid::Vector{Float64}
    itp_Tfreeze::Interpolations.Extrapolation
    itp_Tmax::Interpolations.Extrapolation
    itp_2D_rho::Interpolations.Extrapolation
    itp_2D_lambda::Interpolations.Extrapolation
    itp_2D_cp::Interpolations.Extrapolation
    itp_2D_mu::Interpolations.Extrapolation
end

# ╔═╡ 7571b2d2-489b-4c5e-b3c8-6f352fcc3830
"""
    init_database(db_path::String="")

Initialize the DuckDB connection and set up necessary tables.
If db_path is empty, uses an in-memory database.
Returns true if successful, false otherwise.
"""
function init_database(db_path::String="")
    if !USE_DUCKDB[]
        return false
    end
    
    try
        if isnothing(DB_CONN[])
            # Use in-memory DB if no path provided
            conn_str = isempty(db_path) ? ":memory:" : db_path
            DB_CONN[] = DuckDB.DB(conn_str)
        end
        return true
    catch e
        @warn "Failed to initialize DuckDB: $e"
        return false
    end
end

# ╔═╡ d16ef005-fa61-4456-862f-0ec0cb5287a0
"""
    list_available_fluids(data_dir=nothing)

List all fluid datasets available in the data directory.
"""
function list_available_fluids(data_dir=nothing)
    if isnothing(data_dir)
        data_dir = joinpath(@__DIR__, "..", "data", "fluids")
    end
    
    if !isdir(data_dir)
        @warn "Data directory does not exist: $data_dir"
        return String[]
    end
    
    # Find subdirectories in the data directory - each represents a fluid
    fluids = filter(f -> isdir(joinpath(data_dir, f)), readdir(data_dir))
    
    if isempty(fluids)
        @info "No fluid datasets found in $data_dir. Use regenerate_fluid_data() to create some."
    end
    
    return fluids
end

# ╔═╡ 4eaac3ee-290b-4603-a30b-cefd15d528bd
"""
    default_fluid_approximations(fluid_name::String)

Provide simple approximate models for common fluids when data is unavailable.
These are less accurate but allow basic functionality without data files.
"""
function default_fluid_approximations(fluid_name::String)
    if lowercase(fluid_name) == "water"
        # Very simplified water properties (at atmospheric pressure)
        # Just for basic functionality when data is unavailable
        return Dict(
            "name" => "water",
            "approximation" => true,
            "ρ" => (T, χ) -> 1000.0 - 0.07 * (T - 273.15)^1.7, # kg/m³
            "λ" => (T, χ) -> 0.6, # W/(m·K)
            "cp" => (T, χ) -> 4200.0, # J/(kg·K)
            "μ" => (T, χ) -> 0.001 * exp(-0.021 * (T - 293.15)), # Pa·s
            "β" => (T, χ) -> 2.14e-4, # 1/K
            "T_freeze" => (χ) -> 273.15,
            "T_max" => (χ) -> 373.15
        )
    elseif lowercase(fluid_name) == "mea"
        # Very simplified MEA properties 
        return Dict(
            "name" => "mea",
            "approximation" => true,
            "ρ" => (T, χ) -> 1000.0 + 95.0 * χ - 0.5 * (T - 273.15), # kg/m³
            "λ" => (T, χ) -> 0.5 - 0.1 * χ, # W/(m·K)
            "cp" => (T, χ) -> 4200.0 - 950.0 * χ, # J/(kg·K)
            "μ" => (T, χ) -> 0.001 * (1.0 + 0.5 * χ) * exp(-0.02 * (T - 293.15)), # Pa·s
            "β" => (T, χ) -> 2.0e-4 * (1.0 - 0.2 * χ), # 1/K
            "T_freeze" => (χ) -> 273.15 - 12.0 * χ,
            "T_max" => (χ) -> 373.15 - 10.0 * χ
        )
    else
        return nothing
    end
end

# ╔═╡ 7f3d7705-8e7f-46dc-9c52-9fc6978f38d5
"""
    generate_fluid_data(fluid_name::String; n_chi=100, n_temp=50)

Generate raw data for interpolation for the specified fluid.
Returns a Dict with necessary data for building interpolations.
"""
function generate_fluid_data(fluid_name::String; n_chi=100, n_temp=50)
    # Determine the CoolProp fluid string based on name
    fluid_string = if lowercase(fluid_name) == "mea"
        "INCOMP::MEA"
    elseif lowercase(fluid_name) == "water"
        "INCOMP::Water"
    else
        error("Unsupported fluid: $fluid_name. Currently supported: mea, water")
    end
    
    # Sample over the χ dimension
    chi_vec = range(0.0, stop=1.0, length=n_chi)
    
    # Prepare arrays for data
    T_freeze_vals = Float64[]
    T_max_vals = Float64[]
    chi_valid = Float64[]
    
    # Collect freezing and max temperatures
    for χ in chi_vec
        try
            fluid_str = "$fluid_string[$χ]"
            T_f = CoolProp.PropsSI("T_freeze", fluid_str)
            T_max = CoolProp.PropsSI("T_max", fluid_str)
            if !isnan(T_f) && !isnan(T_max)
                push!(chi_valid, χ)
                push!(T_freeze_vals, T_f)
                push!(T_max_vals, T_max)
            end
        catch e
            # Skip invalid concentrations
        end
    end
    
    # Prepare for property interpolations
    T_tilde = range(0.0, stop=1.0, length=n_temp)
    
    # Create DataFrames for temperature limits
    df_freezing = DataFrame(chi = chi_valid, T_freeze = T_freeze_vals)
    df_max_temps = DataFrame(chi = chi_valid, T_max = T_max_vals)
    
    # Initialize property data in long format
    property_data = Vector{NamedTuple{(:chi, :T_tilde, :property, :value), 
                                      Tuple{Float64, Float64, String, Float64}}}()
    
    # Collect property data
    for (i, χ) in enumerate(chi_valid)
        T_f = T_freeze_vals[i]
        T_M = T_max_vals[i]
        fluid_str = "$fluid_string[$χ]"
        
        for (j, t̃) in enumerate(T_tilde)
            T_val = t̃ * (T_M - T_f) + T_f
            
            try
                # Density
                ρ = CoolProp.PropsSI("D", "P", 1e5, "T", T_val, fluid_str)
                if !isnan(ρ)
                    push!(property_data, (chi=χ, T_tilde=t̃, property="density", value=ρ))
                end
                
                # Thermal conductivity
                λ = CoolProp.PropsSI("L", "P", 1e5, "T", T_val, fluid_str)
                if !isnan(λ)
                    push!(property_data, (chi=χ, T_tilde=t̃, property="thermal_conductivity", value=λ))
                end
                
                # Specific heat
                cp = CoolProp.PropsSI("C", "P", 1e5, "T", T_val, fluid_str)
                if !isnan(cp)
                    push!(property_data, (chi=χ, T_tilde=t̃, property="specific_heat", value=cp))
                end
                
                # Viscosity
                μ = CoolProp.PropsSI("VISCOSITY", "P", 1e5, "T", T_val, fluid_str)
                if !isnan(μ)
                    push!(property_data, (chi=χ, T_tilde=t̃, property="viscosity", value=μ))
                end
            catch e
                # Skip problematic points
            end
        end
    end
    
    # Convert to DataFrame for saving
    df_properties = DataFrame(property_data)
    
    return Dict(
        "name" => fluid_name,
        "freezing" => df_freezing,
        "max_temps" => df_max_temps,
        "properties" => df_properties
    )
end

# ╔═╡ 7a09b7e1-5a6c-4d03-8b55-36d1d376531f
"""
    save_fluid_data(data_dict::Dict, output_dir::String)

Save fluid interpolation data to Parquet files in the specified directory using long format.
"""
function save_fluid_data(data_dict::Dict, output_dir::String)
    fluid_name = data_dict["name"]
    fluid_dir = joinpath(output_dir, lowercase(fluid_name))
    
    try
        mkpath(fluid_dir)
        
        # Save temperature limits
        write_parquet(joinpath(fluid_dir, "freezing.parquet"), data_dict["freezing"])
        write_parquet(joinpath(fluid_dir, "max_temps.parquet"), data_dict["max_temps"])
        
        # Save properties in long format
        write_parquet(joinpath(fluid_dir, "properties.parquet"), data_dict["properties"])
        
        # Initialize DuckDB and register the parquet files if available
        if USE_DUCKDB[] && init_database()
            con = DuckDB.connect(DB_CONN[])
            
            # Register the parquet files with DuckDB
            DuckDB.execute(con, """
                CREATE OR REPLACE VIEW $(lowercase(fluid_name))_freezing AS 
                SELECT * FROM parquet_scan('$(joinpath(fluid_dir, "freezing.parquet"))');
            """)
            
            DuckDB.execute(con, """
                CREATE OR REPLACE VIEW $(lowercase(fluid_name))_max_temps AS 
                SELECT * FROM parquet_scan('$(joinpath(fluid_dir, "max_temps.parquet"))');
            """)
            
            DuckDB.execute(con, """
                CREATE OR REPLACE VIEW $(lowercase(fluid_name))_properties AS 
                SELECT * FROM parquet_scan('$(joinpath(fluid_dir, "properties.parquet"))');
            """)
        end
        
        return fluid_dir
    catch e
        @error "Failed to save fluid data: $e"
        rethrow()
    end
end

# ╔═╡ aa0c2982-4be1-4ff9-ab86-5708d16ed087
"""
    load_fluid_data(fluid_name::String, data_dir::String)

Load fluid interpolation data from Parquet files in the specified directory.
Uses DuckDB if available, otherwise falls back to DataFrames.
"""
function load_fluid_data(fluid_name::String, data_dir::String)
    fluid_dir = joinpath(data_dir, lowercase(fluid_name))
    
    if !isdir(fluid_dir)
        error("Data for fluid '$fluid_name' not found in $data_dir. Use regenerate_fluid_data(\"$fluid_name\") to create it.")
    end
    
    try
        # Check if files exist
        freezing_path = joinpath(fluid_dir, "freezing.parquet")
        max_temps_path = joinpath(fluid_dir, "max_temps.parquet")
        properties_path = joinpath(fluid_dir, "properties.parquet")
        
        for path in [freezing_path, max_temps_path, properties_path]
            if !isfile(path)
                error("Required data file not found: $path. The fluid data may be corrupted.")
            end
        end
        
        # Use DuckDB if available, otherwise use direct DataFrame loading
        if USE_DUCKDB[] && init_database()
            con = DuckDB.connect(DB_CONN[])
            
            # Register the parquet files with DuckDB
            DuckDB.execute(con, """
                CREATE OR REPLACE VIEW $(lowercase(fluid_name))_freezing AS 
                SELECT * FROM parquet_scan('$(freezing_path)');
            """)
            
            DuckDB.execute(con, """
                CREATE OR REPLACE VIEW $(lowercase(fluid_name))_max_temps AS 
                SELECT * FROM parquet_scan('$(max_temps_path)');
            """)
            
            DuckDB.execute(con, """
                CREATE OR REPLACE VIEW $(lowercase(fluid_name))_properties AS 
                SELECT * FROM parquet_scan('$(properties_path)');
            """)
            
            # Query the tables to get DataFrames
            freezing = DuckDB.execute(con, "SELECT * FROM $(lowercase(fluid_name))_freezing") |> DataFrame
            max_temps = DuckDB.execute(con, "SELECT * FROM $(lowercase(fluid_name))_max_temps") |> DataFrame
            
            # Query and pivot the properties data
            density = DuckDB.execute(con, """
                SELECT chi, T_tilde, value 
                FROM $(lowercase(fluid_name))_properties 
                WHERE property = 'density'
            """) |> DataFrame
            
            thermal_conductivity = DuckDB.execute(con, """
                SELECT chi, T_tilde, value 
                FROM $(lowercase(fluid_name))_properties 
                WHERE property = 'thermal_conductivity'
            """) |> DataFrame
            
            specific_heat = DuckDB.execute(con, """
                SELECT chi, T_tilde, value 
                FROM $(lowercase(fluid_name))_properties 
                WHERE property = 'specific_heat'
            """) |> DataFrame
            
            viscosity = DuckDB.execute(con, """
                SELECT chi, T_tilde, value 
                FROM $(lowercase(fluid_name))_properties 
                WHERE property = 'viscosity'
            """) |> DataFrame
        else
            # Fallback to direct DataFrame loading using Parquet.jl
            freezing = read_parquet(freezing_path) |> DataFrame
            max_temps = read_parquet(max_temps_path) |> DataFrame
            
            # Load properties and filter by property type
            properties = read_parquet(properties_path) |> DataFrame
            
            density = filter(row -> row.property == "density", properties)[:, [:chi, :T_tilde, :value]]
            thermal_conductivity = filter(row -> row.property == "thermal_conductivity", properties)[:, [:chi, :T_tilde, :value]]
            specific_heat = filter(row -> row.property == "specific_heat", properties)[:, [:chi, :T_tilde, :value]] 
            viscosity = filter(row -> row.property == "viscosity", properties)[:, [:chi, :T_tilde, :value]]
        end
        
        return Dict(
            "name" => fluid_name,
            "freezing" => freezing,
            "max_temps" => max_temps,
            "rho" => density,
            "lambda" => thermal_conductivity,
            "cp" => specific_heat,
            "mu" => viscosity
        )
    catch e
        @error "Error loading fluid data for $fluid_name: $e"
        
        # Try to use approximations if available
        approx = default_fluid_approximations(fluid_name)
        if !isnothing(approx)
            @info "Using simplified approximations for $fluid_name properties. For better accuracy, regenerate the fluid data."
            return approx
        else
            rethrow()
        end
    end
end

# ╔═╡ dcea5e9a-0968-48b6-9be1-e5f6f8fe5c2a
"""
    build_interpolants(data::Dict)

Build interpolation functions from the loaded data.
Handles both full datasets and approximation functions.
"""
function build_interpolants(data::Dict)
    fluid_name = data["name"]
    
    # Check if we're using approximation functions
    if get(data, "approximation", false)
        # Create a simplified version of FluidInterpolants that uses the approximation functions
        return FluidInterpolants(
            fluid_name,
            [0.0, 1.0],  # Placeholder for chi_valid
            # Create simple wrapper functions that behave like the interpolants
            build_approximate_interp(data["T_freeze"]), 
            build_approximate_interp(data["T_max"]),
            build_approximate_interp(data["ρ"]),
            build_approximate_interp(data["λ"]), 
            build_approximate_interp(data["cp"]),
            build_approximate_interp(data["μ"])
        )
    end
    
    # Normal interpolation based on data
    # Extract data for 1D interpolations
    freezing_df = data["freezing"]
    max_temps_df = data["max_temps"]
    
    chi_valid = freezing_df.chi
    T_freeze_vals = freezing_df.T_freeze
    T_max_vals = max_temps_df.T_max
    
    # Build 1D interpolants
    itp_Tfreeze = LinearInterpolation(chi_valid, T_freeze_vals, extrapolation_bc=Interpolations.Flat())
    itp_Tmax = LinearInterpolation(chi_valid, T_max_vals, extrapolation_bc=Interpolations.Flat())
    
    # Build 2D interpolants for properties
    function build_2D_interpolant(prop_df)
        # Group data by chi and T_tilde to reconstruct the 2D grid
        chi_unique = sort(unique(prop_df.chi))
        t_tilde_unique = sort(unique(prop_df.T_tilde))
        
        # Create a matrix to hold values
        values_matrix = Array{Float64}(undef, length(chi_unique), length(t_tilde_unique))
        values_matrix .= NaN
        
        # Map from chi, T_tilde to indices
        chi_to_idx = Dict(chi => i for (i, chi) in enumerate(chi_unique))
        t_to_idx = Dict(t => i for (i, t) in enumerate(t_tilde_unique))
        
        # Fill the matrix with values
        for row in eachrow(prop_df)
            i = chi_to_idx[row.chi]
            j = t_to_idx[row.T_tilde]
            values_matrix[i, j] = row.value
        end
        
        # Find rows with valid data
        good_rows = [i for i in 1:size(values_matrix, 1) 
                     if !all(isnan, view(values_matrix, i, :))]
        
        if isempty(good_rows)
            error("No valid data for interpolation")
        end
        
        # Extract valid data
        chi_good = chi_unique[good_rows]
        values_good = values_matrix[good_rows, :]
        
        # Create interpolation
        itp_2D = Interpolations.interpolate((chi_good, t_tilde_unique), values_good, Gridded(Linear()))
        return extrapolate(itp_2D, Interpolations.Flat())
    end
    
    # Build interpolants for each property
    itp_2D_rho = build_2D_interpolant(data["rho"])
    itp_2D_lambda = build_2D_interpolant(data["lambda"])
    itp_2D_cp = build_2D_interpolant(data["cp"])
    itp_2D_mu = build_2D_interpolant(data["mu"])
    
    return FluidInterpolants(
        fluid_name,
        chi_valid,
        itp_Tfreeze,
        itp_Tmax,
        itp_2D_rho,
        itp_2D_lambda,
        itp_2D_cp,
        itp_2D_mu
    )
end

# ╔═╡ 1b699afb-5dec-4941-b8dd-c7a7e08817df
"""
    build_approximate_interp(func)

Create an approximation wrapper that behaves like an interpolation object
but uses the provided function directly.
"""
function build_approximate_interp(func)
    # For 1D functions (T_freeze, T_max)
    if applicable(func, 0.0)
        return x -> func(x)
    else
        # For 2D functions (properties)
        return (x, y) -> func(y, x)  # Note: we swap arguments since our property funcs are (T, χ)
    end
end

# ╔═╡ c43bc4b9-eaba-4204-a0b4-5e000b66e44e
"""
    load_fluid(name::String; data_dir=nothing, regenerate=false)

Load a fluid's interpolants, either from cache or from saved data.
If regenerate is true, regenerate the data from scratch.
"""
function load_fluid(name::String; data_dir=nothing, regenerate=false)
    # Check if already loaded
    if !regenerate && haskey(LOADED_FLUIDS, lowercase(name))
        return LOADED_FLUIDS[lowercase(name)]
    end
    
    # Default to package data directory if not specified
    if isnothing(data_dir)
        data_dir = joinpath(@__DIR__, "..", "data", "fluids")
    end
    
    local interpolants
    
    if regenerate
        # Generate fresh data
        @info "Generating fluid data for $name..."
        fluid_data = generate_fluid_data(name)
        
        # Ensure directory exists
        if !isdir(data_dir)
            try
                mkpath(data_dir)
            catch e
                @warn "Cannot create data directory $data_dir: $e. Will continue with in-memory data only."
            end
        end
        
        # Try to save, but continue even if it fails
        try
            save_fluid_data(fluid_data, data_dir)
        catch e
            @warn "Could not save fluid data: $e. Will continue with in-memory data only."
        end
        
        interpolants = build_interpolants(fluid_data)
    else
        try
            # Try to load from saved data
            @info "Loading fluid data for $name..."
            fluid_data = load_fluid_data(name, data_dir)
            interpolants = build_interpolants(fluid_data)
        catch e
            if occursin("not found", "$e")
                @warn "Data for $name not found. Generating new data..."
                fluid_data = generate_fluid_data(name)
                
                # Try to save, but continue even if it fails
                try
                    save_fluid_data(fluid_data, data_dir)
                catch save_error
                    @warn "Could not save fluid data: $save_error. Will continue with in-memory data only."
                end
                
                interpolants = build_interpolants(fluid_data)
            else
                # For other errors, try using approximations
                @warn "Error loading/generating data for $name: $e"
                approx = default_fluid_approximations(name)
                if !isnothing(approx)
                    @info "Using simplified approximations for $name properties."
                    interpolants = build_interpolants(approx)
                else
                    rethrow()
                end
            end
        end
    end
    
    # Cache the interpolants
    LOADED_FLUIDS[lowercase(name)] = interpolants
    
    # Update the actual property functions with implementations using the loaded data
    setup_property_functions()
    
    return interpolants
end

# ╔═╡ c9adc754-772a-45d6-9f4e-4578b60ad7df
"""
    regenerate_fluid_data(name::String; n_chi=100, n_temp=50, output_dir=nothing)

Regenerate fluid data with specified sampling parameters.
"""
function regenerate_fluid_data(name::String; n_chi=100, n_temp=50, output_dir=nothing)
    if isnothing(output_dir)
        output_dir = joinpath(@__DIR__(), "..", "data", "fluids")
    end
    
    @info "Regenerating fluid data for $name with n_chi=$n_chi, n_temp=$n_temp..."
    fluid_data = generate_fluid_data(name, n_chi=n_chi, n_temp=n_temp)
    
    # Ensure directory exists
    if !isdir(output_dir)
        try
            mkpath(output_dir)
        catch e
            @warn "Cannot create data directory $output_dir: $e. Will continue with in-memory data only."
        end
    end
    
    # Try to save, but continue even if it fails
    try
        save_fluid_data(fluid_data, output_dir)
    catch e
        @warn "Could not save fluid data: $e. Will continue with in-memory data only."
    end
    
    interpolants = build_interpolants(fluid_data)
    
    # Update the cache
    LOADED_FLUIDS[lowercase(name)] = interpolants
    
    # Update the actual property functions with implementations using the loaded data
    setup_property_functions()
    
    return interpolants
end

# ╔═╡ fe3a9822-bcd9-4291-a2c0-3dc1a46ba2cd
"""
    setup_property_functions()

Define the implementations of the fluid property functions using the loaded interpolants.
"""
function setup_property_functions()
    # Define the fluid property functions to use the loaded interpolants
    # These implementations preserve compatibility with Symbolics

    # Define ρ_b_
    @eval function ρ_b_(T::Number, χ::Number; fluid="mea")
        interp = get(LOADED_FLUIDS, lowercase(fluid), nothing)
        if isnothing(interp)
            interp = load_fluid(fluid)
        end
        
        T_f = interp.itp_Tfreeze(χ)
        T_M = interp.itp_Tmax(χ)
        T_tilde_val = (T - T_f) / (T_M - T_f)
        return interp.itp_2D_rho(χ, T_tilde_val)
    end

    # Define λ_b_
    @eval function λ_b_(T::Number, χ::Number; fluid="mea")
        interp = get(LOADED_FLUIDS, lowercase(fluid), nothing)
        if isnothing(interp)
            interp = load_fluid(fluid)
        end
        
        T_f = interp.itp_Tfreeze(χ)
        T_M = interp.itp_Tmax(χ)
        T_tilde_val = (T - T_f) / (T_M - T_f)
        return interp.itp_2D_lambda(χ, T_tilde_val)
    end

    # Define cp_b_
    @eval function cp_b_(T::Number, χ::Number; fluid="mea")
        interp = get(LOADED_FLUIDS, lowercase(fluid), nothing)
        if isnothing(interp)
            interp = load_fluid(fluid)
        end
        
        T_f = interp.itp_Tfreeze(χ)
        T_M = interp.itp_Tmax(χ)
        T_tilde_val = (T - T_f) / (T_M - T_f)
        return interp.itp_2D_cp(χ, T_tilde_val)
    end

    # Define μ_b_
    @eval function μ_b_(T::Number, χ::Number; fluid="mea")
        interp = get(LOADED_FLUIDS, lowercase(fluid), nothing)
        if isnothing(interp)
            interp = load_fluid(fluid)
        end
        
        T_f = interp.itp_Tfreeze(χ)
        T_M = interp.itp_Tmax(χ)
        T_tilde_val = (T - T_f) / (T_M - T_f)
        return interp.itp_2D_mu(χ, T_tilde_val)
    end

    # Define β_b_
    @eval function β_b_(T::Number, χ::Number; fluid="mea")
        interp = get(LOADED_FLUIDS, lowercase(fluid), nothing)
        if isnothing(interp)
            interp = load_fluid(fluid)
        end
        
        # Get the temperature limits for concentration χ
        T_f = interp.itp_Tfreeze(χ)
        T_M = interp.itp_Tmax(χ)
        
        # Map the physical temperature T to the nondimensional coordinate
        T_tilde_val = (T - T_f) / (T_M - T_f)
        
        # Evaluate the density at (T,χ)
        ρ_val = interp.itp_2D_rho(χ, T_tilde_val)
        
        # Compute the gradient of the density interpolation
        grad = Interpolations.gradient(interp.itp_2D_rho, χ, T_tilde_val)
        
        # Extract the derivative with respect to the nondimensional temperature T̃
        dρ_dTtilde = grad[2]
        
        # Convert derivative to physical temperature using the chain rule
        dρ_dT = dρ_dTtilde / (T_M - T_f)
        
        # Return the thermal expansion coefficient
        # Note the negative sign to ensure β is positive
        return -dρ_dT / ρ_val
    end

    # Define single-argument version of β_b_
    @eval β_b_(T::Number; fluid="mea") = β_b_(T, 0.0, fluid=fluid)
end

# ╔═╡ a0bd5972-be0e-4a3d-a2f9-6db2d24a1c4d
"""
    query_fluid_property(fluid_name::String, property_name::String; where_clause::String="")

Query the database to get specific property data with optional filtering.
If DuckDB is not available, uses DataFrame filtering instead.

# Examples
```julia
# Get density values for water at chi=0.5
df = query_fluid_property("water", "density", where_clause="chi=0.5")

# Get viscosity values in a temperature range (T_tilde between 0.2 and 0.4)
df = query_fluid_property("mea", "viscosity", where_clause="T_tilde BETWEEN 0.2 AND 0.4")
```
"""
function query_fluid_property(fluid_name::String, property_name::String; where_clause::String="")
    # Map property aliases to standard names
    property_mapping = Dict(
        "density" => "density",
        "rho" => "density",
        "thermal_conductivity" => "thermal_conductivity",
        "lambda" => "thermal_conductivity",
        "specific_heat" => "specific_heat",
        "cp" => "specific_heat",
        "viscosity" => "viscosity",
        "mu" => "viscosity"
    )
    
    property_key = get(property_mapping, lowercase(property_name), property_name)
    
    # Ensure we have the fluid data loaded
    if !haskey(LOADED_FLUIDS, lowercase(fluid_name))
        load_fluid(fluid_name)
    end
    
    # If DuckDB is available, use it for querying
    if USE_DUCKDB[] && init_database()
        con = DuckDB.connect(DB_CONN[])
        
        query = """
            SELECT chi, T_tilde, value 
            FROM $(lowercase(fluid_name))_properties 
            WHERE property = '$(property_key)'
        """
        
        if !isempty(where_clause)
            query *= " AND $(where_clause)"
        end
        
        result = DuckDB.execute(con, query) |> DataFrame
        return result
    else
        # Fallback to DataFrame filtering - we need to reload the data
        data_dir = joinpath(@__DIR__, "..", "data", "fluids")
        fluid_path = joinpath(data_dir, lowercase(fluid_name), "properties.parquet")
        
        if !isfile(fluid_path)
            error("Cannot query properties: data file not found at $fluid_path")
        end
        
        properties = read_parquet(fluid_path) |> DataFrame
        
        # Filter by property
        result = filter(row -> row.property == property_key, properties)
        
        # Handle where_clause using simple parsing (limited functionality)
        if !isempty(where_clause)
            if occursin("=", where_clause)
                # Simple equality
                field, value = split(where_clause, "=")
                field = strip(field)
                value = parse(Float64, strip(value))
                result = filter(row -> getproperty(row, Symbol(field)) ≈ value, result)
            elseif occursin("BETWEEN", uppercase(where_clause))
                # Simple BETWEEN
                parts = split(uppercase(where_clause), "BETWEEN")
                field = strip(lowercase(parts[1]))
                range_parts = split(parts[2], "AND")
                lower = parse(Float64, strip(range_parts[1]))
                upper = parse(Float64, strip(range_parts[2]))
                result = filter(row -> lower ≤ getproperty(row, Symbol(field)) ≤ upper, result)
            else
                @warn "Complex where_clause not supported with DataFrame filtering: $where_clause"
            end
        end
        
        return result[:, [:chi, :T_tilde, :value]]
    end
end

# ╔═╡ b4b76a5a-1d65-4cd9-ac86-6dc4995bd02e
"""
    close_database()

Close the DuckDB connection if it's open.
"""
function close_database()
    if USE_DUCKDB[] && !isnothing(DB_CONN[])
        try
            close(DB_CONN[])
            DB_CONN[] = nothing
        catch e
            @warn "Error closing database connection: $e"
        end
    end
end

# Preload MEA as the default fluid at module initialization

# ╔═╡ 2070c1dc-ac03-493a-a7aa-8eb7ed8538e3
# ╠═╡ disabled = true
#=╠═╡
begin
	try
        # Initialize database if DuckDB is available
        if USE_DUCKDB[]
			println("!oiuvhnfsdoiv")
            init_database()
        end
        
        # Initialize property functions even before loading any data
        setup_property_functions()
			println("!oiuvhnfsdoiv")
        
        # Check if data directory exists, create if needed
        data_dir = joinpath("C:\\Users\\adrien.vautrin\\Documents\\GitHub\\OptimalGshpFlow.jl", "data", "fluids")
			println(data_dir)
        if !isdir(data_dir)
            try
                mkpath(data_dir)
                @info "Created fluid data directory: $data_dir"
            catch e
                @warn "Could not create data directory: $e. Will use in-memory data only."
            end
        end
        
        # Load the default fluid - use approximation if data loading fails
        try
            load_fluid("mea")
        catch e
            @warn "Could not load MEA fluid data: $e. Using approximations instead."
            approx = default_fluid_approximations("mea")
            if !isnothing(approx)
                LOADED_FLUIDS["mea"] = build_interpolants(approx)
                setup_property_functions()
            end
        end
    catch e
        @warn "Error during FluidProperties initialization: $e"
    end
close_database()
end
  ╠═╡ =#

# ╔═╡ a8114d21-19e1-46d9-befe-eb62f6de6d24
#=╠═╡
using DataFrames
  ╠═╡ =#

# ╔═╡ b1df0328-5f22-45af-8b9c-d30328d4f771
#=╠═╡
using DataFrames
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═14abfdab-a03a-4649-8aa6-14e8a4cea7bc
# ╠═a8114d21-19e1-46d9-befe-eb62f6de6d24
# ╠═fedaaebd-cf14-48fa-a2c9-a528f0c65033
# ╠═f177d2ef-ff07-444e-b353-e2484e7e8b00
# ╠═b1df0328-5f22-45af-8b9c-d30328d4f771
# ╠═122c8036-63f1-4deb-8e1d-4dcf471bd298
# ╠═924b8d55-a8f9-4f1a-b1e0-31601ff811c9
# ╠═ee9f209e-8f89-487d-9755-67d06b42308a
# ╠═5fb03f74-ad6a-4d0b-9c6e-96c9d96132ec
# ╠═2dc9195f-aebc-4407-9802-35cc85cfef85
# ╠═e20e2ca5-9ad4-4eda-8242-229496d0c244
# ╠═7af52b81-cdaa-4947-b684-3cd0b658b0a4
# ╠═7fdb2dc7-af0b-4035-bd9d-15af34f27fd1
# ╠═8c986457-9d64-4be5-b5b0-d39b9702e05c
# ╠═3a793d39-1f74-43f6-950a-b46dcb8a3de8
# ╠═806ebe0f-7be4-4f35-8827-eb111d3a8891
# ╠═d370d15e-07a4-49a0-a6b2-b18bb760d22f
# ╠═64dcc6ea-d858-4c29-ad36-534afdf3c8d8
# ╠═1a4a329b-9002-471e-9013-653d14044bf1
# ╠═b49b1c08-a158-4721-bb2e-d2d2fc90201d
# ╠═2009ccd1-5cec-4a3e-8865-876e0a4e50f8
# ╠═ea8f6eb0-49c2-41a2-b5da-56b59ec858bd
# ╠═e9540434-f7c5-486a-af94-0eda47ddae0d
# ╠═1d1cdfa3-a995-4286-bcf7-8f1f741810df
# ╠═69f773f7-3602-4215-9caa-019406b19268
# ╠═739d9847-5a6a-4f2e-8143-167269d44be4
# ╠═44ab047d-d1c0-4165-9525-b2c80572d796
# ╠═b9689240-5e8b-4e9e-8b9c-fef6c20ddfe8
# ╠═0335e609-cfa7-4961-b85b-7214aeecf92e
# ╠═bc04af84-5e7c-46c1-ab08-0f234d17e440
# ╠═20433cc4-97b7-44d0-958a-aff03fdd12ba
# ╠═1022217d-2e51-4a88-ae93-f7bf47f0959d
# ╠═72aef11f-3a07-43ae-821b-87f48dec2874
# ╠═4d01cde3-4840-4e9c-ac3a-8c36c8a14f7d
# ╠═432e6cac-14f9-4d09-9864-6500b306d5d0
# ╠═426d5e79-f04c-4518-8faf-9832bb896d4f
# ╠═7571b2d2-489b-4c5e-b3c8-6f352fcc3830
# ╠═d16ef005-fa61-4456-862f-0ec0cb5287a0
# ╠═4eaac3ee-290b-4603-a30b-cefd15d528bd
# ╠═7f3d7705-8e7f-46dc-9c52-9fc6978f38d5
# ╠═7a09b7e1-5a6c-4d03-8b55-36d1d376531f
# ╠═aa0c2982-4be1-4ff9-ab86-5708d16ed087
# ╠═dcea5e9a-0968-48b6-9be1-e5f6f8fe5c2a
# ╠═1b699afb-5dec-4941-b8dd-c7a7e08817df
# ╠═c43bc4b9-eaba-4204-a0b4-5e000b66e44e
# ╠═c9adc754-772a-45d6-9f4e-4578b60ad7df
# ╠═fe3a9822-bcd9-4291-a2c0-3dc1a46ba2cd
# ╠═a0bd5972-be0e-4a3d-a2f9-6db2d24a1c4d
# ╠═b4b76a5a-1d65-4cd9-ac86-6dc4995bd02e
# ╠═2070c1dc-ac03-493a-a7aa-8eb7ed8538e3
