module FluidProperties

using Interpolations
using CoolProp
using Parquet
using DataFrames
using ModelingToolkit
using FileIO

# Optional DuckDB support
const USE_DUCKDB = Ref(false)
try
    using DuckDB
    global USE_DUCKDB = Ref(true)
catch
    @warn "DuckDB not available. Long-format queries will use DataFrames filtering instead."
end

export FluidInterpolants, load_fluid, regenerate_fluid_data
export ρ_b_, λ_b_, cp_b_, μ_b_, β_b_
export list_available_fluids, query_fluid_property

# Global variable to store loaded fluid interpolants
const LOADED_FLUIDS = Dict{String, Any}()
# DuckDB connection for fluid data querying
const DB_CONN = Ref{Any}(nothing)

# Pre-declare the symbolic versions of the fluid property functions
# These will be properly defined once data is loaded
"""
    ρ_b_(T, χ)

Calculate the density of the fluid at temperature T and concentration χ.
"""
function ρ_b_ end

"""
    λ_b_(T, χ)

Calculate the thermal conductivity of the fluid at temperature T and concentration χ.
"""
function λ_b_ end

"""
    cp_b_(T, χ)

Calculate the specific heat of the fluid at temperature T and concentration χ.
"""
function cp_b_ end

"""
    μ_b_(T, χ)

Calculate the viscosity of the fluid at temperature T and concentration χ.
"""
function μ_b_ end

"""
    β_b_(T, χ)

Calculate the thermal expansion coefficient of the fluid at temperature T and concentration χ.
"""
function β_b_ end

"""
    β_b_(T)

Calculate the thermal expansion coefficient with default concentration (0.0).
"""
β_b_(T) = β_b_(T, 0.0)

# Register the symbolic functions
ModelingToolkit.@register_symbolic ρ_b_(T, χ)
ModelingToolkit.@register_symbolic λ_b_(T, χ)
ModelingToolkit.@register_symbolic cp_b_(T, χ)
ModelingToolkit.@register_symbolic μ_b_(T, χ)
ModelingToolkit.@register_symbolic β_b_(T, χ)
ModelingToolkit.@register_symbolic β_b_(T)

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

"""
    regenerate_fluid_data(name::String; n_chi=100, n_temp=50, output_dir=nothing)

Regenerate fluid data with specified sampling parameters.
"""
function regenerate_fluid_data(name::String; n_chi=100, n_temp=50, output_dir=nothing)
    if isnothing(output_dir)
        output_dir = joinpath(@__DIR__, "..", "data", "fluids")
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
function __init__()
    try
        # Initialize database if DuckDB is available
        if USE_DUCKDB[]
            init_database()
        end
        
        # Initialize property functions even before loading any data
        setup_property_functions()
        
        # Check if data directory exists, create if needed
        data_dir = joinpath(@__DIR__, "..", "data", "fluids")
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
    
    # Set up finalizer to close database connection
    atexit(() -> close_database())
end

end # module
