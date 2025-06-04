### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 92fc6816-440c-43ca-be74-e2df5dd428e6
using Printf

# ╔═╡ 2d2cdfcb-8944-4327-908e-be40a4ebd7c4
using Parquet

# ╔═╡ f1e7350d-8d20-4165-ae23-2a9fe1d9ef4b
using FileIO

# ╔═╡ f789f2bd-46e5-4c3c-8994-21e2641bd990
using DuckDB

# ╔═╡ 9573d9fc-09d4-42b1-b15a-ef24de6fe683
using SpecialFunctions

# ╔═╡ 5293784d-1b83-4438-846b-09b501cb0130
using Interpolations

# ╔═╡ 8aa581a6-41e1-4c9c-b354-1240b0905ec3
using Clapeyron

# ╔═╡ f382d83f-c022-4ed8-a2ba-c4a6d5eaec7c
using ModelingToolkit

# ╔═╡ 5c89d2d9-0b23-4f2d-a11f-e99d2e9fccb5
using LinearAlgebra

# ╔═╡ 3ddca4e8-cb44-4009-92ba-1b6d4e84bd46
using DataFrames

# ╔═╡ e8d6642a-d5dc-4687-9e33-f94589bc89db
using CoolProp

# Get the path from the parent module if possible

# ╔═╡ aafe8e6b-a7db-4ef4-a359-e990ff367a04
# Get the path from the parent module if possible
const REF_DB_PATH = @isdefined(Utils) && isdefined(Utils, :REFRIGERANT_DATA_PATH) ? 
                   Utils.REFRIGERANT_DATA_PATH : 
                   joinpath(dirname(dirname(dirname(@__FILE__))), "data", "refrigerants")

# ╔═╡ d713fddc-7512-4aed-b484-086c1f019198
const REFRIGERANT_METADATA_FILE = joinpath(REF_DB_PATH, "refrigerants_metadata.parquet")

# ╔═╡ e2f17d36-ec9f-4ae1-b42e-1dabc14bfb0e
const REFRIGERANT_PROPERTIES_FILE = joinpath(REF_DB_PATH, "refrigerants_properties.parquet")

# ╔═╡ 42393a53-fad6-436a-9658-c57e5504c9ab
const DEFAULT_REFRIGERANT = "R134a"

# ╔═╡ 551fdee9-bfc7-41cb-b34f-55f0c791961e
# Cache for loaded refrigerant interpolants - similar to LOADED_FLUIDS
const LOADED_REFRIGERANTS = Dict{String, Dict{String, Any}}()

# ╔═╡ efabaacf-c212-4880-b484-108658efbe10
function mkdir_if_not_exists(path)
    if !isdir(path)
        mkpath(path)
    end
end

# ╔═╡ 0c8c016c-e80e-11ee-1fee-50bdde504b06
function ensure_refrigerant_directories()
    # Make sure the main data directory exists
    mkdir_if_not_exists(REF_DB_PATH)
    
    # Create any necessary subdirectories here if needed
    
    return REF_DB_PATH
end

# Initialize directories when loading the module
ensure_refrigerant_directories()

# ╔═╡ e69cc47c-5860-4f42-afb8-caf81d136ed4
function get_property_ranges(refrigerant::String)
    # Retrieve critical properties
    T_crit = CoolProp.PropsSI("Tcrit", refrigerant)
    p_crit = CoolProp.PropsSI("pcrit", refrigerant)
    T_triple = CoolProp.PropsSI("Ttriple", refrigerant)
    
    # Define ranges
    T_range = LinRange(T_triple, 0.95 * T_crit, 100)
    p_range = 10.0 .^ LinRange(log10(1e4), log10(0.95 * p_crit), 100)
    
    return T_range, p_range
end

# ╔═╡ 8cdf46e3-8331-4a48-9c74-8ae7f38b11de
function create_refrigerant_interpolant(refrigerant::String, property::String)
    T_range, p_range = get_property_ranges(refrigerant)
    
    # Create a grid of property values
    prop_data = zeros(length(T_range), length(p_range))
    
    # For each temperature and pressure, compute the property
    for (i, T) in enumerate(T_range)
        for (j, p) in enumerate(p_range)
            try
                prop_data[i, j] = CoolProp.PropsSI(property, "T", T, "P", p, refrigerant)
            catch e
                # Handle errors (e.g., beyond saturation curve)
                prop_data[i, j] = NaN
            end
        end
    end
    
    # Replace NaNs with valid neighboring values
    for i in 1:size(prop_data, 1)
        for j in 1:size(prop_data, 2)
            if isnan(prop_data[i, j])
                # Find nearest non-NaN value
                valid_indices = findall(!isnan, prop_data)
                if !isempty(valid_indices)
                    # Fixed: Properly access components of CartesianIndex
                    distances = [sqrt((i-idx[1])^2 + (j-idx[2])^2) for idx in valid_indices]
                    _, min_idx = findmin(distances)
                    best_idx = valid_indices[min_idx]
                    prop_data[i, j] = prop_data[best_idx]  # Use CartesianIndex directly
                else
                    prop_data[i, j] = 0.0  # Fallback value if all are NaN
                end
            end
        end
    end
    
    return prop_data, T_range, p_range
end

# ╔═╡ 7e9019cc-6c93-4d1c-b13e-a612ea7e4dac
function save_to_consolidated_database(metadata_df::DataFrame, properties_df::DataFrame)
    # Check if files exist and load existing data if they do
    if isfile(REFRIGERANT_METADATA_FILE)
        existing_metadata = Parquet.read_parquet(REFRIGERANT_METADATA_FILE) |> DataFrame
        # Remove any existing data for this refrigerant to avoid duplication
        filter!(row -> row.refrigerant != metadata_df.refrigerant[1], existing_metadata)
        # Append new metadata
        metadata_df = vcat(existing_metadata, metadata_df)
    end
    
    if isfile(REFRIGERANT_PROPERTIES_FILE)
        existing_properties = Parquet.read_parquet(REFRIGERANT_PROPERTIES_FILE) |> DataFrame
        # Remove any existing data for this refrigerant to avoid duplication
        filter!(row -> row.refrigerant != properties_df.refrigerant[1], existing_properties)
        # Append new properties
        properties_df = vcat(existing_properties, properties_df)
    end
    
    # Save consolidated data
    Parquet.write_parquet(REFRIGERANT_METADATA_FILE, metadata_df)
    Parquet.write_parquet(REFRIGERANT_PROPERTIES_FILE, properties_df)
    
    println("Saved refrigerant data to consolidated database")
end

# ╔═╡ 819722ed-55c3-4383-bc03-8445b5b22291
function generate_refrigerant_data(refrigerant::String)
    mkdir_if_not_exists(REF_DB_PATH)
    
    properties = [
        "D", # Density
        "H", # Enthalpy
        "S", # Entropy
        "C", # Specific heat capacity
        "L", # Thermal conductivity
        "M"  # Dynamic viscosity
    ]
    
    property_data = Dict()
    T_range = nothing
    p_range = nothing
    
    for prop in properties
        println("Generating data for $refrigerant property: $prop")
        data, T_range_new, p_range_new = create_refrigerant_interpolant(refrigerant, prop)
        
        # Use the first property's ranges for all
        if isnothing(T_range)
            T_range = T_range_new
            p_range = p_range_new
        end
        
        property_data[prop] = data
    end
    
    # Also generate saturation curves
    T_sat_range = LinRange(CoolProp.PropsSI("Ttriple", refrigerant), 
                         0.99 * CoolProp.PropsSI("Tcrit", refrigerant), 100)
    
    p_sat_liq = zeros(length(T_sat_range))
    p_sat_vap = zeros(length(T_sat_range))
    
    for (i, T) in enumerate(T_sat_range)
        try
            p_sat_liq[i] = CoolProp.PropsSI("P", "T", T, "Q", 0, refrigerant)
            p_sat_vap[i] = CoolProp.PropsSI("P", "T", T, "Q", 1, refrigerant)
        catch e
            p_sat_liq[i] = NaN
            p_sat_vap[i] = NaN
        end
    end
    
    # Create metadata dataframe
    metadata_df = DataFrame(
        refrigerant = refrigerant,
        T_min = minimum(T_range),
        T_max = maximum(T_range),
        T_count = length(T_range),
        p_min = minimum(p_range),
        p_max = maximum(p_range),
        p_count = length(p_range),
        T_sat_min = minimum(T_sat_range),
        T_sat_max = maximum(T_sat_range),
        T_sat_count = length(T_sat_range)
    )
    
    # Create properties dataframe in tall format
    property_rows = []
    
    # Process 2D properties (D, H, S, C, L, M)
    for property in properties
        for (i, T) in enumerate(T_range)
            for (j, p) in enumerate(p_range)
                value = property_data[property][i, j]
                push!(property_rows, (
                    refrigerant = refrigerant,
                    property = property,
                    T = T,
                    p = p,
                    value = value
                ))
            end
        end
    end
    
    # Process saturation data
    for (i, T) in enumerate(T_sat_range)
        push!(property_rows, (
            refrigerant = refrigerant,
            property = "p_sat_liq", 
            T = T,
            p = NaN,
            value = p_sat_liq[i]
        ))
        
        push!(property_rows, (
            refrigerant = refrigerant,
            property = "p_sat_vap",
            T = T,
            p = NaN,
            value = p_sat_vap[i]
        ))
    end
    
    properties_df = DataFrame(property_rows)
    
    # Save the data to consolidated files
    save_to_consolidated_database(metadata_df, properties_df)
    
    # Return in-memory format for immediate use
    return (
        name = refrigerant,
        T_range = T_range, 
        p_range = p_range, 
        T_sat_range = T_sat_range,
        data_D = property_data["D"],
        data_H = property_data["H"],
        data_S = property_data["S"],
        data_C = property_data["C"],
        data_L = property_data["L"],
        data_M = property_data["M"],
        p_sat_liq = p_sat_liq,
        p_sat_vap = p_sat_vap
    )
end

# ╔═╡ cbe11304-22ee-4eea-8d73-0ec9464f020a
function load_refrigerant_data(refrigerant::String)
    refrigerant = lowercase(refrigerant)
    
    if !isfile(REFRIGERANT_METADATA_FILE) || !isfile(REFRIGERANT_PROPERTIES_FILE)
        println("Generating data for $refrigerant...")
        return generate_refrigerant_data(refrigerant)
    else
        println("Loading data for $refrigerant from consolidated database...")
        
        # Load metadata
        all_metadata = Parquet.read_parquet(REFRIGERANT_METADATA_FILE) |> DataFrame
        refrigerant_metadata = filter(row -> lowercase(row.refrigerant) == refrigerant, all_metadata)
        
        if isempty(refrigerant_metadata)
            println("Refrigerant $refrigerant not found in database, generating data...")
            return generate_refrigerant_data(refrigerant)
        end
        
        metadata = first(refrigerant_metadata)
        
        # Reconstruct ranges
        T_range = LinRange(metadata.T_min, metadata.T_max, metadata.T_count)
        p_range = exp10.(LinRange(log10(metadata.p_min), log10(metadata.p_max), metadata.p_count))
        T_sat_range = LinRange(metadata.T_sat_min, metadata.T_sat_max, metadata.T_sat_count)
        
        # Load properties
        all_properties = Parquet.read_parquet(REFRIGERANT_PROPERTIES_FILE) |> DataFrame
        refrigerant_properties = filter(row -> lowercase(row.refrigerant) == refrigerant, all_properties)
        
        # Reconstruct data matrices and vectors
        property_data = Dict()
        
        # Initialize matrices for each property
        for prop in ["D", "H", "S", "C", "L", "M"]
            property_data[prop] = zeros(length(T_range), length(p_range))
        end
        
        p_sat_liq = zeros(length(T_sat_range))
        p_sat_vap = zeros(length(T_sat_range))
        
        # Fill in data
        for row in eachrow(refrigerant_properties)
            if row.property in ["D", "H", "S", "C", "L", "M"]
                # 2D property - find indices
                i = findfirst(T -> isapprox(T, row.T), T_range)
                j = findfirst(p -> isapprox(p, row.p), p_range)
                if !isnothing(i) && !isnothing(j)
                    property_data[row.property][i, j] = row.value
                end
            elseif row.property == "p_sat_liq"
                i = findfirst(T -> isapprox(T, row.T), T_sat_range)
                if !isnothing(i)
                    p_sat_liq[i] = row.value
                end
            elseif row.property == "p_sat_vap"
                i = findfirst(T -> isapprox(T, row.T), T_sat_range)
                if !isnothing(i)
                    p_sat_vap[i] = row.value
                end
            end
        end
        
        return (
            name = metadata.refrigerant,
            T_range = T_range, 
            p_range = p_range, 
            T_sat_range = T_sat_range,
            data_D = property_data["D"],
            data_H = property_data["H"],
            data_S = property_data["S"],
            data_C = property_data["C"],
            data_L = property_data["L"],
            data_M = property_data["M"],
            p_sat_liq = p_sat_liq,
            p_sat_vap = p_sat_vap
        )
    end
end

# ╔═╡ 0c8c016c-e80e-11ee-1fee-50bdde504b06
"""
    initialize_refrigerant(refrigerant::String)

Load or generate refrigerant data and create interpolants for various properties.
Returns a dictionary of interpolant functions for the specified refrigerant.
"""
function initialize_refrigerant(refrigerant::String)
    refrigerant = lowercase(refrigerant)
    
    # Check if already loaded
    if haskey(LOADED_REFRIGERANTS, refrigerant)
        return LOADED_REFRIGERANTS[refrigerant]
    end
    
    # Load raw refrigerant data
    data = load_refrigerant_data(refrigerant)
    
    # Create interpolants
    # 1. Create interpolation objects for 2D properties (D, H, S, C, L, M)
    T_range = data.T_range
    p_range = data.p_range
    
    # Create interpolation objects for each property
    D_interp = LinearInterpolation((T_range, p_range), data.data_D)
    H_interp = LinearInterpolation((T_range, p_range), data.data_H)
    S_interp = LinearInterpolation((T_range, p_range), data.data_S)
    C_interp = LinearInterpolation((T_range, p_range), data.data_C)
    L_interp = LinearInterpolation((T_range, p_range), data.data_L)
    M_interp = LinearInterpolation((T_range, p_range), data.data_M)
    
    # 2. Create interpolation objects for saturation curves
    T_sat_range = data.T_sat_range
    p_sat_liq_interp = LinearInterpolation(T_sat_range, data.p_sat_liq)
    p_sat_vap_interp = LinearInterpolation(T_sat_range, data.p_sat_vap)
    
    # Store interpolants in cache
    interp_dict = Dict(
        "D" => (T, p) -> D_interp(T, p),
        "H" => (T, p) -> H_interp(T, p),
        "S" => (T, p) -> S_interp(T, p),
        "C" => (T, p) -> C_interp(T, p),
        "L" => (T, p) -> L_interp(T, p),
        "M" => (T, p) -> M_interp(T, p),
        "p_sat_liq" => (T) -> p_sat_liq_interp(T),
        "p_sat_vap" => (T) -> p_sat_vap_interp(T)
    )
    
    LOADED_REFRIGERANTS[refrigerant] = interp_dict
    
    println("Initialized refrigerant: $refrigerant")
    
    return interp_dict
end

# ╔═╡ 0e194b86-afcf-4390-a32d-cafda8a29905
function list_available_refrigerants()
    if isfile(REFRIGERANT_METADATA_FILE)
        metadata = Parquet.read_parquet(REFRIGERANT_METADATA_FILE) |> DataFrame
        return unique(metadata.refrigerant)
    else
        return String[]
    end
end

# ╔═╡ 52e0f9f1-e6f0-4fd5-8067-d23ce2fa430a
# ╠═╡ disabled = true
#=╠═╡
list_available_refrigerants() |> length
  ╠═╡ =#

# ╔═╡ b41c17c5-0325-424e-b968-2b1c2a7246a9
"""
    setup_refrigerant_functions()

Define the implementations of the refrigerant property functions using the loaded interpolants.
"""
function setup_refrigerant_functions()
    # Define the ρ_r_ function (density)
    @eval function ρ_r_(T::Number, p::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["D"](T, p)
    end

    # Define h_r_ function (enthalpy)
    @eval function h_r_(T::Number, p::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["H"](T, p)
    end

    # Define s_r_ function (entropy)
    @eval function s_r_(T::Number, p::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["S"](T, p)
    end

    # Define cp_r_ function (specific heat)
    @eval function cp_r_(T::Number, p::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["C"](T, p)
    end

    # Define λ_r_ function (thermal conductivity)
    @eval function λ_r_(T::Number, p::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["L"](T, p)
    end

    # Define μ_r_ function (dynamic viscosity)
    @eval function μ_r_(T::Number, p::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["M"](T, p)
    end

    # Define p_sat_liq_r_ function (saturation pressure - liquid)
    @eval function p_sat_liq_r_(T::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["p_sat_liq"](T)
    end

    # Define p_sat_vap_r_ function (saturation pressure - vapor)
    @eval function p_sat_vap_r_(T::Number; refrigerant=DEFAULT_REFRIGERANT)
        interp = get(LOADED_REFRIGERANTS, lowercase(refrigerant), nothing)
        if isnothing(interp)
            interp = initialize_refrigerant(refrigerant)
        end
        return interp["p_sat_vap"](T)
    end
end

# ╔═╡ 77d4e548-3823-4e69-b347-a977aa09b2da
"""
    initialize_refrigerant_properties()

Initialize refrigerant property functions and load the default refrigerant.
"""
function initialize_refrigerant_properties()
    # Define the refrigerant property functions implementations
    setup_refrigerant_functions()
    
    # Register with ModelingToolkit
    #register_refrigerant_symbols()
    
    # Load default refrigerant
    if !haskey(LOADED_REFRIGERANTS, lowercase(DEFAULT_REFRIGERANT))
        try
            initialize_refrigerant(DEFAULT_REFRIGERANT)
        catch e
            @warn "Could not load refrigerant data: $e"
        end
    end
end

# ╔═╡ 2c78b4c3-e214-4a8e-bbcb-6dca96a5409a
# Register density function
ModelingToolkit.@register_symbolic ρ_r_(T, p)

# ╔═╡ c47149cc-e74c-4c0f-b1a7-3dc85f90be58
# Register enthalpy function
ModelingToolkit.@register_symbolic h_r_(T, p)

# ╔═╡ 91578bc9-0d3b-4c17-aaec-e591ec1d6804
# Register entropy function
ModelingToolkit.@register_symbolic s_r_(T, p)

# ╔═╡ d1bb4eea-9d01-4858-b71a-906ef576a753
# Register specific heat function
ModelingToolkit.@register_symbolic cp_r_(T, p)

# ╔═╡ c0073343-c62a-4e8f-b4a3-717113312601
# Register thermal conductivity function
ModelingToolkit.@register_symbolic λ_r_(T, p)

# ╔═╡ 7a46985b-7239-42eb-abe2-6f39fdac46a6
# Register dynamic viscosity function
ModelingToolkit.@register_symbolic μ_r_(T, p)

# ╔═╡ 4d5b8247-73ab-4006-80c9-a82a2a2609cf
# Register saturation pressure functions
ModelingToolkit.@register_symbolic p_sat_liq_r_(T)

# ╔═╡ 38b5ea3f-12da-4851-8379-2f5d0ea641b2
ModelingToolkit.@register_symbolic p_sat_vap_r_(T)

# ╔═╡ 0bccffd2-5e41-40d4-b088-a15a4cd20a3e
initialize_refrigerant_properties()
