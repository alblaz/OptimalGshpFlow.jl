module Utils

# Define common data paths
const DATA_ROOT = joinpath(dirname(dirname(dirname(@__FILE__))), "data")
const INCOMPRESSIBLE_DATA_PATH = joinpath(DATA_ROOT, "fluids")
const REFRIGERANT_DATA_PATH = joinpath(DATA_ROOT, "refrigerants")

# Create directories if they don't exist
function ensure_data_dirs()
    for dir in [DATA_ROOT, INCOMPRESSIBLE_DATA_PATH, REFRIGERANT_DATA_PATH]
        if !isdir(dir)
            mkpath(dir)
        end
    end
end

# Initialize data directories
ensure_data_dirs()

include("pluto_incompressible_interpolants.jl")
include("pluto_refrigerant_interpolants.jl")

# Re-export specific functions from both modules
export 
    # Incompressible fluids functions
    ρ_b_, λ_b_, cp_b_, μ_b_, β_b_,
    load_fluid, generate_fluid_data, list_available_fluids,
    # Refrigerant functions
    ρ_r_, λ_r_, cp_r_, μ_r_, h_r_, s_r_,
    p_sat_liq_, p_sat_vap_, T_sat_,
    load_refrigerant, generate_refrigerant_data, list_available_refrigerants

end # module Utils