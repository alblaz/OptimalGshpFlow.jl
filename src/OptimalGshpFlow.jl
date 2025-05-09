module OptimalGshpFlow

# Export fluid property functions
export ρ_b_, λ_b_, cp_b_, μ_b_, β_b_
export load_fluid, regenerate_fluid_data

# Include the fluid properties module
include("fluid_properties.jl")

using .FluidProperties

# Re-export the fluid property functions
ρ_b_ = FluidProperties.ρ_b_
λ_b_ = FluidProperties.λ_b_
cp_b_ = FluidProperties.cp_b_
μ_b_ = FluidProperties.μ_b_
β_b_ = FluidProperties.β_b_

# Include test functions
include("test_pluto.jl")

# Default concentration for MEA
const χ_ = 0.285

end # module OptimalGshpFlow
    