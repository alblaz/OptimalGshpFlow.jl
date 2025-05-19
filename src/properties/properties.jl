module Properties

# Import utility functions from Utils
import ..Utils: ρ_b_, λ_b_, cp_b_, μ_b_, β_b_

# Include utility functions
include("pluto_incompressible.jl")

# Export all functions that should be accessible outside the module
export *

end # module Properties