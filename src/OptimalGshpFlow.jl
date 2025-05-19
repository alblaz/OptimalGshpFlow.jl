module OptimalGshpFlow

# Include and export the Utils module
include("utils/Utils.jl")
#export Utils
Utils.initialize_properties()
# Include and export the Properties module
include("properties/Properties.jl")
include("hydraulics/Hydraulics.jl")
#export Properties

# Include other modules as needed
#include("hydraulics/hydraulics.jl")
include("thermodynamic_cycle/Thermodynamic_cycle.jl")

end # module OptimalGshpFlow