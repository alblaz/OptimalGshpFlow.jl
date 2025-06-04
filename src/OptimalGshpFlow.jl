module OptimalGshpFlow

# Include and export the Utils module
include("utils/Utils.jl")
Utils.initialize_properties()
Utils.initialize_refrigerant_properties()
# Include and export the Properties module
include("properties/Properties.jl")
include("hydraulics/Hydraulics.jl")
#export Properties

# Include other modules as needed
#include("hydraulics/hydraulics.jl")
#include("thermodynamic_cycle/Thermodynamic_cycle.jl")

end # module OptimalGshpFlow