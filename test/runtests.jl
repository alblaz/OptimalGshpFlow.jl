using OptimalGshpFlow
using Test
using CoolProp

include("./helpers.jl")        # shared ForwardDiff tag (compile time)
include("./utilities.jl")
include("./transforms.jl")
include("./fluidproperties.jl")
include("./pipeflow.jl")
include("./pressuredrop.jl")
include("./refrigerant_pressuredrop.jl")
include("./borehole.jl")
include("./borehole_coupling.jl")
include("./ground.jl")
include("./properties_clapeyron.jl")
include("./property_interface.jl")
include("./tv_state.jl")
include("./components.jl")
include("./operating_envelope.jl")
include("./solutiontypes.jl")
include("./compressor.jl")
# include("./models.jl")

@testset "OptimalGshpFlow.jl" begin
    # Write your tests here
    
    
end
