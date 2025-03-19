module OptimalGshpFlow

using CoolProp, Interpolations
# Write your package code here.
using StaticArrays, Parameters, ExponentialUtilities
using RequiredInterfaces


include("utilities.jl")
export linearmap_01_to_ab, linearmap_ab_to_01
export tanhmap_01_to_R, tanhmap_R_to_01, tanhmap_ab_to_R, tanhmap_R_to_ab
export sigmoid_01_to_R, sigmoid_R_to_01

include("fluidproperties.jl")
# export interpolate_sf_properties, interpolate_ref_properties_twophaseregion, interpolate_ref_properties_superheatedregion
export SecondaryFluid 
# export secondary_fluid
export Refrigerant
#export  RefrigerantSaturatedVapor, RefrigerantSuperHeatedVapor, refrigerant_twophaseregion, refrigerant_superheated

include("pipeflow.jl")
export nusselt_circularpipe, h_circularpipe 

include("plateheatexchanger.jl")
export PlateHeatExchanger, h_platehex_wanniarachchi

include("evaporatordesign.jl")
export EvaporatorDesign 
#export evaporator_design

include("borehole.jl")
export Borehole, SingleUPipeBorehole
#export uniform_Tb_coeffs

include("solutiontypes.jl")
export HeatExchangerSolution, BoreholeSolution, EvaporatorSolution, EvaporatorBoreholeSolution

include("models.jl")
export evaporator_model!, evaporator_model, evaporatorborehole_model!, evaporatorborehole_model, borehole_model!
export EvaporatorParameters, BoreholeParameters, EvaporatorBoreholeParameters

end
