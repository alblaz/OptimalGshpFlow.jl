module OptimalGshpFlow

using CoolProp, Interpolations
# Write your package code here.
using StaticArrays, Parameters, ExponentialUtilities
using RequiredInterfaces
using LambertW
using SpecialFunctions: expint
import Clapeyron, EntropyScaling, ForwardDiff

include("utilities.jl")
export linearmap_01_to_ab, linearmap_ab_to_01
export tanhmap_01_to_R, tanhmap_R_to_01, tanhmap_ab_to_R, tanhmap_R_to_ab
export sigmoid_01_to_R, sigmoid_R_to_01, sigmoid_ab_to_R, sigmoid_R_to_ab

include("transforms.jl")
export VariableTransform, Free, Linear, Unit, Box, LogBox, Above, from_R, to_R, bounds

include("properties/interface.jl")
export AbstractRefrigerant, critical_pressure, molar_mass, pressure_transform
export T_dew, T_bubble, h_dew, h_bubble, ρ_dew, ρ_bubble, s_dew, cp_dew, μ_dew, λ_dew, surface_tension
export h_pT, s_pT, ρ_pT, cp_pT, μ_pT, λ_pT
export h_superheated, s_superheated, ρ_superheated, cp_superheated, μ_superheated, λ_superheated

include("properties/coolprop_state.jl")
include("fluidproperties.jl")
# export interpolate_sf_properties, interpolate_ref_properties_twophaseregion, interpolate_ref_properties_superheatedregion
export SecondaryFluid, freezing_point
# export secondary_fluid
export Refrigerant
#export  RefrigerantSaturatedVapor, RefrigerantSuperHeatedVapor, refrigerant_twophaseregion, refrigerant_superheated

include("properties/clapeyron_refrigerant.jl")
include("properties/tv_state.jl")
export volume_pT, p_Tv, h_Tv, s_Tv, cp_Tv, ρ_Tv, μ_Tv, λ_Tv, pressure_residual
export volume_transform, volume_to_R, volume_from_R, saturation_T, p_sat
export ClapeyronRefrigerant, TransportBackend, EntropyScalingTransport, CoolPropTransport
export cp_twophase, T_ps, h_ps

include("evaporatordesign.jl")
include("pipeflow.jl")        # before borehole.jl: defines NusseltTransition
include("borehole.jl")
include("plateheatexchanger.jl")
include("solutiontypes.jl")
include("heatexchanger.jl")
include("evaporator_formulations.jl")
include("models.jl")
include("compressor.jl")
include("pressuredrop.jl")
include("refrigerant_pressuredrop.jl")
include("ground.jl")
include("components/core.jl")
include("components/evaporator.jl")

# from "solutiontypes.jl"
export HeatExchangerSolution, BoreholeSolution, EvaporatorSolution, EvaporatorBoreholeSolution

# from "borehole.jl"
export Borehole, SingleUPipeBorehole

# from "pipeflow.jl"
export nusselt_circularpipe, h_circularpipe
export NusseltTransition, SmoothTransition, SharpTransition, GnielinskiTransition

# from "plateheatexchanger.jl"
export PlateHeatExchanger, h_platehex_wanniarachchi
export counterflow_heatexchanger

#from "heatexchanger.jl"
export temperatureprofiles

# from "evaporatordesign.jl"
export EvaporatorDesign

# from "pressuredrop.jl"
export  pressure_drop_evap,  pressure_drop, pumping_power

# include("models.jl")
export evaporator_model!, evaporator_model, evaporatorborehole_model!, evaporatorborehole_model, borehole_model!
export EvaporatorParameters, BoreholeParameters, EvaporatorBoreholeParameters
export PropertyFormulation, NestedPressure, SaturationTemperature, ExplicitTv
export n_extra_unknowns, temperature_bounds, temperature_transform, evaporator_transforms
export evaporator_unknowns, nested_unknowns

export compressor_power, condenser_power, compressor_condenser_power

# from "components/"
export AbstractComponent, unknowns, ports, residual_scales, residuals!, guess
export ComponentSystem, residual!, initial_guess, physical, port_values, validity, validity_issues, continuation, set_fixed
export Evaporator, BoreholeHX, default_formulation, default_Te_min, nominal_duty, evaporator_solution, borehole_solution
export pressure_drop, pumping_power
export TwoPhaseFriction, Amalfi2016, KineticEnergyFriction, LONGO_FITS, twophase_pressure_drop, friction_pressure_drop
export acceleration_pressure_drop, gravity_pressure_drop, port_pressure_drop, vapour_pressure_drop, homogeneous_density, channel_mass_flux
export ILS, InfiniteLineSource, borehole_wall_temperature, UniformLoadHistory, push_load!

__init__() = _init_coolprop_indices()

end
