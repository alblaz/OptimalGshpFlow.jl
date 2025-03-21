"""
    `compressor_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::Refrigerant)``

Returns the power delivered by the compressor in W.

`Tc`is the condensing temperature in Kelvin,
`ηis` is the isotropic efficiency of the compressor, 
`evaporator_solution::EvaporatorSolution` contains the results from the evaporator simulation
`refrigerant::Refrigerant` 
"""
function compressor_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::Refrigerant)
    
    Vref = evaporator_solution.twophase_region.Vc
    T_compressor_inlet = evaporator_solution.superheat_region.Tch

    pc = PropsSI("P", "T",  Tc, "Q", 0, refrigerant.name)  #condensation pressure
    pe = PropsSI("P", "T",  evaporator_solution.twophase_region.Tcl, "Q", 1, refrigerant.name) #evaporation pressure

    h_compressor_inlet = PropsSI("H", "T",  T_compressor_inlet, "P", pe, refrigerant.name)
    s_compressor_inlet = PropsSI("S", "T",  T_compressor_inlet, "P", pe, refrigerant.name)
    his_compressor_outlet = PropsSI("H", "S",  s_compressor_inlet, "P", pc, refrigerant.name)
    Δh_compressor = (his_compressor_outlet - h_compressor_inlet)/ηis

    ρ =  PropsSI("D", "T",  T_compressor_inlet, "P", pe, refrigerant.name) 
    
    W = ρ * Vref * Δh_compressor
    return W
end

"""
    `condenser_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::Refrigerant)``

Returns the power delivered at the condenser in W.

`Tc`is the condensing temperature in Kelvin,
`ηis` is the isotropic efficiency of the compressor, 
`evaporator_solution::EvaporatorSolution` contains the results from the evaporator simulation
`refrigerant::Refrigerant` 
"""
function condenser_power(Tc, ηis, evaporator_solution::EvaporatorSolution, refrigerant::Refrigerant)
    
    Vref = evaporator_solution.twophase_region.Vc
    T_compressor_inlet = evaporator_solution.superheat_region.Tch

    pc = PropsSI("P", "T",  Tc, "Q", 0, refrigerant.name)  #condensation pressure
    pe = PropsSI("P", "T",  evaporator_solution.twophase_region.Tcl, "Q", 1, refrigerant.name) #evaporation pressure

    h_compressor_inlet = PropsSI("H", "T",  T_compressor_inlet, "P", pe, refrigerant.name)
    s_compressor_inlet = PropsSI("S", "T",  T_compressor_inlet, "P", pe, refrigerant.name)
    his_compressor_outlet = PropsSI("H", "S",  s_compressor_inlet, "P", pc, refrigerant.name)
    Δh_compressor = (his_compressor_outlet - h_compressor_inlet)/ηis

    h_compressor_outlet =  h_compressor_inlet + Δh_compressor
    h_condenser_outlet =  PropsSI("H", "T",  Tc, "Q", 0, refrigerant.name) # add subcooling
    
    ρ =  PropsSI("D", "T",  T_compressor_inlet, "P", pe, refrigerant.name) 
    Q1 = ρ * Vref * (h_compressor_outlet - h_condenser_outlet) # add volumetric efficiency

    return Q1
end