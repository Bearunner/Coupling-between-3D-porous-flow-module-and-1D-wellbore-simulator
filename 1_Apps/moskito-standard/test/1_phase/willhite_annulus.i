# Over-all heat transfer coefficients in Steam and hot water injection wells
# Willhite, G. P. 1967
# Appendix: Sample Calculation
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 1
[]

[UserObjects]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 1.14e-3
  [../]
  [./eos]
    type = MoskitoEOS1P_IdealFluid
    reference_density = 1000
    specific_heat = 4190
    thermal_conductivity = 0.59
    thermal_expansion_0 = 0
    thermal_expansion_1 = 0
    bulk_modulus = 2e15
    reference_temperature = 273.15
    reference_enthalpy = 0
  [../]
  [./annulus]
    type = MoskitoAnnulus
    density = 0.6215
    viscosity = 0.00002853
    thermal_expansion = 1.755e-3
    thermal_conductivity = 0.0441
    heat_capacity = 1025.77
    emissivity_inner = 0.9
    emissivity_outer = 0.9
    convective_heat_method = Dropkin_Sommerscales
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '310.928'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = 1e5
    flowrate = 0
    well_direction = x
    well_diameter = 0.06096
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0 0 0'
    well_type = injection
  [../]
  [./Lateral]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.06096 0.089 0.2164 0.2438 0.3048'
    # high thermal conductivities used to ignor the effect of these layer to match the example
    conductivities = '1e5 0 1e5 0.3461'
    convective_thermal_resistance = false
    annulus_uo = annulus
    # Rock parameters
    formation_density = 2468.34
    formation_thermal_conductivity = 1.7307
    formation_heat_capacity = 950
    # Configuration of material
    formation_temperature_function = 310.928
    nondimensional_time_function = Hasan_Kabir_2012
    manual_time = 1814400
    output_properties = 'total_thermal_resistivity well_thermal_resistivity'
    outputs = exodus
  [../]
[]

[Variables]
  [./T]
    initial_condition = 588.706
  [../]
[]

[Kernels]
  [./Tkernel]
    type = NullKernel
    variable = T
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 50
  l_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
