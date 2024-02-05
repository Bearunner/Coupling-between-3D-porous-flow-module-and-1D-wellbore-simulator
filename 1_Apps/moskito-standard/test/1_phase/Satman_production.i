# Validation of Heat exchange using the paper from Satman & Tureyen 2016
# Actually ramey analytical solution
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 2500
  nx = 400
[]

[UserObjects]
  [./viscosity]
    type = MoskitoViscosityWaterSmith
  [../]
  [./eos]
    type = MoskitoEOS1P_IdealFluid
    specific_heat = 3160
    reference_density = 1000
    reference_temperature = 273.15
    reference_enthalpy = 0
    thermal_expansion_0 = 0
    thermal_expansion_1 = 0
    bulk_modulus = 2e15
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '293.15 + 0.09 * x'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell_1p1c
    temperature = T
    pressure = p
    flowrate = q
    well_direction = x
    well_diameter = 0.3
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '0.0 0 0'
    well_type = production
  [../]
  [./Lateral]
    type = MoskitoLatHeat_Inc_Formation_1p
     temperature_inner = T
     outer_diameters = '0.3 0.30000000001'
     conductivities = '2.92'
     convective_thermal_resistance = true
     # Rock parameters
     formation_density = 2650
     formation_thermal_conductivity = 2.92
     formation_heat_capacity = 1000
     # Configuration of material
     formation_temperature_function = grad_func
     nondimensional_time_function = Ramey_1981_BF
     output_properties = 'total_thermal_resistivity'
     outputs = exodus
   [../]
[]

[Variables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = grad_func
    [../]
  [../]
[]

[AuxVariables]
  [./p]
    initial_condition = 101325
  [../]
  [./q]
    initial_condition = 0.02
  [../]
[]

[BCs]
  [./Tbc]
    type = DirichletBC
    variable = T
    boundary = right
    value =  518.15
  [../]
[]

[Kernels]
  [./Tkernel]
    type = MoskitoEnergy_1p1c
    variable = T
    pressure = p
    flowrate = q
  [../]
  [./Ttkernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
    pressure = p
    flowrate = q
  [../]
  [./heat]
    type = MoskitoLatHeatIncFormation_1p
    variable = T
  [../]
[]

[Preconditioning]
  active = p3
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
  [../]
[]

[Executioner]
  type = Transient
  end_time = 2592000
  l_max_its = 50
  l_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
  solve_type = NEWTON
  nl_abs_tol = 1e-7
  [./TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '0 8640 86400 864000 2592000'
  [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
[]
