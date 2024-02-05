
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1500
  nx = 100
[]

# [MeshModifiers]
#   [./rotate]
#     type = Transform
#     transform = ROTATE
#     vector_value = '3 0 0'
#   [../]
# []

[UserObjects]
  [./viscosity]
    type = MoskitoViscositySmith
  [../]
  [./eos]
    type = MoskitoEOSIdealFluid
    thermal_expansion_0 = 4e-04
    thermal_expansion_1 = 7e-7
    reference_density = 798.92
    reference_enthalpy = 1085800
    reference_pressure = 4.0e6
    reference_temperature = 523.15
    specific_heat = 4800
  [../]
  # [./df]
  #   type = MoskitoDFHK
  #   surface_tension = 0.07
  # [../]
  # [./viscosity_gas]
  #   type = MoskitoViscosityConst
  #   viscosity = 0.00001
  # [../]
  # [./viscosity_liquid]
  #   type = MoskitoViscosityConst
  #   viscosity = 0.001
  # [../]
  # [./viscosity]
  #   type = MoskitoViscosity2P
  #    ve_uo_gas = viscosity_gas
  #    ve_uo_liquid = viscosity_liquid
  # [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '0.03513 * x'
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    well_diameter = 0.076
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = 'temperature density well_velocity specific_heat well_reynolds_no well_moody_friction viscosity diameter'
  [../]
  [./Lateral]
    type = MoskitoLatHeatIterationXiongSI
     radius_wellbore = 0.124
     radius_tubbing_outer = 0.04445
     radius_casing_inner = 0.081
     radius_cement = 0.089
     conductivity_cement = 0.35
     conductivity_tubing = 43.2639
     conductivity_casing = 43.2639
     # Rock parameters
     Surface_temperature = 285.15
     thermal_diffusivity_rock = 0.738063e-6
     conductivity_rock = 1.73
     # Annulus parameters
     density_annulus = 800
     conductivity_annulus = 0.58
     dyn_viscosity_annulus = 0.0001055
     capacity_annulus = 4870
     thermal_expansion_annulus = 0.001494
     emissivity_annulus_outer = 0.8
     emissivity_annulus_inner = 0.8
     # Configuration of material
     geothermal_gradient = grad_func
     hc_calucation_model = Dropkin_Sommerscales
     time_model = simulation_time
     # internal_solve_full_iteration_history = true
     output_properties = 'Temperature_Rankine temperature_well_formation_interface heat_loss formation_temperature'
     outputs = exodus
   [../]
[]

[Variables]
  [./h]
    scaling = 1e-6
    # initial_condition = 50350
    [./InitialCondition]
      type = FunctionIC
      variable = h
      function = ' 135 *x + 85350'
    [../]
  [../]
  [./p]
    initial_condition = 4.1e6
    scaling = 1e1
  [../]
  [./q]
    initial_condition = 0.004050926
    scaling = 1e-2
  [../]
[]

# [ICs]
#   [./hic]
#     type = FunctionIC
#     variable = h
#     function = grad_func
#   [../]
# []

[BCs]
  # [./hbc]
  #   type = DirichletBC
  #   variable = h
  #   boundary = left
  #   value = 83950
  # [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = right
    value = 0.004050926
  [../]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 4.1e6
  [../]
  [./hbc]
    type = MoskitoTemperatureToEnthalpy1P
    variable = h
    pressure = p
    boundary = left
    temperature = 523.45
    eos_uo = eos
  [../]
[../]

[Kernels]
  [./hkernel]
    type = MoskitoEnergy
    variable = h
    pressure = p
    flowrate = q
  [../]
  [./h_time]
    type = MoskitoTimeEnergy
    variable = h
    flowrate = q
    pressure = p
  [../]
  [./qkernel1]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./p_time]
    type = MoskitoTimeMass
    variable = p
    enthalpy = h
  [../]
  [./pkernel1]
    type = MoskitoMass
    variable = p
    enthalpy = h
    flowrate = q
  [../]
  [./q_time]
    type = MoskitoTimeMomentum
    variable = q
    enthalpy = h
    pressure = p
  [../]
  [./heat]
    type = MoskitoHeat
    variable = h
  [../]
[]

[Preconditioning]
  active = p3
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
  [./p2]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -sub_pc_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu newtonls basic NONZERO 51'
  [../]
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
  l_max_its = 50
  l_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
  solve_type = NEWTON
  # end_time = 2592000
  end_time = 31536000
  nl_abs_tol = 1e-7
  dtmax = 30000
 [./TimeStepper]
   type = SolutionTimeAdaptiveDT
   dt = 10
   percent_change = 0.2
 [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  # [./test]
  #   type = VariableResidualNormsDebugOutput
  # [../]
  console = true
[]
