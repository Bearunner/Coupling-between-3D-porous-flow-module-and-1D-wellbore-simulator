[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0.123
  xmax = 0.25
  nx = 10
[]

[UserObjects]
  [./rock_uo]
    type =  TigerPermeabilityConst
    permeability_type = isotropic
    k0 = '1.1973e-12'
  [../]
[]

[Modules]
  [./FluidProperties]
    [./water_uo]
      type = TigerIdealWater
      viscosity = 0.001
      bulk_modulus = 2e+09
      reference_pressure = 0
      reference_temperature = 293.15
      reference_density = 998.29
    [../]
    # [./water_uo]
    #   type = TigerWaterConst
    #   viscosity = 0.001
    #   density = 1000
    # [../]
  [../]
[]

[Materials]
  [./matrix_g]
    type = TigerGeometryMaterial
    porosity = 1.0
    scale_factor = 7.854e-3
  [../]
  [./matrix_f]
    type = TigerFluidMaterial
    fp_uo = water_uo
    pressure = p
    temperature = 303.15
    output_properties = 'fluid_density'
    outputs = exodus
  [../]
  [./matrix_h]
    type = TigerHydraulicMaterialH
    pressure = p
    has_gravity = true
    # to make the gravity along +x dir
    gravity_acceleration = -9.81
    compressibility = 1.0e-9
    kf_uo = rock_uo
    output_properties = 'darcy_velocity'
    outputs = exodus
  [../]
[]

[BCs]
  [./pbc_s1]
    type = DirichletBC
    variable = p
    boundary = right
    value = 0
  [../]
  # another way of coupling with Neumann BC
  # [./qbc_s]
  #   type = PostprocessorNeumannBC
  #   variable = p
  #   boundary = left
  #   postprocessor = q_from_sub
  # [../]
  # boundary for constant head
  # [./qbc_s]
  #   type = NeumannBC
  #   variable = p
  #   boundary = left
  #   value = 2.2e-7
  # [../]
[]

[Variables]
  [./p]
    scaling = 1e6
  [../]
[]

[DiracKernels]
  [./inj]
    type = TigerHydraulicPointSourceH
    point = '0.123 0 0'
    variable = p
    mass_flux_function = massrate
  [../]
[]

[Functions]
  [./massrate]
    type = ParsedFunction
    vars = 'q'
    vals = 'q_from_sub'
    value = 'q*996.15'
  [../]
[]

[Kernels]
  [./pkernels]
    type = TigerHydraulicKernelH
    variable = p
  [../]
  # [./pkernels_t]
  #   type = TigerHydraulicTimeKernelH
  #   variable = p
  # [../]
[]

[Preconditioning]
  [./p2]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
[]

[Executioner]
  type = Transient
  dt = 10
  end_time = 122
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
  picard_max_its = 50
  picard_abs_tol = 1e-12
  picard_rel_tol = 1e-05
[]

[Postprocessors]
  [./p_left_master]
    type = SideAverageValue
    boundary = left
    variable = p
    # execute_on = NONLINEAR
  [../]
  [./q_from_sub]
    type = Receiver
    # execute_on = NONLINEAR
  [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]

[MultiApps]
  [./pipe_sand]
    type = TransientMultiApp
    input_files = pipe_sub.i
    positions = '0 0 0'
    # relative to place of execution
    library_path = ../../lib
    app_type = MoskitoApp
    execute_on = TIMESTEP_BEGIN
  [../]
[]

[Transfers]
  [./pressure_transfer]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    multi_app = pipe_sand
    from_postprocessor = p_left_master
    to_postprocessor = p_from_master
    execute_on = SAME_AS_MULTIAPP
  [../]
  [./flowrate_transfer]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = pipe_sand
    from_postprocessor = q_right_sub
    to_postprocessor = q_from_sub
    reduction_type = average
    execute_on = SAME_AS_MULTIAPP
  [../]
[]
