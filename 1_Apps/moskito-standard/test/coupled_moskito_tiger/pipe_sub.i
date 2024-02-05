[Mesh]
  type = GeneratedMesh
  dim = 1
  xmax = 0.123
  nx = 10
[]

[UserObjects]
  [./eos]
    type = MoskitoEOSIdealFluid
    bulk_modulus = 2e+9
    reference_pressure = 0
    reference_temperature = 293.15
    reference_density = 998.29
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
[]

[Materials]
  [./wells]
    type = MoskitoFluidWell1P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.012
    manual_friction_factor = 0
    roughness_type = smooth
    gravity = '9.81 0 0'
    output_properties = 'density well_velocity'
    outputs = exodus
  [../]
[]

[BCs]
  [./pbc_p]
    type = PostprocessorDirichletBC
    variable = p
    boundary = right
    postprocessor = p_from_master
  [../]
  # [./pbc_p]
  #   type = DirichletBC
  #   variable = p
  #   boundary = right
  #   value = 0
  # [../]
  [./qbc_p]
    type = FunctionDirichletBC
    variable = q
    boundary = left
    function = flowrate
  [../]
[]

[Functions]
  [./flowrate]
    type = ParsedFunction
    vars = 'K A a L ho to'
    vals = '0.0117 7854 113 127 1200 122'
    value = '-K*A*ho/a/L*exp(-K*A/a/L*t)*a*1e-9'
  [../]
[]

[Variables]
  [./h]
  [../]
  [./p]
    scaling = 1e-2
  [../]
  [./q]
  [../]
[]

[Kernels]
  [./pkernelp]
    type = MoskitoMass
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernelp]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
  [../]
  [./hkernelp]
    type = NullKernel
    variable = h
  [../]
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
  # dt = 10
  # end_time = 2
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Postprocessors]
  [./p_from_master]
    type = Receiver
    # execute_on = NONLINEAR
  [../]
  [./q_right_sub]
    type = SideAverageValue
    boundary = right
    variable = q
    # execute_on = NONLINEAR
  [../]
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
