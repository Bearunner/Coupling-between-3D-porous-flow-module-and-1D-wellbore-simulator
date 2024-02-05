[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 901
  nx = 901
[]

[MeshModifiers]
  [./rotate]
    type = Transform
    transform = ROTATE
    vector_value = '2.87 0 0'
  [../]
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 1.3e-5
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 1.44e-4
  [../]
  [./viscosity_2p]
    type = MoskitoViscosity2P
    ve_uo_gas = viscosity_gas
    ve_uo_liquid = viscosity_liqid
  [../]
  [./df]
    type = MoskitoDFHK
    surface_tension = 0.0288
  [../]
  [./eos]
    type = MoskitoPureWater2P
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    pressure = p
    enthalpy = h
    flowrate = q
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    well_diameter = 0.224
    roughness_type = smooth
    gravity = '3 0 0'
    outputs = exodus
    output_properties = 'profile_mixture_density gas_velocity liquid_velocity void_fraction flow_pattern current_phase gas_density liquid_density mass_fraction density specific_heat temperature'
  [../]
[]

[BCs]
  [./pbc_1]
    type = DirichletBC
    variable = p
    boundary = right
    value = 3.292e6
  [../]
  [./qbc_1]
    type = MoskitoMassFlowRate
    variable = q
    boundary = left
    mass_flowrate = -24.95
    mixture_density = 817.42
  [../]
  # [./hbc]
  #   type = DirichletBC
  #   variable = h
  #   boundary = right
  #   value = 0.995e6
  # [../]
  # [./pbc_1]
  #   type = NeumannBC
  #   variable = p
  #   boundary = right
  #   value = -24.95
  # [../]
  # [./qbc_2]
  #   type = NeumannBC
  #   variable = q
  #   boundary = right
  #   value = -3.292e6
  # [../]
[]

# [DiracKernels]
#   [./flow]
#     type = ConstantPointSource
#     point = '0 0 0'
#     variable = p
#     value = -0.000001
#   [../]
# []

[Variables]
  [./h]
    initial_condition = 0.995e6
    # [./InitialCondition]
    #   type = FunctionIC
    #   variable = h
    #   function = 630000+411*x
    # [../]
  [../]
  [./p]
    initial_condition = 1e5
  [../]
  [./q]
    scaling = 1e-6
    # initial_condition = -0.046296296
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
  [./pkernel]
    type = MoskitoMass
    variable = p
    flowrate = q
    enthalpy = h
  [../]
  [./qkernel]
    type = MoskitoMomentum
    variable = q
    pressure = p
    enthalpy = h
  [../]
[]

[Preconditioning]
  active = p2
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
  [./p4]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
[]

[Executioner]
  type = Steady
  l_tol = 1e-12
  l_max_its = 50
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
