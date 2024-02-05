[Mesh]
  type = FileMesh
  file = well.msh
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 0.0001
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 0.001
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

[GlobalParams]
  pressure = p
  enthalpy = h
  flowrate = q
  well_direction = -y
  eos_uo = eos
  viscosity_uo = viscosity_2p
  drift_flux_uo = df
  roughness_type = smooth
  gravity = '0 -9.8 0'
  outputs = exodus
  output_properties = 'gas_velocity liquid_velocity void_fraction mass_fraction flow_pattern current_phase gas_density liquid_density density temperature well_velocity'
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell2P
    well_diameter = 0.45
    block = 3
  [../]
  [./area1]
    type = MoskitoFluidWell2P
    well_diameter = 0.3
    block = 4
  [../]
  [./area2]
    type = MoskitoFluidWell2P
    well_diameter = 0.15
    block = 5
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = top
    value = 1e5
  [../]
  [./qbc]
    type = DirichletBC
    variable = q
    boundary = top
    value = -0.015
  [../]
[]

[Variables]
  [./h]
    [./InitialCondition]
      type = FunctionIC
      variable = h
      # function = 8e5-y*125
      function = 1.5e6-25*y
    [../]
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      variable = p
      function = 1e5
    [../]
  [../]
  [./q]
    scaling = 1e-6
    initial_condition = -0.015
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
  # [./pkernel1]
  #   type = NullKernel
  #   variable = p
  # [../]
  # [./qkernel1]
  #   type = NullKernel
  #   variable = q
  # [../]
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
  # l_tol = 1e-13
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  # [./test]
  #   type = VariableResidualNormsDebugOutput
  # [../]
[]
