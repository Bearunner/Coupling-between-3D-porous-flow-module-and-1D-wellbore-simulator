[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 600
  nx = 300
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
    type = MoskitoDFShi
    surface_tension = 0.0288
    Pan_param_cMax = 1.2
    Shi_param_Fv = 0.3
  [../]
  [./eos]
    type = MoskitoPureWater2P
  [../]
[]

[Materials]
  [./area]
    type = MoskitoFluidWell_2p1c
    well_diameter = 0.1
    pressure = p
    enthalpy = h
    massrate = m
    well_direction = x
    well_type = production
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    roughness_type = smooth
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = 'well_velocity void_fraction density temperature'
  [../]
[]

[BCs]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = 1e5
  [../]
  [./mbc]
    type = FunctionDirichletBC
    variable = m
    boundary = right
    function = 'if(t<15, 0.06*t, 0.9)'
  [../]
  [./hbc]
    type = FunctionDirichletBC
    variable = h
    boundary = right
    function = 'if(t<20, 1e5, if(t>40, 5.5e5, 1e5+2.25e4*(t-20)))'
  [../]
[]

[Variables]
  [./p]
    initial_condition = 1e5
  [../]
  [./m]
    initial_condition = 1e-5
  [../]
  [./h]
    initial_condition = 5.5e5
  [../]
[]

[Kernels]
  [./mkernel]
    type = MoskitoMass_2p1c
    variable = m
  [../]
  [./pkernel]
    type = MoskitoMomentum_2p1c
    variable = p
    massrate = m
    enthalpy = h
  [../]
  [./hkernel]
    type = MoskitoEnergy_2p1c
    variable = h
    massrate = m
    pressure = p
  [../]
[]

[Preconditioning]
  active = pn1
  [./p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                 '
  [../]
  [./pn1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type  -snes_type -snes_linesearch_type'
    petsc_options_value = ' bjacobi  ilu                           newtonls   basic               '
  [../]
[]

[Dampers]
  [./p]
    type = BoundingValueNodalDamper
    variable = p
    min_value = 1e4
    max_value = 6e6
    min_damping = 1e-10
  [../]
  [./h]
    type = BoundingValueNodalDamper
    variable = h
    min_value = 1e4
    max_value = 2e6
    min_damping = 1e-10
  [../]
[]

[Executioner]
  type = Transient
  dt = 5
  end_time = 50
  l_max_its = 50
  nl_max_its = 100
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-7
  solve_type = NEWTON
  automatic_scaling = true
  compute_scaling_once = false
  [./Quadrature]
    type = GAUSS
    element_order = FIRST
  [../]
[]


[Outputs]
  # execute_on = 'INITIAL NONLINEAR TIMESTEP_END'
  print_linear_residuals = false
  exodus = true
[]
