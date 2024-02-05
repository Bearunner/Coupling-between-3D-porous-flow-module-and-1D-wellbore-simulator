[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 300
  nx = 100
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
    output_properties = 'drift_velocity well_velocity void_fraction density temperature gas_velocity liquid_velocity '
  [../]
[]

[BCs]
  [./pbc]
    type = FunctionDirichletBC
    variable = p
    boundary = left
    function = 'if(t<=10,2e5-1e4*t,1e5)'
  [../]
  [./mbc]
    type = FunctionDirichletBC
    variable = m
    boundary = left
    function = 'if(t<10,0,0.15*(t-10))'
  [../]
[]

[Variables]
  [./p]
    initial_condition = 2e5
  [../]
  [./m]
    initial_condition = 0.0
  [../]
[]

[AuxVariables]
  [./h]
    order = FIRST
    family = LAGRANGE
    initial_condition = 5e5
  [../]
[]

[Kernels]
  [./pkernel]
    type = MoskitoMomentum_2p1c
    variable = p
    massrate = m
    enthalpy = h
  [../]
  [./mkernel]
    type = MoskitoMass_2p1c
    variable = m
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
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -snes_type -snes_linesearch_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                   newtonls   basic               '
  [../]
[]

[Dampers]
  [./p]
    type = BoundingValueNodalDamper
    variable = p
    min_value = 1e4
    max_value = 5e6
    min_damping = 1e-5
  [../]
[]

[Executioner]
  type = Transient
  dt = 5
  end_time = 20
  l_max_its = 50
  nl_max_its = 50
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  solve_type = NEWTON
  automatic_scaling = true
  compute_scaling_once = false
  [./Quadrature]
    type = GAUSS
    element_order = FIRST
  [../]
[]


[Outputs]
  exodus = true
[]
