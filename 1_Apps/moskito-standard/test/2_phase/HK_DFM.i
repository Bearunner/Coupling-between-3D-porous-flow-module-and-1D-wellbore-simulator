[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 200
  nx = 200
[]

[UserObjects]
  [./viscosity_gas]
    type = MoskitoViscosityConst
    viscosity = 1e-3
  [../]
  [./viscosity_liqid]
    type = MoskitoViscosityConst
    viscosity = 1e-3
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
    type = MoskitoFluidWell_2p1c
    pressure = p
    enthalpy = h
    massrate = m
    well_direction = x
    well_type = production
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    well_diameter = 0.05
    gravity = '9.8 0 0'
    outputs = exodus
    output_properties = 'well_velocity gas_velocity liquid_velocity void_fraction flow_type_c0 drift_velocity flow_pattern'
  [../]
[]

[Variables]
  [./h]
    initial_condition = 4.25e5
  [../]
  [./p]
    initial_condition = 1e5
  [../]
  [./m]
    [./InitialCondition]
      type = FunctionIC
      function = 'if(x<100,1e-2*x,(x-99)*2)*0.57'
    [../]
  [../]
[]

[Kernels]
  [./hkernel]
    type = NullKernel
    variable = h
  [../]
  [./pkernel]
    type = NullKernel
    variable = p
  [../]
  [./mkernel]
    type = NullKernel
    variable = m
  [../]
[]

[Executioner]
  type = Steady
  l_tol = 1e-10
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-9
  nl_max_its = 50
  solve_type = NEWTON
[]

[Outputs]
  exodus = true

[]
