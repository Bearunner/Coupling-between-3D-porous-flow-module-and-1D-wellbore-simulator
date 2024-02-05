[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 3000
  nx = 200
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
    type = MoskitoMixture2P
  [../]
[]

[Functions]
  [./grad_func]
    type = ParsedFunction
    value = '293.15 + 0.03 * x'
  [../]
  [./dts]
    type = PiecewiseLinear
    x = '0 3153600 31536000'
    y = '3600 86400 864000'
  [../]
[]

[Materials]
  [./area]
    type = MoskitoFluidWell_2p1c
    well_diameter = 0.2
    pressure = p
    enthalpy = h
    massrate = m
    concentration = c
    c2 = c2
    well_direction = x
    well_type = production
    eos_uo = eos
    viscosity_uo = viscosity_2p
    drift_flux_uo = df
    roughness_type = smooth
    gravity = '9.8 0 0'
  [../]
  [./Lateral]
    type = MoskitoLatHeat_Inc_Formation_2p
     outer_diameters = '0.2 0.24 0.3'
     conductivities = '50.0 0.7'
     formation_density = 2650
     formation_thermal_conductivity = 2.5
     formation_heat_capacity = 1000
     # Configuration of material
     formation_temperature_function = grad_func
     nondimensional_time_function = Ramey_1981_BF
     #output_properties = 'total_thermal_resistivity'
     #outputs = exodus
   [../]
   [./sumgrad]
     type = MoskitoGrad
      concentration = c
      c2 = c2
      #output_properties = 'drho_dc'
      #outputs = exodus
    [../]
[]

[BCs]
  [./pbc]
    type = FunctionDirichletBC
    variable = p
    boundary = left
    function =  '2e7'
  [../]
  [./mbc]
    type = FunctionDirichletBC
    variable = m
    boundary = right
    function = '40'
  [../]
  [./hbc]
    type = FunctionDirichletBC
    variable = h
    boundary = right
    function = '352634'
  [../]
  [./testibc]
    type = FunctionDirichletBC
    variable = c
    boundary = right
    function = '0.3'
  [../]
  [./testibc2]
    type = FunctionDirichletBC
    variable = c2
    boundary = right
    function = '0.01'
  [../]
[]




[Variables]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      function = '20000000+1000*9.8*x'
    [../]
  [../]
  [./m]
    initial_condition = 40
  [../]
  [./h]
    [./InitialCondition]
      type = FunctionIC
      function = '168050+x*(352634-168050)/3000'
    [../]
  [../]
  [./c]
    initial_condition = 0.3
  [../]
  [./c2]
    initial_condition = 0.01
  [../]
[]



[Kernels]
  [./mkernel]
    type = MoskitoMass_2p1c
    variable = m
  [../]
  [./mkernelt]
    type = MoskitoTimeMass_2p1c
    variable = m
    enthalpy = h
    pressure = p
    molarity = c
  [../]
  [./pkernel]
    type = MoskitoMomentum_2p1c
    variable = p
    enthalpy = h
    massrate = m
    molarity = c
  [../]
  [./pkernelt]
    type = MoskitoTimeMomentum_2p1c
    variable = p
    massrate = m
  [../]
  [./hkernel]
    type = MoskitoEnergy_2p1c
    variable = h
    massrate = m
    pressure = p
    molarity = c
  [../]
  [./do]
    type = MoskitoTimeEnergy_2p1c
    variable = h
    massrate = m
    pressure = p
    molarity = c
  [../]
  [./heat]
    type = MoskitoLatHeatIncFormation_2p
    variable = h
    pressure = p
    molarity = c
  [../]

[]

[convectiondiffusion]
    variables = 'p h m c c2'
    index_CO2 = 4
    index_presure = 1
    index_enthalpy = 2
    index_massrate = 3
    index_CH4 = 5
    #index_CH4 = 5
[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 315360000
  nl_max_its = 15
  solve_type = NEWTON
   dtmin = 2e-14
  [./TimeStepper]
    type =  FunctionDT
    function = dts
  [../]
[]


[Outputs]
  print_linear_residuals = false
  exodus = true
[]
