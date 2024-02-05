
[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 100
  nx = 50
[]

[Modules]
  [./FluidProperties]
    [./water_uo]
      type = Water97FluidProperties
      execute_on = TIMESTEP_BEGIN
    [../]
  [../]
[]

[UserObjects]
  [./eos]
    type = MoskitoEOS1P_FPModule
    SinglePhase_fp = water_uo
  [../]
  [./viscosity]
    type = MoskitoViscosityWaterSmith
  [../]
[]

[Materials]
  [./area0]
    type = MoskitoFluidWell_1p1c
    pressure = p
    temperature = T
    flowrate = q
    well_type = production
    well_direction = x
    well_diameter = 0.0890016
    eos_uo = eos
    viscosity_uo = viscosity
    roughness_type = smooth
    gravity = '10.0 0 0'
    output_properties = 'density well_velocity'
    outputs = exodus
  [../]
[]

[Variables]
  [./T]
    initial_condition = 320
  [../]
  [./p]
    [./InitialCondition]
      type = FunctionIC
      function = '2e6'
    [../]
  [../]
  [./q]
  [../]
[]

[BCs]
  [./Tbc_top]
    type = DirichletBC
    variable = T
    boundary = right
    value = 400
  [../]
  [./pbc]
    type = DirichletBC
    variable = p
    boundary = left
    value = '2e6'
  [../]
  [./qbc]
    type = FunctionDirichletBC
    variable = q
    boundary = left
    function = 'if(t>20,0.01,0.00001)'
  [../]
[]

[Kernels]
  [./Tkernel]
    type = MoskitoEnergy_1p1c
    variable = T
    pressure = p
    flowrate = q
    gravity_energy = true
  [../]
  [./Ttimekernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
    pressure = p
    flowrate = q
  [../]
  [./pkernel]
    type = MoskitoMass_1p1c
    variable = p
    flowrate = q
    temperature = T
  [../]
  [./ptimekernel]
    type = MoskitoTimeMass_1p1c
    variable = p
    temperature = T
  [../]
  [./qkernel]
    type = MoskitoMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
  [../]
  [./qtimekernel]
    type = MoskitoTimeMomentum_1p1c
    variable = q
    pressure = p
    temperature = T
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

[Executioner]
  type = Transient
  dt = 1
  end_time = 120
  l_max_its = 50
  nl_max_its = 50
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  solve_type = NEWTON
  automatic_scaling = true
[]

[Outputs]
  exodus = true
[]
