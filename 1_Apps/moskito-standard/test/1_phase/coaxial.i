[Mesh]
  type = FileMesh
  file = coaxial.msh
  allow_renumbering = false
[]

[UserObjects]
  [./eos]
    type = MoskitoEOS1P_IdealFluid
    bulk_modulus = 2e+012
    reference_temperature = 275
    reference_enthalpy = 0
    reference_density = 1000
  [../]
  [./viscosity]
    type = MoskitoViscosityConst
    viscosity = 0.001
  [../]
[]

[Materials]
  [./Ti_mat]
    type = MoskitoFluidWell_1p1c
    pressure = Pi
    temperature = Ti
    flowrate = Qi
    well_direction = x
    well_type = injection
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.1
    roughness_type = smooth
    gravity = '9.8 0 0'
    block = inner
    outputs = exodus
    output_properties = 'well_reynolds_no'
  [../]
  [./To1_mat]
    type = MoskitoFluidWell_1p1c
    pressure = Po
    temperature = To
    flowrate = Qo
    well_direction = x
    well_type = injection
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.2
    manual_cross_section_area = 0.0235619
    manual_wetted_perimeter = 0.9424778
    roughness_type = smooth
    gravity = '9.8 0 0'
    block = outer
    outputs = exodus
    output_properties = 'well_reynolds_no'
  [../]
  [./To2_mat]
    type = MoskitoFluidWell_1p1c
    pressure = Po
    temperature = To
    flowrate = Qo
    well_direction = x
    well_type = injection
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.15
    roughness_type = smooth
    gravity = '9.8 0 0'
    block = well
  [../]
  [./coaxial]
    type = MoskitoCoaxialHeat_1p
    eos_uo = eos
    viscosity_uo = viscosity
    pressure_inner_pipe = Pi
    pressure_outer_pipe = Po
    temperature_inner_pipe = Ti
    temperature_outer_pipe = To
    flowrate_inner_pipe = Qi
    flowrate_outer_pipe = Qo
    inner_pipe_outer_radius = 0.05
    inner_pipe_wall_thickness = 0.01
    outer_pipe_inner_radius = 0.1
    conductivity_inner_pipe = 80.0
    block = 'inner outer'
    outputs = exodus
    output_properties = 'overall_heat_transfer_coeff Nusselt_number_inner_pipe Nusselt_number_outer_pipe'
  [../]
[]

[BCs]
  [./Po_bc]
    type = DirichletBC
    variable = Po
    boundary = left
    value = 2e6
  [../]
  [./To_bc]
    type = FunctionDirichletBC
    variable = To
    boundary = left
    function = 'if(t<1000, 283.15, if(t>2000, 313.15, 283.15))'
  [../]
  [./Qo_bc]
    type = DirichletBC
    variable = Qo
    boundary = left
    value = 0.02
  [../]
[]

[Variables]
  [./Ti]
    [./InitialCondition]
      type = FunctionIC
      function = '293.15+0.03*x'
    [../]
    scaling = 1e-2
    block = 'outer inner'
  [../]
  [./To]
    [./InitialCondition]
      type = FunctionIC
      function = '293.15+0.03*x'
    [../]
    scaling = 1e-2
  [../]
  [./Po]
    [./InitialCondition]
      type = FunctionIC
      function = '2e6+9.8e3*x'
    [../]
    scaling = 1e5
  [../]
  [./Qo]
    initial_condition = 1e-5
    scaling = 1e-1
  [../]
[]

[AuxVariables]
  [./Qi]
    initial_condition = 1e-10
    block = 'outer inner'
  [../]
  [./Pi]
    [./InitialCondition]
      type = FunctionIC
      function = '2e6+9.8e3*x'
    [../]
    block = 'outer inner'
  [../]
[]

[Kernels]
  [./Ti_kernel]
    type = MoskitoEnergy_1p1c
    variable = Ti
    pressure = Pi
    flowrate = Qi
    block = inner
  [../]
  [./Ti_kernel1]
    type = MoskitoTimeEnergy_1p1c
    variable = Ti
    pressure = Pi
    flowrate = Qi
    block = inner
  [../]
  [./Ti_kernel2]
    type = MoskitoLatHeatCoaxial_1p
    variable = Ti
    Tsur = To
    # factor = 1.24
    block = inner
  [../]
  [./To_kernel2]
    type = MoskitoLatHeatCoaxial_1p
    variable = To
    Tsur = Ti
    # factor = 1.24
    block = outer
  [../]
  [./To_kernel]
    type = MoskitoEnergy_1p1c
    variable = To
    pressure = Po
    flowrate = Qo
    block = 'outer well'
  [../]
  [./To_kernel1]
    type = MoskitoTimeEnergy_1p1c
    variable = To
    pressure = Po
    flowrate = Qo
    block = 'outer well'
  [../]
  [./Po_kernel]
    type = MoskitoMass_1p1c
    variable = Po
    flowrate = Qo
    temperature = To
    block = 'outer well'
  [../]
  [./Qo_kernel]
    type = MoskitoMomentum_1p1c
    variable = Qo
    pressure = Po
    temperature = To
    block = 'outer well'
  [../]
  [./Po_kernel1]
    type = MoskitoTimeMass_1p1c
    variable = Po
    temperature = To
    block = 'outer well'
  [../]
  [./Qo_kernel1]
    type = MoskitoTimeMomentum_1p1c
    variable = Qo
    pressure = Po
    temperature = To
    block = 'outer well'
  [../]
[]

[Preconditioning]
  [./pn1]
    type = SMP
    full = true
    coupled_groups = 'To,Qo,Po Ti'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -snes_type -snes_linesearch_type'
    petsc_options_value = ' bjacobi  ilu          NONZERO                   newtonls   basic               '
  [../]
[]

[Executioner]
  type = Transient
  dt = 100
  end_time = 3000
  l_tol = 1e-8
  l_max_its = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_max_its = 50
  solve_type = NEWTON
  # automatic_scaling = true
[]

[Postprocessors]
  [./node]
    type = NodalVariableValue
    nodeid = 3
    variable = To
    outputs = csv
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
