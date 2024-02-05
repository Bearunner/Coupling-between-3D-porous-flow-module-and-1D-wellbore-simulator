5_day = 129600   ##5 day
starting_time = 0
1st_step = 1
end_time = 9.461e+8  ## 30 years

WI_mw = 4.78E-14
WI_sw = 4.38E-14

k_reservoir = 3.45e-13
n_reservoir = 0.1

k_matrix = 2.6e-17
n_matrix = 0.05

lambda_solid = 2.67
lambda_fluid = 0.6

[Mesh]
  [./total]
 type = GeneratedMeshGenerator
 dim = 3
 nx = 24
 ny = 31
 nz = 31

 xmin = -510
 xmax = -390

 ymin = -77.5
 ymax = 77.5

 zmin = -77.5
 zmax = 77.5
[../]

[reservoir]
  type = ParsedSubdomainMeshGenerator
  input = total
  combinatorial_geometry = 'x >= -500 & x <= -400'
  block_id = 1
  block_name = reservoir
[]

[caprock]
  type = ParsedSubdomainMeshGenerator
  input = reservoir
  combinatorial_geometry = 'x >= -400 & x <= -300'
  block_id = 2
  block_name = caprock
[]

[basement]
  type = ParsedSubdomainMeshGenerator
  input = caprock
  combinatorial_geometry = 'x >= -600 & x <= -500'
  block_id = 3
  block_name = basement
[]

[./mw]
 type = ParsedSubdomainMeshGenerator
 combinatorial_geometry = 'x >= -500 & x <= -400 & y >= -2.5 & y <= 2.5 & z >= -2.5 & z <= 2.5'
 block_id  = '10'
 block_name = mw
 input = 'basement'
 excluded_subdomains = 'caprock basement'
[../]

[./sw1]
 type = ParsedSubdomainMeshGenerator
 combinatorial_geometry = 'x >= -500 & x <= -400 & y >= 47.5 & y <= 52.5 & z >= -2.5 & z <= 2.5'
 block_id  = '11'
 input = 'mw'
 excluded_subdomains = 'caprock basement'
[../]

[./sw2]
 type = ParsedSubdomainMeshGenerator
 combinatorial_geometry = 'x >= -500 & x <= -400 & y >= -52.5 & y <= -47.5 & z >= -2.5 & z <= 2.5'
 block_id  = '12'
 input = 'sw1'
 excluded_subdomains = 'caprock basement'
[../]

[./sw3]
 type = ParsedSubdomainMeshGenerator
 combinatorial_geometry = 'x >= -500 & x <= -400 & z >= 47.5 & z <= 52.5 & y >= -2.5 & y <= 2.5'
 block_id  = '13'
 input = 'sw2'
 excluded_subdomains = 'caprock basement'
[../]

[./sw4]
 type = ParsedSubdomainMeshGenerator
 combinatorial_geometry = 'x >= -500 & x <= -400 & z >= -52.5 & z <= -47.5 & y >= -2.5 & y <= 2.5'
 block_id  = '14'
 input = 'sw3'
 excluded_subdomains = 'caprock basement'
[../]

[rename]
 type = RenameBlockGenerator
 old_block = '11 12 13 14'
 new_block = 'sw sw sw sw'
 input = 'sw4'
[]
[]

[GlobalParams]
  PorousFlowDictator = dictator
  multiply_by_density = true
  porepressure = porepressure
  temperature = temperature
  gravity = '-9.81 0 0'
  pp =  porepressure
  # execute_on = 'initial timestep_begin'
[]

[FluidProperties]
  [water_uo]
    type = SimpleFluidProperties_new
    bulk_modulus = 2E9
    density0 = 1000.0
    thermal_expansion = 0.0
    cv = 4180.0
    cp = 4180.0
  []
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure temperature' ## do not lose anyone!!!
    number_fluid_phases = 1
    number_fluid_components = 1 ## component 0 is C1,component 1 is water. Water is always the last component!!!
  [../]
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
  []
  [./simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = water_uo
    phase = 0
  [../]
  [diff_aquifer]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '0'
    tortuosity = 0.1
    block = ' reservoir mw sw'
  []
  [diff_matrix]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '0'
    tortuosity = 0.05
    block = 'caprock basement'
  []
  [porosity_aquifer]
    type = PorousFlowPorosityConst
    porosity =  '${n_reservoir}'
    block = ' reservoir mw sw'
  []
  [porosity_matrix]
    type = PorousFlowPorosityConst
    porosity =  '${n_matrix}'
    block = 'caprock basement'
  []
  [permeability_aquifer]
    type = PorousFlowPermeabilityConst
    permeability = '${k_reservoir} 0 0   0 ${k_reservoir} 0   0 0 ${k_reservoir}'
    block = ' reservoir mw sw'
  []
  [permeability_matrix]
    type = PorousFlowPermeabilityConst
    permeability = '${k_matrix} 0 0   0 ${k_matrix} 0   0 0 ${k_matrix}'
    block = 'caprock basement'
  []
  [relp]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 750
    density = 2535
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityFromPorosity
    lambda_s = '${lambda_solid} 0 0   0 ${lambda_solid} 0   0 0 ${lambda_solid}'
    lambda_f = '${lambda_fluid} 0 0   0 ${lambda_fluid} 0   0 0 ${lambda_fluid}'
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
[]

[BCs]
[]

[ICs]
  [./p_ic]
   type = FunctionIC
   variable = porepressure
   function = hydrostatic
  [../]
  [./T_ic]
    type = FunctionIC
    variable = temperature
    function = thermalgradient
   [../]
[]

[Functions]
  [./hydrostatic]
    type = ParsedFunction
    expression = '3e5-1000*9.81*x'
  [../]
  [./thermalgradient]
    type = ParsedFunction
    expression = '283.15-0.03*x'
  [../]
[]

[Variables]
  [porepressure]
    scaling = 1e-1
  []
  [temperature]
    scaling = 1e-7
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = porepressure
  []
  [adv0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = porepressure
  []
  [diff0]
    type = PorousFlowDispersiveFlux
    fluid_component = 0
    variable = porepressure
    disp_trans = 0
    disp_long = 0
  []
  [EnergyTransient]
    type = PorousFlowEnergyTimeDerivative
    variable = temperature
  []
  [EnergyConduciton]
    type = PorousFlowHeatConduction
    variable = temperature
  []
  [EnergyAdvection]
    type = PorousFlowHeatAdvection
    variable = temperature
  []

  [./source_mass_1st_mw]
    type = PorousFlowPeacemanBoreholeMass
    variable = porepressure
    WI = ${WI_mw}
    p_well = 1
    function = '25e5-1000*9.81*x'
    rho_well = 1000
    mu_well = 1e-3
    block = 'mw'
  [../]

  [./source_mass_mw]
    type = PorousFlowPeacemanBoreholeMass
    variable = porepressure
    WI = ${WI_mw}
    p_well = pp_mw_well
    rho_well = rho_mw_well
    mu_well = mu_mw_well
    block = 'mw'
  [../]

  [./source_energy_1st_mw]
    type = PorousFlowPeacemanBoreholeEnergy
    variable = temperature
    WI = ${WI_mw}
    p_well = 1
    function = '25e5-1000*9.81*x'
    rho_well = 1000
    mu_well = 1e-3
    h_well = 1280000
    block = 'mw'
  [../]

  [./source_energy_mw]
    type = PorousFlowPeacemanBoreholeEnergy
    variable = temperature
    WI = ${WI_mw}
    p_well = pp_mw_well
    rho_well = rho_mw_well
    mu_well = mu_mw_well
    h_well = h_mw_well
    block = 'mw'
  [../]

  [./sink_mass_1st_sw]
    type = PorousFlowPeacemanBoreholeMassPro
    variable = porepressure
    WI = ${WI_sw}
    p_well = 1
    function = '1e5-1000*9.81*x'
    block = 'sw'
  [../]

  [./sink_mass_sw]
    type = PorousFlowPeacemanBoreholeMassPro
    variable = porepressure
    WI = ${WI_sw}
    p_well = pp_sw_well
    block = 'sw'
  [../]

  [./sink_energy_1st_sw]
      type = PorousFlowPeacemanBoreholeEnergyPro
      variable = temperature
      WI = ${WI_sw}
      p_well = 1
      function = '1e5-1000*9.81*x'
      block = 'sw'
  [../]

  [./sink_energy_sw]
    type = PorousFlowPeacemanBoreholeEnergyPro
    variable = temperature
    WI = ${WI_sw}
    p_well = pp_sw_well
    block = 'sw'
  [../]

[]

[Controls]
  [./ini_kernel_1st]
    type = TimePeriod
    enable_objects = 'Kernels::source_mass_1st_mw Kernels::source_energy_1st_mw Kernels::sink_mass_1st_sw Kernels::sink_energy_1st_sw'
    disable_objects = 'Kernels::source_mass_mw Kernels::source_energy_mw Kernels::sink_mass_sw Kernels::sink_energy_sw'
    start_time = '${starting_time}'
    end_time = '${1st_step}'
    implicit = true
    reverse_on_false = true
    # execute_on = 'initial timestep_begin'
  [../]
[]

[AuxVariables]
  [velocity_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [velocity_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [velocity_z]
    family = MONOMIAL
    order = CONSTANT
  []
  [./density]
   family = MONOMIAL
   order = FIRST
  [../]

  [./viscosity]
   family = MONOMIAL
   order = FIRST
  [../]

  [./enthalpy]
   family = MONOMIAL
   order = FIRST
  [../]

  [./pp_mw_well]
  [../]

  [./rho_mw_well]
  [../]

  [./mu_mw_well]
  [../]

  [./h_mw_well]
  [../]


  [./pp_sw_well]
  [../]

  [./rho_sw_well]
  [../]

  [./mu_sw_well]
  [../]

  [./h_sw_well]
  [../]
[]

[AuxKernels]
  [velocity_x]
    type = PorousFlowDarcyVelocityComponent
    variable = velocity_x
    component = x
  []
  [velocity_y]
    type = PorousFlowDarcyVelocityComponent
    variable = velocity_y
    component = y
  []
  [velocity_z]
    type = PorousFlowDarcyVelocityComponent
    variable = velocity_z
    component = z
  []
  [./fluid_density]
   type = PorousFlowPropertyAux
   property = density
   variable = density
  [../]

  [./fluid_viscosity]
   type = PorousFlowPropertyAux
   property = viscosity
   variable = viscosity
  [../]

  [./enthalpy]
   type = PorousFlowPropertyAux
   property = enthalpy
   variable = enthalpy
  [../]
[]

[Preconditioning]
  active = 'p1'
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
  [smp]
  type = SMP
  full = true
  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
  petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
[]
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
  type = FSP
  full = true
  topsplit = pT
  [./pT]
    splitting = 'p T'
    splitting_type = multiplicative
    petsc_options_iname = '-ksp_type -pc_type -snes_type -snes_linesearch_type'
    petsc_options_value = 'fgmres lu newtonls basic'
  [../]
  [./p]
    vars = 'pressure'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type'
    petsc_options_value = 'fgmres asm ilu'
  [../]
  [./T]
    vars = 'temperature'
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type'
    petsc_options_value = 'preonly hypre boomeramg'
  [../]
[../]
[typically_efficient]
  type = SMP
  full = true
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = ' hypre    boomeramg'
[]
[strong]
  type = SMP
  full = true
  petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
  petsc_options_value = ' asm      ilu           NONZERO                   2'
[]
[probably_too_strong]
  type = SMP
  full = true
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'
[]
[]

[Executioner]
  type = Transient
  start_time =  ${starting_time}
  end_time = ${end_time}
  dtmax = ${5_day}
  l_tol = 1e-12
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_max_its = 100
  nl_max_its = 100
  solve_type = NEWTON
  fixed_point_rel_tol = 1e-5
  fixed_point_max_its = 50
  fixed_point_abs_tol = 1e-5
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 1.25
    [../]
 []

 [Outputs]
print_linear_residuals = false
[exodus]
 type = Exodus
 interval = 5
[]
# [csv]
#  type = CSV
# []
[./console]
type = Console
interval = 50
max_rows = 1
[../]
[]

# [Debug]
#   show_var_residual_norms = true
# []

[MultiApps]
  [./pipe_mw]
    type = TransientMultiApp
    input_files = pipe_mw.i
    positions = '0 0 0'
    library_path = /home/kit/agw/qz9211/projects/moskito-standard/lib
    app_type = MoskitoApp
    # execute_on = 'INITIAL TIMESTEP_BEGIN'
    # catch_up = true
    max_procs_per_app = 1
    # sub_cycling = true
    # verbose_multiapps = true
  [../]

  [./pipe_sw]
    type = TransientMultiApp
    input_files = pipe_sw.i
    positions = '0 50 0
                 0 -50 0
                 0 0 50
                 0 0 -50'
    library_path = /home/kit/agw/qz9211/projects/moskito-standard/lib
    app_type = MoskitoApp
    # execute_on = 'INITIAL TIMESTEP_BEGIN'
    # catch_up = true
    max_procs_per_app = 1
    #sub_cycling = true
    # verbose_multiapps = true
    # global_time_offset =  ${5_day}
  [../]
[]


[Transfers]
  [./to_mw_P]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = pipe_mw
    source_variable = porepressure
    variable = pp_mw_res
    from_blocks = mw
    to_blocks = 1
    execute_on = SAME_AS_MULTIAPP
    greedy_search = true
    use_nearest_app = true
    error_on_miss = true
  [../]

  [./from_mw_P]
    type = MultiAppGeneralFieldNearestNodeTransfer
    from_multi_app = pipe_mw
    source_variable = p
    variable = pp_mw_well
    from_blocks = 1
    to_blocks = mw
    execute_on = SAME_AS_MULTIAPP
  [../]

  [./from_mw_rho]
    type = MultiAppGeneralFieldNearestNodeTransfer
    from_multi_app = pipe_mw
    source_variable = density
    variable = rho_mw_well
    from_blocks = 1
    to_blocks = mw
    execute_on = SAME_AS_MULTIAPP
  [../]

  [./from_mw_mu]
    type = MultiAppGeneralFieldNearestNodeTransfer
    from_multi_app = pipe_mw
    source_variable = viscosity
    variable = mu_mw_well
    from_blocks = 1
    to_blocks = mw
    execute_on = SAME_AS_MULTIAPP
  [../]

  [./from_mw_h]
    type = MultiAppGeneralFieldNearestNodeTransfer
    from_multi_app = pipe_mw
    source_variable = h_from_p_T
    variable = h_mw_well
    from_blocks = 1
    to_blocks = mw
    execute_on = SAME_AS_MULTIAPP
  [../]


  [./to_sw_P]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = pipe_sw
    source_variable = porepressure
    variable = pp_sw_res
    from_blocks = 'sw'
    to_blocks = 1
    execute_on = SAME_AS_MULTIAPP
    greedy_search = true
    use_nearest_app = true
    error_on_miss = true
  [../]


  [./to_sw_rho]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = pipe_sw
    source_variable = density
    variable = rho_sw_res
    from_blocks = sw
    to_blocks = 1
    execute_on = SAME_AS_MULTIAPP
    greedy_search = true
    use_nearest_app = true
    error_on_miss = true
  [../]

  [./to_sw_mu]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = pipe_sw
    source_variable = viscosity
    variable = mu_sw_res
    from_blocks = sw
    to_blocks = 1
    execute_on = SAME_AS_MULTIAPP
    greedy_search = true
    use_nearest_app = true
    error_on_miss = true
  [../]

  [./to_sw_h]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    to_multi_app = pipe_sw
    source_variable = enthalpy
    variable = h_sw_res
    from_blocks = sw
    to_blocks = 1
    execute_on = SAME_AS_MULTIAPP
    greedy_search = true
    use_nearest_app = true
    error_on_miss = true
  [../]

  [./from_sw_P]
    type = MultiAppGeneralFieldNearestNodeTransfer
    from_multi_app = pipe_sw
    source_variable = p
    variable = pp_sw_well
    from_blocks = 1
    to_blocks = sw
    execute_on = SAME_AS_MULTIAPP
  [../]
[]
