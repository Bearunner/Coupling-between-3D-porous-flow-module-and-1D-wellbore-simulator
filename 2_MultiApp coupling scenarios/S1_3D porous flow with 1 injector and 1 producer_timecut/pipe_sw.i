5_day = 129600   ##5 day
starting_time = 0
1st_step = 1
end_time = 9.461e+8  ## 30 years

WI = 1.87E-11

[Mesh]
  [./file]
   type = GeneratedMeshGenerator
   dim = 1
   xmin = -500
   xmax = 0
   nx = 100
  [../]
  [./source_area]
   type = SubdomainBoundingBoxGenerator
   input = file
   block_id = '1'
   bottom_left = '-400 -0.1 -0.1'
   top_right = '-500 0.1 0.1'
  [../]
[]

[GlobalParams]
  temperature = T
  pressure = p
  flowrate = q
  # execute_on = 'initial timestep_begin'
  gravity = '-9.81 0 0'
[]

[UserObjects]
  [./eos]
    # type = MoskitoEOS1P_Brine
    # specific_heat = 4186
    # type = MoskitoEOS1P_Const
    # density = 1000
    type = MoskitoEOS1P_SimpleFluidProperties
  [../]
  [./viscosity]
    type = MoskitoViscosityWaterSmith
    # type = MoskitoViscosityConst
  [../]
[]

[Materials]
  [./sw]
    type = MoskitoFluidWell_1p1c
    well_type = production
    well_direction = x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.27305
    manual_friction_factor = 0
    roughness_type = smooth
    output_properties = 'well_velocity well_P well_T well_q'
    outputs = exodus
  [../]

  [./lateral]
    type = MoskitoLatHeat_Inc_Formation_1p
    temperature_inner = T
    outer_diameters = '0.27305 0.28685 0.37465'
    conductivities = '45.0 1.1'
    convective_thermal_resistance = true
    formation_density = 2535      # Rock parameters
    formation_thermal_conductivity = 2.67
    formation_heat_capacity = 750
    formation_temperature_function = temperature_ini      # Configuration of material
    nondimensional_time_function = Ramey_1981_BF
    output_properties = 'total_thermal_resistivity'
    outputs = exodus
 [../]

[]

[BCs]
  [./top_q]
    type = MoskitoMomentumFreeBC
    variable = q
    boundary = 'right'
  [../]

  [./bot_q]
    type = DirichletBC
    variable = q
    boundary = 'left'
    value = 0
  [../]

  [./top_T]
    type = MoskitoEnergyFreeBC
    variable = T
    boundary = 'right'
  [../]

  [./bot_T]
    type = MoskitoEnergyFreeBC
    variable = T
    boundary = 'left'
  [../]

  [./top_p]
    type = FunctionDirichletBC
    variable = p
    function = hydrostatic
    boundary = 'right'
  [../]

  [./bot_p]
    type = MoskitoMassFreeBC
    variable = p
    boundary = 'left'
  [../]
[]

[ICs]
  [./flowrate_ic_sw]
    type = FunctionIC
    variable = q
    function = flowrateic_sw
  [../]
  [./pressure_ic]
    type = FunctionIC
    variable = p
    function = hydrostatic
  [../]
  [./temperature_ic]
    type = FunctionIC
    variable = T
    function = temperature_ini
  [../]
[]

[Functions]
  [./flowrateic_sw]
    type = ParsedFunction
    expression = '0'
  [../]
  [./hydrostatic]
    type = ParsedFunction
    expression = '1e5-9.81*1000*x' ##density should be close to real value!
  [../]

  [./temperature_ini]
    type = ParsedFunction
    expression = '283.15-0.03*x'
  [../]

[]

[Variables]
  [./q]
    scaling = 1e-2
  [../]
  [./p]
    scaling = 1e-1
  [../]
  [./T]
    scaling = 1e-8
  [../]
[]

[Kernels]
  [./pkernel]
    type = MoskitoMass_1p1c
    variable = p
  [../]
  [./ptimekernel]
    type = MoskitoTimeMass_1p1c
    variable = p
  [../]
  [./qkernel]
    type = MoskitoMomentum_1p1c
    variable = q
  [../]
  [./qtimekernel]
    type = MoskitoTimeMomentum_1p1c
    variable = q
  [../]
  [./Tkernel]
    type = MoskitoEnergy_1p1c
    variable = T
    gravity_energy = true
  [../]
  [./Ttimekernel]
    type = MoskitoTimeEnergy_1p1c
    variable = T
  [../]
  [./lateral_heat]
  type = MoskitoLatHeatIncFormation_1p
  variable = T
  [../]

  [./sw_mass_1st]
    type = MoskitoPeacemanBoreholeMassPro
    variable = p
    WI = ${WI}
    p_res = 1
    function = '3e5-1000*9.81*x'
    rho_res = 1000
    mu_res = 1e-3
    block = 1
  [../]

  [./sw_mass]
    type = MoskitoPeacemanBoreholeMassPro
    variable = p
    WI = ${WI}
    p_res = pp_sw_res
    rho_res = rho_sw_res
    mu_res = mu_sw_res
    block = 1
  [../]

  [./sw_energy_1st]
    type = MoskitoPeacemanBoreholeEnergyPro
    variable = T
    WI = ${WI}
    p_res = 1
    function = '3e5-1000*9.81*x'
    rho_res = 1000
    mu_res = 1e-3
    h_res = 1240000
    block = 1
  [../]

  [./sw_energy]
    type = MoskitoPeacemanBoreholeEnergyPro
    variable = T
    WI = ${WI}
    p_res = pp_sw_res
    rho_res = rho_sw_res
    mu_res = mu_sw_res
    h_res = h_sw_res
    block = 1
  [../]

[]

[Controls]
  [./ini_kernel_1st]
    type = TimePeriod
    enable_objects = 'Kernels::sw_mass_1st Kernels::sw_energy_1st'
    disable_objects = 'Kernels::sw_mass Kernels::sw_energy'
    start_time = '${starting_time}'
    end_time = '${1st_step}'
    implicit = true
    reverse_on_false = true
    # execute_on = 'initial timestep_begin'
  [../]

[]

[AuxVariables]
  [./density]
    order = FIRST
    family = MONOMIAL
  [../]
  [./viscosity]
    order = FIRST
    family = MONOMIAL
  [../]
  [./h_from_p_T]
    order = FIRST
    family = MONOMIAL
  [../]

  [./pp_sw_res]
  [../]

  [./rho_sw_res]
  [../]

  [./mu_sw_res]
  [../]

  [./h_sw_res]
  [../]

[]

[AuxKernels]
  [./rho]
   type = MaterialRealAux
   property = density
   variable = density
  [../]
  [./mu]
   type = MaterialRealAux
   property = viscosity
   variable = viscosity
  [../]
  [./h]
   type = MaterialRealAux
   property = h_from_p_T
   variable = h_from_p_T
  [../]
[]

[Preconditioning]
  active = pn1
  [./p1]
  type = SMP
  full = true
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = ' bjacobi ilu NONZERO'
  [../]
  [./pn1]
  type = SMP
  full = true
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -snes_type -snes_linesearch_type'
  petsc_options_value = ' bjacobi ilu NONZERO newtonls basic'
  [../]
  [./p2]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu NONZERO 51'
  [../]
  [./p3]
  type = SMP
  full = true
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[../]
[./p4]
  type = SMP
  full = true
  #petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
  petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
[../]
[./p5]
  type = SMP
  full = true
  petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'asm gmres lu 2 NONZERO 51'
[../]
[./p6]
  type = SMP
  full = true
  #petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -sub_pc_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'asm lu newtonls basic NONZERO 51'
[../]
[./p7]
  type = SMP
  full = true
  #petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
[../]
[]

[Executioner]
  type = Transient
  start_time =  ${starting_time}
  end_time = ${end_time}
  dtmax = ${5_day}
  # dtmin = 2
  l_tol = 1e-12
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9
  l_max_its = 100
  nl_max_its = 100
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  # automatic_scaling = true
  # compute_scaling_once = false
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
