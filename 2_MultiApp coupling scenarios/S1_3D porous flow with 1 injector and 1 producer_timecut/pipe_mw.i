5_day = 129600   ##5 day
starting_time = 0
1st_step = 1
end_time = 9.461e+8  ## 30 years

WI = 1.46E-11

T_inj_hot = 363.15

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
  [./mw]
    type = MoskitoFluidWell_1p1c
    well_type = injection
    well_direction = -x
    eos_uo = eos
    viscosity_uo = viscosity
    well_diameter = 0.32269545
    manual_friction_factor = 0
    roughness_type = smooth
    output_properties = 'well_velocity well_P well_T well_q'
    outputs = exodus
  [../]

  [./lateral]
   type = MoskitoLatHeat_Inc_Formation_1p
   temperature_inner = T
   outer_diameters = '0.32269545 0.339004545454545 0.442768181818182'
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
    type = DirichletBC
    variable = T
    value = '${T_inj_hot}'
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
  [./flowrate_ic_mw]
    type = FunctionIC
    variable = q
    function = flowrateic_mw
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
  [./flowrateic_mw]
    type = ParsedFunction
    expression = '0'
  [../]
  [./hydrostatic]
    type = ParsedFunction
    expression = '25e5-9.81*1000*x'
  [../]
  [./temperature_ini]
    type = ParsedFunction
    expression = '283.15-0.03*x'
  [../]
[]

[Variables]
  [./q]
    scaling = 1e-3
  [../]
  [./p]
    scaling = 1e-2
  [../]
  [./T]
    scaling = 1e-9
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

  [./mw_mass_1st]
    type = MoskitoPeacemanBoreholeMass
    variable = p
    WI = ${WI}
    p_res = 1
    function = '3e5-1000*9.81*x'
    block = 1
  [../]

  [./mw_mass]
    type = MoskitoPeacemanBoreholeMass
    variable = p
    WI = ${WI}
    p_res = pp_mw_res
    block = 1
  [../]

  [./mw_energy_1st]
      type = MoskitoPeacemanBoreholeEnergy
      variable = T
      WI = ${WI}
      p_res = 1
      function = '3e5-1000*9.81*x'
      block = 1
    [../]

  [./mw_energy]
    type = MoskitoPeacemanBoreholeEnergy
    variable = T
    WI = ${WI}
    p_res = pp_mw_res
    block = 1
  [../]
[]

[Controls]
  [./ini_kernel_1st]
    type = TimePeriod
    enable_objects = 'Kernels::mw_mass_1st Kernels::mw_energy_1st'
    disable_objects = 'Kernels::mw_mass Kernels::mw_energy'
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

  [./pp_mw_res]
  [../]

  [./rho_mw_res]
  [../]

  [./mu_mw_res]
  [../]

  [./h_mw_res]
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

[./console]
type = Console
interval = 50
max_rows = 1
[../]
[]

# [Debug]
#   show_var_residual_norms = true
# []
