/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#include "MoskitoLatHeat_Inc_Formation_2p.h"
#include "Conversion.h"

registerMooseObject("MoskitoApp", MoskitoLatHeat_Inc_Formation_2p);

// template <>
// InputParameters
// validParams<MoskitoLatHeat_Inc_Formation_2p>()
// {
//     InputParameters params = validParams<Material>();

InputParameters
MoskitoLatHeat_Inc_Formation_2p::validParams()
{
  InputParameters params = Material::validParams();

  //  params += validParams<NewtonIteration>();
    params += NewtonIteration::validParams();
    params.addClassDescription("Material for claculating lateral heat resistivity "
          "and it should be coupled with TIGER to get temperature at well-"
          "formation interface");
    //params.addRequiredCoupledVar("temperature_inner", "Fluid temperature variable (K)");
    params.addParam<Real>("derivative_tolerance", 0.001,
          "Tolerance to calculate derivatives based on numerical differentiation");
    params.addRequiredParam<std::vector<Real>>("conductivities",
          "Vector containing conductivity values of full well completion (W/(m*K))."
          "for considering an annulus, zero should be set and annulus_uo should "
          "be provided. Only one annulus can be modelled!!!");
    params.addRequiredParam<std::vector<Real>>("outer_diameters",
          "Vector containing outer diameter values of full well completion [m];"
          "including first tubing inner diameter");
    params.addParam<UserObjectName>("annulus_uo", "",
          "The name of the userobject for annulus");
    params.addParam<bool>("convective_thermal_resistance", false, "Consider thermal"
          " resistance caused by convective term at the wall of the inner tubing");
    params.addRequiredParam<Real>("formation_heat_capacity",
          "Specific heat capacity of the formation (J/(kg*K))");
    params.addRequiredParam<Real>("formation_density",
          "Density of the formation (kg/m^3))");
    params.addRequiredParam<Real>("formation_thermal_conductivity",
          "Thermal conductivity of the formation (W/(m*K))");
    params.addRequiredParam<FunctionName>("formation_temperature_function",
          "Static formation temperature as a function of depth and time (K)");
    params.addParam<Real>("manual_time", 86400,
          "Time defined by the user for steady state simulations (s)");
    MooseEnum dt_model
          ("Kutasov2003_constHF Kutasov2005_constWFT Ramey_1981_BF Hasan_Kabir_2012", "Ramey_1981_BF");
    params.addParam<MooseEnum>("nondimensional_time_function", dt_model,
          "Select a function from the list (Kutasov2003 with constant heat flux,"
          " Kutasov2005 with constant well-formation temperature, Ramey 1981 "
          "best fit, Hasan&Kabir2012)");

    return params;
}

MoskitoLatHeat_Inc_Formation_2p::MoskitoLatHeat_Inc_Formation_2p(const InputParameters & parameters)
  : Material(parameters),
    NewtonIteration(parameters),
    _Tf(getMaterialProperty<Real>("temperature")),
    _Uto(declareProperty<Real>("well_thermal_resistivity")),
    _lambda_t(declareProperty<Real>("total_thermal_resistivity")),
    _tol(getParam<Real>("derivative_tolerance")),
    _OD_vec(getParam<std::vector<Real>>("outer_diameters")),
    _lambda_vec(getParam<std::vector<Real>>("conductivities")),
    _Tform_func(getFunction("formation_temperature_function")),
    _Tform(declareProperty<Real>("formation_temperature")),
    _lambda_form(getParam<Real>("formation_thermal_conductivity")),
    _rho_form(getParam<Real>("formation_density")),
    _cp_form(getParam<Real>("formation_heat_capacity")),
    _Dti(getMaterialProperty<Real>("well_diameter")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
    _Re(getMaterialProperty<Real>("well_reynolds_no")),
    _add_hf(getParam<bool>("convective_thermal_resistance")),
    _time_func_user(getParam<MooseEnum>("nondimensional_time_function"))
{
  for (int i = 0; i < _OD_vec.size()-1 ; i++)
    if (_OD_vec[i] >= _OD_vec[i+1])
      mooseError(name()," Layers with zero thikness detected, check outer_dimeters.\n");

  if (_OD_vec.size() != _lambda_vec.size() + 1)
    mooseError(name(),"Vector conductivities has to be in size 1 element "
              "smaller then vector of the completion diameters.\n");

  for (int i = 0; i < _lambda_vec.size(); i++)
    if (_lambda_vec[i] == 0.0)
    {
      _annulus_ind = true;
      _annulus_loc = i;
      _rai = _OD_vec[i] / 2.0;
      _rao = _OD_vec[i+1] / 2.0;
    }

  if (parameters.isParamSetByUser("annulus_uo")*_annulus_ind)
    _annulus_uo = &getUserObject<MoskitoAnnulus>("annulus_uo");
  else if (_annulus_ind || parameters.isParamSetByUser("annulus_uo"))
    mooseError(name()," Either annulus_uo or zero conductivity is not provided to model annulus.\n");
  else
    _annulus_uo = NULL;

  if (parameters.isParamSetByUser("manual_time"))
  {
    _mt = true;
    _time = getParam<Real>("manual_time");
  }

  _alpha_form = _lambda_form / _rho_form / _cp_form;
  _Rto = _OD_vec[1] / 2.0;
  _Rwf = _OD_vec.back() / 2.0;

  _hf = _add_hf ? &getMaterialProperty<Real>("convective_heat_factor") : NULL;
}

void
MoskitoLatHeat_Inc_Formation_2p::computeQpProperties()
{
  if (_Dti[_qp] != _OD_vec[0])
    mooseError(name(),"The diameter provided in fluid flow material did not ",
              "match with the first tubing inner diameter provided here.\n");

  if(!_mt) _time = _t;
  _Tform[_qp] = _Tform_func.value(_t, _q_point[_qp]);
  _ft = nonDtimefunction();

  if (_annulus_ind)
    returnNewtonSolve(_Uto[_qp], _Uto[_qp], _console);
  else
    _Uto[_qp] = ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf);

  _lambda_t[_qp]  = _lambda_form  * _Rto * _Uto[_qp];
  _lambda_t[_qp] /= _Rto * _Uto[_qp] * _ft + _lambda_form;
}

Real
MoskitoLatHeat_Inc_Formation_2p::ResistivityNoAnnulus(const int & begin, const int & end, const bool & hf)
{
  Real Uto = 0.0;
  for (int i = begin; i < end; ++i)
    if (_lambda_vec[i] != 0.0)
      Uto += std::log(_OD_vec[i+1] / _OD_vec[i]) / _lambda_vec[i];
  if (hf && _Re[_qp]>0.0)
    Uto += 1.0 / (_OD_vec[0] * 0.5 * (*_hf)[_qp]);
  Uto *= _Rto;
  Uto = 1.0 / Uto;

  return Uto;
}

Real
MoskitoLatHeat_Inc_Formation_2p::TemperatureWFinterface(const Real & Uto)
{
  Real Twf = 0.0;
  Twf += _Rto * Uto * _ft * _Tf[_qp] + _lambda_form * _Tform[_qp];
  Twf /= _Rto * Uto * _ft + _lambda_form;

  return Twf;
}

Real
MoskitoLatHeat_Inc_Formation_2p::nonDtimefunction()
{
  Real ft = 0.0;

  Real t_fac = _alpha_form * _time / _Rwf / _Rwf;

  switch (_time_func_user)
  {
  case time_func_cases::Kutasov2003_constHF:
    if(t_fac<0.194598)
      ft = 2.0 * std::sqrt(t_fac / PI);
    else
      ft = std::log(1.0 + (1.781 - 1.0 / (2.701 + std::sqrt(t_fac))) * std::sqrt(t_fac));
  break;

  case time_func_cases::Kutasov2005_constWFT:
    ft = std::log(1.0 + (1.571 - 1.0 / (4.959 + std::sqrt(t_fac))) * std::sqrt(t_fac));
  break;

  case time_func_cases::Ramey_1981_BF:
    // fit with 1% deviation from Ramey_1981
    ft = std::log(1.0 + 1.7 * std::sqrt(t_fac));
  break;

  case time_func_cases::Hasan_Kabir_2012:
    ft = std::log(std::exp(-0.2*t_fac) + (1.5 - 0.3719 * std::exp(-t_fac)) * std::sqrt(t_fac));
  break;
  }

  return ft;
}

Real
MoskitoLatHeat_Inc_Formation_2p::computeReferenceResidual(const Real trial_value, const Real scalar)
{
  return ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf) - scalar;
}

Real
MoskitoLatHeat_Inc_Formation_2p::computeResidual(const Real trial_value, const Real scalar)
{
  Real Uto = 0.0;
  // calculating all resisitivity except annulus; made it inverse to be able
  // to add up other commponents
  Uto = 1.0 / ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf);

  Real Twf = TemperatureWFinterface(scalar);
  Real Ri = -1.0 / ResistivityNoAnnulus(0, _annulus_loc, _add_hf) * scalar;
  Real Ro = 1.0 / ResistivityNoAnnulus(_annulus_loc+1, _lambda_vec.size(), false) * scalar;
  Real Ti = _annulus_uo->SurfaceTemperature(_Tf[_qp], Ri, _Tf[_qp]-Twf);
  Real To = _annulus_uo->SurfaceTemperature(Twf, Ro, _Tf[_qp]-Twf);
  Real hr = _annulus_uo->RadiativeHTCoefficient(_rai, _rao, Ti, To);
  Real hc = _annulus_uo->ConvectiveHTCoefficient(_rai, _rao, Ti, To);

  if (hc!=0.0 || hr!=0.0)
    Uto += _Rto / (_rai * (hc + hr));

  Uto = 1.0 / Uto;

  if(std::abs((scalar-Uto)/Uto)<1e-3)
    _annulus_uo->CheckValidity(_rai, _rao, Ti, To);

  return Uto - scalar;
}

Real
MoskitoLatHeat_Inc_Formation_2p::computeDerivative(const Real trial_value, const Real scalar)
{
  Real Uto_plus_tol, Uto_minus_tol, tol_Uto;
  // Derivation concept according to MoskitoEOS2P userobject
  tol_Uto = scalar * _tol;
  Uto_plus_tol = computeResidual(trial_value, scalar + tol_Uto);
  Uto_minus_tol = computeResidual(trial_value, scalar - tol_Uto);

  return (Uto_plus_tol - Uto_minus_tol) / 2.0 / tol_Uto;
}

Real
MoskitoLatHeat_Inc_Formation_2p::initialGuess(const Real trial_value)
{
  return 0.9*ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf);
}
