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

#include "MoskitoLatHeat_Exc_Formation_1p.h"
#include "Conversion.h"

registerMooseObject("MoskitoApp", MoskitoLatHeat_Exc_Formation_1p);

// template <>
// InputParameters
// validParams<MoskitoLatHeat_Exc_Formation_1p>()
// {
//     InputParameters params = validParams<Material>();

InputParameters
MoskitoLatHeat_Exc_Formation_1p::validParams()
{
  InputParameters params = Material::validParams();

  //  params += validParams<NewtonIteration>();
      params += NewtonIteration::validParams();
    params.addClassDescription("Material for claculating lateral heat resistivity "
          "and it should be coupled with TIGER to get temperature at well-"
          "formation interface");
    params.addRequiredCoupledVar("temperature_inner", "Fluid temperature variable (K)");
    params.addRequiredCoupledVar("temperature_outer", "Formation temperature variable at"
          " the well-formation interface (K)");
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
    return params;
}

MoskitoLatHeat_Exc_Formation_1p::MoskitoLatHeat_Exc_Formation_1p(const InputParameters & parameters)
  : Material(parameters),
    NewtonIteration(parameters),
    _Tf(coupledValue("temperature_inner")),
    _Twf(coupledValue("temperature_outer")),
    _Rto(declareProperty<Real>("tubing_outer_radius")),
    _Uto(declareProperty<Real>("well_thermal_resistivity")),
    _tol(getParam<Real>("derivative_tolerance")),
    _OD_vec(getParam<std::vector<Real>>("outer_diameters")),
    _lambda_vec(getParam<std::vector<Real>>("conductivities")),
    _Dti(getMaterialProperty<Real>("well_diameter")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
    _Re(getMaterialProperty<Real>("well_reynolds_no")),
    _add_hf(getParam<bool>("convective_thermal_resistance"))
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

  _hf = _add_hf ? &getMaterialProperty<Real>("convective_heat_factor") : NULL;
}

void
MoskitoLatHeat_Exc_Formation_1p::computeQpProperties()
{
  if (_Dti[_qp] != _OD_vec[0])
    mooseError(name(),"The diameter provided in fluid flow material did not ",
              "match with the first tubing inner diameter provided here.\n");

  _Rto[_qp] = _OD_vec[1] / 2.0;

  if (_annulus_ind)
    returnNewtonSolve(_Uto[_qp], _Uto[_qp], _console);
  else
    _Uto[_qp] = ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf);
}

Real
MoskitoLatHeat_Exc_Formation_1p::ResistivityNoAnnulus(const int & begin, const int & end, const bool & hf)
{
  Real Uto = 0.0;
  for (int i = begin; i < end; ++i)
    if (_lambda_vec[i] != 0.0)
      Uto += std::log(_OD_vec[i+1] / _OD_vec[i]) / _lambda_vec[i];
  if (hf && _Re[_qp]>0.0)
    Uto += 1.0 / (_OD_vec[0] * 0.5 * (*_hf)[_qp]);
  Uto *= _Rto[_qp];
  Uto = 1.0 / Uto;

  return Uto;
}

Real
MoskitoLatHeat_Exc_Formation_1p::computeReferenceResidual(const Real trail_value, const Real scalar)
{
  return ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf) - scalar;
}

Real
MoskitoLatHeat_Exc_Formation_1p::computeResidual(const Real trail_value, const Real scalar)
{
  Real Uto = 0.0;
  // calculating all resisitivity except annulus; made it inverse to be able
  // to add up other commponents
  Uto = 1.0 / ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf);

  Real Ri = -1.0 / ResistivityNoAnnulus(0, _annulus_loc, _add_hf) * scalar;
  Real Ro = 1.0 / ResistivityNoAnnulus(_annulus_loc+1, _lambda_vec.size(), false) * scalar;
  Real Ti = _annulus_uo->SurfaceTemperature(_Tf[_qp], Ri, _Tf[_qp]-_Twf[_qp]);
  Real To = _annulus_uo->SurfaceTemperature(_Twf[_qp], Ro, _Tf[_qp]-_Twf[_qp]);
  Real hr = _annulus_uo->RadiativeHTCoefficient(_rai, _rao, Ti, To);
  Real hc = _annulus_uo->ConvectiveHTCoefficient(_rai, _rao, Ti, To);

  if (hc!=0.0 || hr!=0.0)
    Uto += _Rto[_qp] / (_rai * (hc + hr));

  Uto = 1.0 / Uto;

  if(std::abs((scalar-Uto)/Uto)<1e-3)
    _annulus_uo->CheckValidity(_rai, _rao, Ti, To);

  return Uto - scalar;
}

Real
MoskitoLatHeat_Exc_Formation_1p::computeDerivative(const Real trail_value, const Real scalar)
{
  Real Uto_plus_tol, Uto_minus_tol, tol_Uto;
  // Derivation concept according to MoskitoEOS2P userobject
  tol_Uto = scalar * _tol;
  Uto_plus_tol = computeResidual(trail_value, scalar + tol_Uto);
  Uto_minus_tol = computeResidual(trail_value, scalar - tol_Uto);

  return (Uto_plus_tol - Uto_minus_tol) / 2.0 / tol_Uto;
}

Real
MoskitoLatHeat_Exc_Formation_1p::initialGuess(const Real trial_value)
{
  return 0.9*ResistivityNoAnnulus(0, _lambda_vec.size(), _add_hf);
}
