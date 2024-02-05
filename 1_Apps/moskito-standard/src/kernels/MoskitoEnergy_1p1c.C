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

#include "MoskitoEnergy_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoEnergy_1p1c);

InputParameters
MoskitoEnergy_1p1c::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addParam<bool>("gravity_energy", false, "Consider potential energy "
          "caused by gravity acceleration");
  params.addClassDescription("Energy conservation equation for 1 phase (either "
          "liquid or gas) pipe flow and it returns enthalpy");
  return params;
}

MoskitoEnergy_1p1c::MoskitoEnergy_1p1c(const InputParameters & parameters)
  : Kernel(parameters),
  _q(coupledValue("flowrate")),
  // _grad_q(coupledGradient("flowrate")),
  _p(coupledValue("pressure")),
  // _grad_p(coupledGradient("pressure")),
  _q_var_number(coupled("flowrate")),
  _p_var_number(coupled("pressure")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  // _lambda(getMaterialProperty<Real>("thermal_conductivity")),
  _cp(getMaterialProperty<Real>("specific_heat")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dT(getMaterialProperty<Real>("drho_dT")),
  _h(getMaterialProperty<Real>("h_from_p_T")),
  _add_g(getParam<bool>("gravity_energy")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity"))
{
  if(_add_g)
    _gfac = 1.0;
  else
    _gfac = 0.0;
}

Real
MoskitoEnergy_1p1c::computeQpResidual()
{
  // r += _grad_test[_i][_qp] * _lambda[_qp] * _grad_u[_qp] / _cp[_qp];

  RealVectorValue r = 0.0;

  // r += _drho_dp[_qp] * _grad_p[_qp] + _drho_dT[_qp] * _grad_u[_qp];
  // r *= _q[_qp] * (_h[_qp] + _q[_qp] * _q[_qp] / 2.0 / _area[_qp] / _area[_qp]);
  // r += _grad_q[_qp] * _rho[_qp] * (_h[_qp] + 1.5 * _q[_qp] * _q[_qp] / _area[_qp] / _area[_qp]);
  // r += _rho[_qp] * _q[_qp] * (_cp[_qp] * _grad_u[_qp] - _gfac * _gravity[_qp]);
  // r *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  r += -(_rho[_qp] * (_well_sign[_qp] * _q[_qp] / _area[_qp]) * (_h[_qp] + _p[_qp] / _rho[_qp] + 0.5 * _q[_qp] * _q[_qp] / (_area[_qp] * _area[_qp]))) * _grad_test[_i][_qp];
  r += - _rho[_qp] * (_well_sign[_qp] * _q[_qp] / _area[_qp]) * _gfac * _gravity[_qp] * _test[_i][_qp];

  return r * _well_dir[_qp];
}

Real
MoskitoEnergy_1p1c::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  // j += _drho_dT[_qp] * _grad_phi[_j][_qp];
  // j *= _q[_qp] * (_h[_qp] + _q[_qp] * _q[_qp] / 2.0 / _area[_qp] / _area[_qp]);
  // j += (_drho_dp[_qp] * _grad_p[_qp] + _drho_dT[_qp] * _grad_u[_qp]) * _q[_qp] * _cp[_qp] * _phi[_j][_qp];
  // j += _grad_q[_qp] * _drho_dT[_qp] * _phi[_j][_qp] * (_h[_qp] + 1.5 * _q[_qp] * _q[_qp] / _area[_qp] / _area[_qp]);
  // j += _grad_q[_qp] * _rho[_qp] * _cp[_qp] * _phi[_j][_qp];
  // j += _drho_dT[_qp] * _phi[_j][_qp] * _q[_qp] * (_cp[_qp] * _grad_u[_qp] - _gfac * _gravity[_qp]);
  // j += _rho[_qp] * _q[_qp] * _cp[_qp] * _grad_phi[_j][_qp];
  // j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  j += -_drho_dT[_qp] * _phi[_j][_qp] * (_well_sign[_qp] * _q[_qp] / _area[_qp]) * (_h[_qp] + _p[_qp] / _rho[_qp] + 0.5 * _q[_qp] * _q[_qp] / (_area[_qp] * _area[_qp])) * _grad_test[_i][_qp];

  j += -_rho[_qp] * ((_well_sign[_qp] * _q[_qp] / _area[_qp]) * (_cp[_qp] * _phi[_j][_qp]) - _p[_qp] / (_rho[_qp] * _rho[_qp]) * _drho_dT[_qp] * _phi[_j][_qp]) * _grad_test[_i][_qp];

  j += -_drho_dT[_qp] * _phi[_j][_qp] * (_well_sign[_qp] * _q[_qp] / _area[_qp]) * _gfac * _gravity[_qp] * _test[_i][_qp];

  return j * _well_dir[_qp];
}

Real
MoskitoEnergy_1p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;

  if (jvar == _q_var_number)
  {
    // j += _drho_dp[_qp] * _grad_p[_qp] + _drho_dT[_qp] * _grad_u[_qp];
    // j *= _phi[_j][_qp] * _h[_qp] + 1.5 * _phi[_j][_qp] * _q[_qp] * _q[_qp] / _area[_qp] / _area[_qp];
    // j += _grad_phi[_j][_qp] * _rho[_qp] * (_h[_qp] + 1.5 * _q[_qp] * _q[_qp] / _area[_qp] / _area[_qp]);
    // j += _grad_q[_qp] * _rho[_qp] * 3.0 * _phi[_j][_qp] * _q[_qp] / _area[_qp] / _area[_qp];
    // j += _rho[_qp] * _phi[_j][_qp] * (_cp[_qp] * _grad_u[_qp] - _gfac * _gravity[_qp]);
    j += -((_rho[_qp] * (_h[_qp] * _phi[_j][_qp] + 1.5 * _q[_qp] * _q[_qp] / (_area[_qp] * _area[_qp]) * _phi[_j][_qp])) + _p[_qp] * _phi[_j][_qp]) * (_grad_test[_i][_qp] * _well_sign[_qp] / _area[_qp]);
    j += -_rho[_qp] * (_well_sign[_qp] * _phi[_j][_qp] / _area[_qp]) * _gfac * _gravity[_qp] * _test[_i][_qp];
  }

  if (jvar == _p_var_number)
  {
    // j += _drho_dp[_qp] * _grad_phi[_j][_qp];
    // j *= _q[_qp] * (_h[_qp] + _q[_qp] * _q[_qp] / 2.0 / _area[_qp] / _area[_qp]);
    // j += _grad_q[_qp] * _drho_dp[_qp] * _phi[_j][_qp] * (_h[_qp] + 1.5 * _q[_qp] * _q[_qp] / _area[_qp] / _area[_qp]);
    // j += _drho_dp[_qp] * _phi[_j][_qp] * _q[_qp] * (_cp[_qp] * _grad_u[_qp] - _gfac * _gravity[_qp]);
    // j *= _grad_test[_i][_qp] * _well_sign[_qp] / _area[_qp];

    j += -((_drho_dp[_qp] * _phi[_j][_qp] * _q[_qp] * (_h[_qp] + _p[_qp] / _rho[_qp] + 0.5 * _q[_qp] * _q[_qp] / (_area[_qp] * _area[_qp]))) + _q[_qp] * _phi[_j][_qp] - _p[_qp] * _q[_qp] / _rho[_qp] * _phi[_j][_qp] * _drho_dp[_qp]) * (_grad_test[_i][_qp] * _well_sign[_qp] / _area[_qp]);
    j += -_drho_dp[_qp] * _phi[_j][_qp] * (_well_sign[_qp] * _q[_qp] / _area[_qp]) * _gfac * _gravity[_qp] * _test[_i][_qp];
  }

  return j * _well_dir[_qp];
}
