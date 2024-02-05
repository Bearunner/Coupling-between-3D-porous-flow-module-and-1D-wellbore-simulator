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

#include "MoskitoTimeEnergy_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoTimeEnergy_1p1c);

// template <>
// InputParameters
// validParams<MoskitoTimeEnergy_1p1c>()
// {
//   InputParameters params = validParams<TimeKernel>();

InputParameters
MoskitoTimeEnergy_1p1c::validParams()
{
  InputParameters params = TimeKernel::validParams();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addClassDescription("Time derivative part of energy conservation equation for "
                  "1 phase (either liquid or gas) pipe flow and it returns enthalpy");

  return params;
}

MoskitoTimeEnergy_1p1c::MoskitoTimeEnergy_1p1c(const InputParameters & parameters)
  : TimeKernel(parameters),
    _q(coupledValue("flowrate")),
    _p_dot(coupledDot("pressure")),
    _q_dot(coupledDot("flowrate")),
    _dp_dot(coupledDotDu("pressure")),
    _dq_dot(coupledDotDu("flowrate")),
    _p_var_number(coupled("pressure")),
    _q_var_number(coupled("flowrate")),
    _area(getMaterialProperty<Real>("well_area")),
    _rho(getMaterialProperty<Real>("density")),
    _h(getMaterialProperty<Real>("h_from_p_T")),
    _cp(getMaterialProperty<Real>("specific_heat")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dT(getMaterialProperty<Real>("drho_dT"))
{
}

Real
MoskitoTimeEnergy_1p1c::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _p_dot[_qp];
  r += _drho_dT[_qp] * _u_dot[_qp];
  r *= _h[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]);
  r += _rho[_qp] * (_cp[_qp] * _u_dot[_qp] + _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]));
  // r -= _p_dot[_qp]; // error?

  return r * _test[_i][_qp];
}

Real
MoskitoTimeEnergy_1p1c::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_dT[_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
  j *= (_h[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]));
  j += (_drho_dp[_qp] * _p_dot[_qp] + _drho_dT[_qp] * _u_dot[_qp]) * _cp[_qp] * _phi[_j][_qp];
  j += _drho_dT[_qp] * _phi[_j][_qp] * (_cp[_qp] * _u_dot[_qp] + _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]));
  j += _rho[_qp] * _cp[_qp] * _phi[_j][_qp] * _du_dot_du[_qp];

  return j * _test[_i][_qp];
}

Real
MoskitoTimeEnergy_1p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _q_var_number)
  {
    j += _drho_dp[_qp] * _p_dot[_qp];
    j += _drho_dT[_qp] * _u_dot[_qp];
    j *= _q[_qp] * _phi[_j][_qp];
    j += _rho[_qp] * _phi[_j][_qp] * _q_dot[_qp];

    j += _rho[_qp] * _q[_qp] * _phi[_j][_qp] * _dq_dot[_qp];

    j *= _test[_i][_qp] / (_area[_qp] * _area[_qp]);
  }

  if (jvar == _p_var_number)
  {
    j += _drho_dp[_qp] * _phi[_j][_qp] * _dp_dot[_qp];
    j *= (_h[_qp] + _q[_qp] * _q[_qp] / (2.0 * _area[_qp] * _area[_qp]));
    j += _drho_dp[_qp] * _phi[_j][_qp] * (_cp[_qp] * _u_dot[_qp] + _q[_qp] * _q_dot[_qp] / (_area[_qp] * _area[_qp]));
    // j -= _phi[_j][_qp] * _dp_dot[_qp]; // error?
    j *= _test[_i][_qp];
  }

  return j;
}
