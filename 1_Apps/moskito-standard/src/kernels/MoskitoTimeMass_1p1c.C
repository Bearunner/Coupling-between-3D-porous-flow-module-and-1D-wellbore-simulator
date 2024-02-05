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

#include "MoskitoTimeMass_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoTimeMass_1p1c);

// template <>
// InputParameters
// validParams<MoskitoTimeMass_1p1c>()
// {
//   InputParameters params = validParams<TimeKernel>();

InputParameters
MoskitoTimeMass_1p1c::validParams()
{
  InputParameters params = TimeKernel::validParams();

  params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable");
  params.addClassDescription("Time derivative part of mass conservation equation for "
                  "1 phase (either liquid or gas) pipe flow and it returns pressure");

  return params;
}

MoskitoTimeMass_1p1c::MoskitoTimeMass_1p1c(const InputParameters & parameters)
  : TimeKernel(parameters),
    _T_dot(coupledDot("temperature")),
    _dT_dot(coupledDotDu("temperature")),
    _T_var_number(coupled("temperature")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dT(getMaterialProperty<Real>("drho_dT"))
{
}

Real
MoskitoTimeMass_1p1c::computeQpResidual()
{
  Real r = 0.0;

  r += _drho_dp[_qp] * _u_dot[_qp];
  r += _drho_dT[_qp] * _T_dot[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoTimeMass_1p1c::computeQpJacobian()
{
  Real j = 0.0;

  j += _drho_dp[_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoTimeMass_1p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _T_var_number)
  {
    j += _drho_dT[_qp] * _phi[_j][_qp] * _dT_dot[_qp];
    j *= _test[_i][_qp];
  }

  return j;
}
