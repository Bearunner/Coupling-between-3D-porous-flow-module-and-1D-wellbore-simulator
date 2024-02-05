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

#include "MoskitoMass_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoMass_1p1c);

InputParameters
MoskitoMass_1p1c::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addClassDescription("Mass conservation equation for pipe flow and "
        "it returns pressure");
  return params;
}

MoskitoMass_1p1c::MoskitoMass_1p1c(const InputParameters & parameters)
  : Kernel(parameters),
  _q(coupledValue("flowrate")),
  // _grad_q(coupledGradient("flowrate")),
  // _grad_T(coupledGradient("temperature")),
  _q_var_number(coupled("flowrate")),
  _T_var_number(coupled("temperature")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dT(getMaterialProperty<Real>("drho_dT"))
{
}

Real
MoskitoMass_1p1c::computeQpResidual()
{
  RealVectorValue r = 0.0;

  // r += _drho_dp[_qp] * _grad_u[_qp];
  // r += _drho_dT[_qp] * _grad_T[_qp];
  // r *= _q[_qp];
  // r += _rho[_qp] * _grad_q[_qp];
  // r *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  r += -_rho[_qp] * _q[_qp] * _well_sign[_qp] / _area[_qp] * _grad_test[_i][_qp];

  return r * _well_dir[_qp];
}

Real
MoskitoMass_1p1c::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  // j += _drho_dp[_qp] * _grad_phi[_j][_qp] * _q[_qp];
  // j += _drho_dp[_qp] * _phi[_j][_qp] * _grad_q[_qp];
  // j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  j += -_drho_dp[_qp] * _phi[_j][_qp] * _q[_qp] * _well_sign[_qp] / _area[_qp] * _grad_test[_i][_qp];

  return j * _well_dir[_qp];
}

Real
MoskitoMass_1p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;
  if (jvar == _q_var_number)
  {
    // j += _drho_dp[_qp] * _grad_u[_qp];
    // j += _drho_dT[_qp] * _grad_T[_qp];
    // j *= _phi[_j][_qp];
    // j += _rho[_qp] * _grad_phi[_j][_qp];
    // j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

    j += -_rho[_qp] * _phi[_j][_qp] * _well_sign[_qp] / _area[_qp] * _grad_test[_i][_qp];
  }

  if (jvar == _T_var_number)
  {
    // j += _drho_dT[_qp] * _grad_phi[_j][_qp] * _q[_qp];
    // j += _drho_dT[_qp] * _phi[_j][_qp] * _grad_q[_qp];
    // j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

    j += -_drho_dT[_qp] * _phi[_j][_qp] * _q[_qp] * _well_sign[_qp] / _area[_qp] * _grad_test[_i][_qp];
  }

  return j * _well_dir[_qp];
}
