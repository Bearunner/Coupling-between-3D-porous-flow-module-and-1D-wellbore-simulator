/**************************************************************************/
/*  Moskito - THMC sImulator for GEoscience Research                        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of Moskito App                                        */
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

#include "MoskitoMassFreeBC.h"

registerMooseObject("MoskitoApp", MoskitoMassFreeBC);

InputParameters
MoskitoMassFreeBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addClassDescription(
    "Implements free advective flow boundary conditions for the Mass equation.");
  return params;
}

MoskitoMassFreeBC::MoskitoMassFreeBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
  _q(coupledValue("flowrate")),
  _T_var_number(coupled("temperature")),
  _q_var_number(coupled("flowrate")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dT(getMaterialProperty<Real>("drho_dT")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  _area(getMaterialProperty<Real>("well_area"))
{
}

Real
MoskitoMassFreeBC::computeQpResidual()
{
  return  _rho[_qp] * (_q[_qp] *  _well_sign[_qp] / _area[_qp]) * _normals[_qp] * _test[_i][_qp] * _well_dir[_qp];
}

Real
MoskitoMassFreeBC::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  j += _drho_dp[_qp] * _phi[_j][_qp] * _q[_qp] * _well_sign[_qp] * _normals[_qp] / _area[_qp] * _test[_i][_qp];

  return j * _well_dir[_qp];
}

Real
MoskitoMassFreeBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;

  if (jvar == _q_var_number)
  {
    j += _rho[_qp] * _phi[_j][_qp] * _well_sign[_qp] * _normals[_qp] / _area[_qp] * _test[_i][_qp];
  }

  if (jvar == _T_var_number)
  {
    j += _drho_dT[_qp] * _phi[_j][_qp] * _q[_qp] * _well_sign[_qp] * _normals[_qp] / _area[_qp] * _test[_i][_qp];
  }

  return j * _well_dir[_qp];
}
