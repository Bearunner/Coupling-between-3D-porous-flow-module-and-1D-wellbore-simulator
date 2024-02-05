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

#include "MoskitoLatHeatCoaxial_1p.h"

registerMooseObject("MoskitoApp", MoskitoLatHeatCoaxial_1p);

// template <>
// InputParameters
// validParams<MoskitoLatHeatCoaxial_1p>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoLatHeatCoaxial_1p::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Heat exchange between inner and outer pipes");
  params.addRequiredCoupledVar("Tsur", "temperature of other pipe");
  return params;
}

MoskitoLatHeatCoaxial_1p::MoskitoLatHeatCoaxial_1p(const InputParameters & parameters)
  : Kernel(parameters),
  _Tsur(coupledValue("Tsur")),
  _Tsur_var_number(coupled("Tsur")),
  _area(getMaterialProperty<Real>("well_area")),
  _ohc(getMaterialProperty<Real>("overall_heat_transfer_coeff"))
  {
  }

Real
MoskitoLatHeatCoaxial_1p::computeQpResidual()
{
  return _test[_i][_qp] * _ohc[_qp] * 2.0 * PI * (_u[_qp] - _Tsur[_qp]) / _area[_qp];
}

Real
MoskitoLatHeatCoaxial_1p::computeQpJacobian()
{
  return _test[_i][_qp] * _ohc[_qp] * 2.0 * PI * _phi[_j][_qp] / _area[_qp];
}

Real
MoskitoLatHeatCoaxial_1p::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real s = 0.0;
  if (jvar == _Tsur_var_number)
  {
     s -= _test[_i][_qp] * _ohc[_qp] * 2.0 * PI * _phi[_j][_qp] / _area[_qp];
  }
  return s;
}
