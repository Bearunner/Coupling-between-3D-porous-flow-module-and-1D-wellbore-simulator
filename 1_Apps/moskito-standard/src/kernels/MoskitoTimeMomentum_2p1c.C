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

#include "MoskitoTimeMomentum_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoTimeMomentum_2p1c);

// template <>
// InputParameters
// validParams<MoskitoTimeMomentum_2p1c>()
// {
//   InputParameters params = validParams<TimeKernel>();

InputParameters
MoskitoTimeMomentum_2p1c::validParams()
{
  InputParameters params = TimeKernel::validParams();

  params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable");
  params.addClassDescription("Time derivative part of momentum conservation equation for "
                  "2 phase pipe flow and it returns pressure");

  return params;
}

MoskitoTimeMomentum_2p1c::MoskitoTimeMomentum_2p1c(const InputParameters & parameters)
  : TimeKernel(parameters),
    _m_dot(coupledDot("massrate")),
    _dm_dot(coupledDotDu("massrate")),
    _m_var_number(coupled("massrate")),
    _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
    _area(getMaterialProperty<Real>("well_area"))
{
}

Real
MoskitoTimeMomentum_2p1c::computeQpResidual()
{
  Real r = 0.0;

  r = _m_dot[_qp] * _well_sign[_qp] / _area[_qp] * _test[_i][_qp];

  return r;
}

Real
MoskitoTimeMomentum_2p1c::computeQpJacobian()
{
  return 0.0;
}

Real
MoskitoTimeMomentum_2p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _m_var_number)
  {
    j = _dm_dot[_qp] * _phi[_j][_qp] * _well_sign[_qp] / _area[_qp] * _test[_i][_qp];
  }

  return j;
}
