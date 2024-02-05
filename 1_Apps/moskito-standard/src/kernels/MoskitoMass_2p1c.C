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

#include "MoskitoMass_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoMass_2p1c);

// template <>
// InputParameters
// validParams<MoskitoMass_2p1c>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoMass_2p1c::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Mass conservation equation for 2 phase pipe flow and "
        "it returns massrate");

  return params;
}

MoskitoMass_2p1c::MoskitoMass_2p1c(const InputParameters & parameters)
  : Kernel(parameters),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign"))
{
}

Real
MoskitoMass_2p1c::computeQpResidual()
{
  Real r = 0.0;

  r += _grad_u[_qp] * _well_dir[_qp];
  r *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  return r;
}

Real
MoskitoMass_2p1c::computeQpJacobian()
{
  Real j = 0.0;

  j += _grad_phi[_j][_qp] * _well_dir[_qp];
  j *= _test[_i][_qp] * _well_sign[_qp] / _area[_qp];

  return j;
}
