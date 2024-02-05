/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*  Co-developed by Sebastian Held                                        */
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

#include "MoskitoLatHeatIncFormation_1p.h"

registerMooseObject("MoskitoApp", MoskitoLatHeatIncFormation_1p);

// template <>
// InputParameters
// validParams<MoskitoLatHeatIncFormation_1p>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoLatHeatIncFormation_1p::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Lateral heat exchange between wellbore "
        "and formation including formation");
  return params;
}

MoskitoLatHeatIncFormation_1p::MoskitoLatHeatIncFormation_1p(const InputParameters & parameters)
  : Kernel(parameters),
  _lambda(getMaterialProperty<Real>("total_thermal_resistivity")),
  _Tform(getMaterialProperty<Real>("formation_temperature")),
  _area(getMaterialProperty<Real>("well_area"))
{
}

Real
MoskitoLatHeatIncFormation_1p::computeQpResidual()
{
  Real r = 0.0;
  r =  2.0 * PI * _lambda[_qp] * (_Tform[_qp] - _u[_qp] );
  r /=  _area[_qp];

  return  -1.0 * r * _test[_i][_qp];
}

Real
MoskitoLatHeatIncFormation_1p::computeQpJacobian()
{
  Real j = 0.0;
  j =  -2.0 * PI * _lambda[_qp] * _phi[_j][_qp];
  j /=  _area[_qp];

  return  -1.0 * j * _test[_i][_qp];
}
