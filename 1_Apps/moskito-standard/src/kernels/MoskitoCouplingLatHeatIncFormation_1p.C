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

#include "MoskitoCouplingLatHeatIncFormation_1p.h"

registerMooseObject("MoskitoApp", MoskitoCouplingLatHeatIncFormation_1p);

// template <>
// InputParameters
// validParams<MoskitoCouplingLatHeatIncFormation_1p>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoCouplingLatHeatIncFormation_1p::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("T_res", "T_res");
  params.addClassDescription("Lateral heat exchange between wellbore "
        "and formation including formation");
  return params;
}

MoskitoCouplingLatHeatIncFormation_1p::MoskitoCouplingLatHeatIncFormation_1p(const InputParameters & parameters)
  : Kernel(parameters),
  _lambda(getMaterialProperty<Real>("total_thermal_resistivity")),
  _T_res(coupledValue("T_res")),
  _area(getMaterialProperty<Real>("well_area"))
{
}

Real
MoskitoCouplingLatHeatIncFormation_1p::computeQpResidual()
{
  Real r = 0.0;
  r =  2.0 * PI * _lambda[_qp] * (_T_res[_qp] - _u[_qp]);
  r /=  _area[_qp];

  return  -1.0 * r * _test[_i][_qp];
}

Real
MoskitoCouplingLatHeatIncFormation_1p::computeQpJacobian()
{
  Real j = 0.0;
  j =  -2.0 * PI * _lambda[_qp] * _phi[_j][_qp];
  j /=  _area[_qp];

  return  -1.0 * j * _test[_i][_qp];
}
