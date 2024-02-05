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

#include "MoskitoLatHeatExcFormation_1p.h"

registerMooseObject("MoskitoApp", MoskitoLatHeatExcFormation_1p);

// template <>
// InputParameters
// validParams<MoskitoLatHeatExcFormation_1p>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoLatHeatExcFormation_1p::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Lateral heat exchange between wellbore "
        "and formation; it should take temperature at the well-formation"
        " interface from TIGER");
  params.addRequiredCoupledVar("temperature_outer", "Formation temperature variable at"
        " the well-formation interface (K)");
  return params;
}

MoskitoLatHeatExcFormation_1p::MoskitoLatHeatExcFormation_1p(const InputParameters & parameters)
  : Kernel(parameters),
  _rto(getMaterialProperty<Real>("radius_tubbing_outer")),
  _Uto(getMaterialProperty<Real>("thermal_resistivity_well")),
  _Twf(coupledValue("temperature_outer")),
  _Twf_var_number(coupled("temperature_outer")),
  _area(getMaterialProperty<Real>("well_area"))
{
}

Real
MoskitoLatHeatExcFormation_1p::computeQpResidual()
{
  Real r = 0.0;
  r =  2.0 * PI * _rto[_qp] * _Uto[_qp] * (_Twf[_qp] - _u[_qp]);
  r /=  _area[_qp];

  return  -1.0 * r * _test[_i][_qp];
}

Real
MoskitoLatHeatExcFormation_1p::computeQpJacobian()
{
  Real j = 0.0;
  j =  -2.0 * PI * _rto[_qp] * _Uto[_qp] * _phi[_j][_qp];
  j /=  _area[_qp];

  return  -1.0 * j * _test[_i][_qp];
}

Real
MoskitoLatHeatExcFormation_1p::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _Twf_var_number)
  {
    j =  2.0 * PI * _rto[_qp] * _Uto[_qp] * _phi[_j][_qp];
    j /=  _area[_qp];
  }

  return  -1.0 * j * _test[_i][_qp];
}
