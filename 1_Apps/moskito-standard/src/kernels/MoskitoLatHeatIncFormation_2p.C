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

#include "MoskitoLatHeatIncFormation_2p.h"

registerMooseObject("MoskitoApp", MoskitoLatHeatIncFormation_2p);

// template <>
// InputParameters
// validParams<MoskitoLatHeatIncFormation_2p>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoLatHeatIncFormation_2p::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("molarity", "gas concentration");
  params.addClassDescription("Lateral heat exchange between wellbore "
        "and formation including formation");
  return params;
}

MoskitoLatHeatIncFormation_2p::MoskitoLatHeatIncFormation_2p(const InputParameters & parameters)
  : Kernel(parameters),
  _p_var_number(coupled("pressure")),
  _c_var_number(coupled("molarity")),
  _lambda(getMaterialProperty<Real>("total_thermal_resistivity")),
  _Tform(getMaterialProperty<Real>("formation_temperature")),
  _area(getMaterialProperty<Real>("well_area")),
  _T(getMaterialProperty<Real>("temperature")),
  _T_h(getMaterialProperty<Real>("T_h")),
  _T_p(getMaterialProperty<Real>("T_p")),
  _T_c(getMaterialProperty<Real>("T_c"))
{
}

Real
MoskitoLatHeatIncFormation_2p::computeQpResidual()
{
  Real r = 0.0;
  r =  2.0 * PI * _lambda[_qp] * (_Tform[_qp] - _T[_qp] );
  //std::cout<<_lambda[_qp]<<std::endl;
  r /=  _area[_qp];

  return  -1.0 * r * _test[_i][_qp];
}

Real
MoskitoLatHeatIncFormation_2p::computeQpJacobian()
{

  Real j = 0.0;
  j =  -2.0 * PI * _lambda[_qp] * _phi[_j][_qp] * _T_h[_qp];
  j /=  _area[_qp];

  return  -1.0 * j * _test[_i][_qp];

}

Real
MoskitoLatHeatIncFormation_2p::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;


  if (jvar == _p_var_number)
  {

    j +=  -2.0 * PI * _lambda[_qp] * _phi[_j][_qp] * _T_p[_qp];
    j /=  _area[_qp];

  }

  if (jvar == _c_var_number)
  {

    j +=  -2.0 * PI * _lambda[_qp] * _phi[_j][_qp] * _T_c[_qp];
    j /=  _area[_qp];

  }

  return -1.0 * j * _test[_i][_qp];
}
