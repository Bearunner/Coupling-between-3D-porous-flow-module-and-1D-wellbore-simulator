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

#include "MoskitoTimeTransport.h"

registerMooseObject("MoskitoApp", MoskitoTimeTransport);

// template <>
// InputParameters
// validParams<MoskitoTimeTransport>()
// {
//   InputParameters params = validParams<TimeKernel>();

InputParameters
MoskitoTimeTransport::validParams()
{
  InputParameters params = TimeKernel::validParams();

  params.addRequiredCoupledVar("enthalpy", "Temperature nonlinear variable");
  params.addRequiredCoupledVar("pressure", "pressure");
  params.addClassDescription("Time derivative part of mass conservation equation for "
                  "1 phase (either liquid or gas) pipe flow and it returns pressure");
  params.addRequiredParam<MaterialPropertyName>("aux1_name", "the name of aux1");
  params.addRequiredParam<MaterialPropertyName>("aux2_name", "the name of aux2");
  params.addRequiredParam<MaterialPropertyName>("aux1_p_name", "the name of aux1_p");
  params.addRequiredParam<MaterialPropertyName>("aux2_p_name", "the name of aux2_p");
  params.addRequiredParam<MaterialPropertyName>("aux1_h_name", "the name of aux1_h");
  params.addRequiredParam<MaterialPropertyName>("aux2_h_name", "the name of aux2_h");
  params.addRequiredParam<MaterialPropertyName>("aux1_c_name", "the name of aux1_c");
  params.addRequiredParam<MaterialPropertyName>("aux2_c_name", "the name of aux2_c");

  return params;
}

MoskitoTimeTransport::MoskitoTimeTransport(const InputParameters & parameters)
  : TimeKernel(parameters),
    _h_dot(coupledDot("enthalpy")),
    _dh_dot(coupledDotDu("enthalpy")),
    _h_var_number(coupled("enthalpy")),
    _p_dot(coupledDot("pressure")),
    _dp_dot(coupledDotDu("pressure")),
    _p_var_number(coupled("pressure")),
    _aux1(getMaterialProperty<Real>("aux1_name")),
    _aux2(getMaterialProperty<Real>("aux2_name")),
    _aux1_p(getMaterialProperty<Real>("aux1_p_name")),
    _aux2_p(getMaterialProperty<Real>("aux2_p_name")),
    _aux1_h(getMaterialProperty<Real>("aux1_h_name")),
    _aux2_h(getMaterialProperty<Real>("aux2_h_name")),
    _aux1_c(getMaterialProperty<Real>("aux1_c_name")),
    _aux2_c(getMaterialProperty<Real>("aux2_c_name"))
{
}

Real
MoskitoTimeTransport::computeQpResidual()
{
  Real r = 0.0;

  r += (_aux1_p[_qp] + _aux2_p[_qp]) * _p_dot[_qp];
  r += (_aux1_h[_qp] + _aux2_h[_qp]) * _h_dot[_qp];
  r += (_aux1_c[_qp] + _aux2_c[_qp]) * _u_dot[_qp];
  r *= _u[_qp];
  r += (_aux1[_qp] + _aux2[_qp]) * _u_dot[_qp];
  r *= _test[_i][_qp];

  return r;
}

Real
MoskitoTimeTransport::computeQpJacobian()
{
  Real j = 0.0;

  j += (_aux1_p[_qp] + _aux2_p[_qp]) * _p_dot[_qp];
  j += (_aux1_h[_qp] + _aux2_h[_qp]) * _h_dot[_qp];
  j += (_aux1_c[_qp] + _aux2_c[_qp]) * _u_dot[_qp];
  j *= _phi[_j][_qp];
  j += (_aux1_c[_qp] + _aux2_c[_qp]) * _phi[_j][_qp] * _du_dot_du[_qp] * _u[_qp];
  j += (_aux1_c[_qp] + _aux2_c[_qp]) * _phi[_j][_qp] * _u_dot[_qp];
  j += (_aux1[_qp] + _aux2[_qp]) * _phi[_j][_qp] * _du_dot_du[_qp];
  j *= _test[_i][_qp];

  return j;
}

Real
MoskitoTimeTransport::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _h_var_number)
  {
    j += (_aux1_h[_qp] + _aux2_h[_qp]) * _phi[_j][_qp] * _dh_dot[_qp] * _u[_qp];
    j += (_aux1_h[_qp] + _aux2_h[_qp]) * _phi[_j][_qp] * _u_dot[_qp];
    j *= _test[_i][_qp];
  }

  if (jvar == _p_var_number)
  {
    j += (_aux1_p[_qp] + _aux2_p[_qp]) * _phi[_j][_qp] * _dp_dot[_qp] * _u[_qp];
    j += (_aux1_p[_qp] + _aux2_p[_qp]) * _phi[_j][_qp] * _u_dot[_qp];
    j *= _test[_i][_qp];
  }

  return j;
}
