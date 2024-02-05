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

#include "MoskitoTransport.h"

registerMooseObject("MoskitoApp", MoskitoTransport);

// template <>
// InputParameters
// validParams<MoskitoTransport>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoTransport::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable");
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addClassDescription("Transport");
  params.addRequiredParam<MaterialPropertyName>("aux1_name", "the name of aux1");
  params.addRequiredParam<MaterialPropertyName>("aux2_name", "the name of aux2");
  params.addRequiredParam<MaterialPropertyName>("aux1_p_name", "the name of aux1_p");
  params.addRequiredParam<MaterialPropertyName>("aux2_p_name", "the name of aux2_p");
  params.addRequiredParam<MaterialPropertyName>("aux1_h_name", "the name of aux1_h");
  params.addRequiredParam<MaterialPropertyName>("aux2_h_name", "the name of aux2_h");
  params.addRequiredParam<MaterialPropertyName>("aux1_c_name", "the name of aux1_c");
  params.addRequiredParam<MaterialPropertyName>("aux2_c_name", "the name of aux2_c");
  params.addRequiredParam<MaterialPropertyName>("u_l_c_name", "the name of u_l_c");
  params.addRequiredParam<MaterialPropertyName>("u_g_c_name", "the name of u_g_c");

  return params;
}

MoskitoTransport::MoskitoTransport(const InputParameters & parameters)
  : Kernel(parameters),
  _grad_m(coupledGradient("massrate")),
  _m_var_number(coupled("massrate")),
  _grad_p(coupledGradient("pressure")),
  _p_var_number(coupled("pressure")),
  _grad_h(coupledGradient("enthalpy")),
  _h_var_number(coupled("enthalpy")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  _aux1(getMaterialProperty<Real>("aux1_name")),
  _aux2(getMaterialProperty<Real>("aux2_name")),
  _aux1_p(getMaterialProperty<Real>("aux1_p_name")),
  _aux2_p(getMaterialProperty<Real>("aux2_p_name")),
  _u_l_p(getMaterialProperty<Real>("u_l_p")),
  _u_g_p(getMaterialProperty<Real>("u_g_p")),
  _aux1_h(getMaterialProperty<Real>("aux1_h_name")),
  _aux2_h(getMaterialProperty<Real>("aux2_h_name")),
  _u_l_h(getMaterialProperty<Real>("u_l_h")),
  _u_g_h(getMaterialProperty<Real>("u_g_h")),
  _aux1_c(getMaterialProperty<Real>("aux1_c_name")),
  _aux2_c(getMaterialProperty<Real>("aux2_c_name")),
  _u_l_c(getMaterialProperty<Real>("u_l_c_name")),
  _u_g_c(getMaterialProperty<Real>("u_g_c_name")),
  _u_l_m(getMaterialProperty<Real>("u_l_m")),
  _u_g_m(getMaterialProperty<Real>("u_g_m")),
  _u_g(getMaterialProperty<Real>("gas_velocity")),
  _u_l(getMaterialProperty<Real>("liquid_velocity"))
{
}

Real
MoskitoTransport::computeQpResidual()
{

  RealVectorValue r = 0.0;
  RealVectorValue e = 0.0;

  r += _aux1_p[_qp] * _grad_p[_qp];
  r += _aux1_h[_qp] * _grad_h[_qp];
  r += _aux1_c[_qp] * _grad_u[_qp];
  r *= _u_g[_qp];
  r += _aux1[_qp] * (_u_g_p[_qp]*_grad_p[_qp] + _u_g_h[_qp]*_grad_h[_qp] + _u_g_c[_qp]*_grad_u[_qp] + _u_g_m[_qp]*_grad_m[_qp]);
  r *= _u[_qp];
  r += _aux1[_qp] * _u_g[_qp] * _grad_u[_qp];
  r *= _test[_i][_qp] * _well_sign[_qp];

  e += _aux2_p[_qp] * _grad_p[_qp];
  e += _aux2_h[_qp] * _grad_h[_qp];
  e += _aux2_c[_qp] * _grad_u[_qp];
  e *= _u_l[_qp];
  e += _aux2[_qp] * (_u_l_p[_qp]*_grad_p[_qp] + _u_l_h[_qp]*_grad_h[_qp] + _u_l_c[_qp]*_grad_u[_qp] + _u_l_m[_qp]*_grad_m[_qp]);
  e *= _u[_qp];
  e += _aux2[_qp] * _u_l[_qp] * _grad_u[_qp];
  e *= _test[_i][_qp] * _well_sign[_qp];

  return (r + e) * _well_dir[_qp];
}

Real
MoskitoTransport::computeQpJacobian()
{

  RealVectorValue j = 0.0;
  RealVectorValue h = 0.0;
  RealVectorValue s = 0.0;

  RealVectorValue a = 0.0;
  RealVectorValue b = 0.0;
  RealVectorValue c = 0.0;

  j += _aux1_p[_qp] * _grad_p[_qp];
  j += _aux1_h[_qp] * _grad_h[_qp];
  j += _aux1_c[_qp] * _grad_u[_qp];
  j *= _phi[_j][_qp];
  j += _aux1_c[_qp] * _grad_phi[_j][_qp] * _u[_qp];
  j *= _u_g[_qp];
  j += (_aux1_p[_qp]*_grad_p[_qp] + _aux1_h[_qp]*_grad_h[_qp] + _aux1_c[_qp]*_grad_u[_qp])*_u_g_c[_qp]*_phi[_j][_qp]*_u[_qp];

  a += _aux2_p[_qp] * _grad_p[_qp];
  a += _aux2_h[_qp] * _grad_h[_qp];
  a += _aux2_c[_qp] * _grad_u[_qp];
  a *= _phi[_j][_qp];
  a += _aux2_c[_qp] * _grad_phi[_j][_qp] * _u[_qp];
  a *= _u_l[_qp];
  a += (_aux2_p[_qp]*_grad_p[_qp] + _aux2_h[_qp]*_grad_h[_qp] + _aux2_c[_qp]*_grad_u[_qp])*_u_l_c[_qp]*_phi[_j][_qp]*_u[_qp];
  //////////////////////////////////////////////////////////////////////////////

  h += _aux1_c[_qp] * (_u_g_p[_qp]*_grad_p[_qp] + _u_g_h[_qp]*_grad_h[_qp] + _u_g_c[_qp]*_grad_u[_qp] + _u_g_m[_qp]*_grad_m[_qp]) * _u[_qp];
  h += _aux1[_qp] * (_u_g_p[_qp]*_grad_p[_qp] + _u_g_h[_qp]*_grad_h[_qp] + _u_g_c[_qp]*_grad_u[_qp] + _u_g_m[_qp]*_grad_m[_qp]);
  h *= _phi[_j][_qp];
  h += _aux1[_qp]*_u_g_c[_qp]*_grad_phi[_j][_qp]*_u[_qp];

  b += _aux2_c[_qp] * (_u_l_p[_qp]*_grad_p[_qp] + _u_l_h[_qp]*_grad_h[_qp] + _u_l_c[_qp]*_grad_u[_qp] + _u_l_m[_qp]*_grad_m[_qp]) * _u[_qp];
  b += _aux2[_qp] * (_u_l_p[_qp]*_grad_p[_qp] + _u_l_h[_qp]*_grad_h[_qp] + _u_l_c[_qp]*_grad_u[_qp] + _u_l_m[_qp]*_grad_m[_qp]);
  b *= _phi[_j][_qp];
  b += _aux2[_qp]*_u_l_c[_qp]*_grad_phi[_j][_qp]*_u[_qp];
  //////////////////////////////////////////////////////////////////////////////

  s += _aux1_c[_qp] * _phi[_j][_qp] * _grad_u[_qp];
  s += _aux1[_qp] * _grad_phi[_j][_qp];
  s *= _u_g[_qp];
  s += _aux1[_qp]*_u_g_c[_qp]*_phi[_j][_qp]*_grad_u[_qp];

  c += _aux2_c[_qp] * _phi[_j][_qp] * _grad_u[_qp];
  c += _aux2[_qp] * _grad_phi[_j][_qp];
  c *= _u_l[_qp];
  c += _aux2[_qp]*_u_l_c[_qp]*_phi[_j][_qp]*_grad_u[_qp];
  //////////////////////////////////////////////////////////////////////////////

  j += h + s;
  j *= _test[_i][_qp] * _well_sign[_qp];

  a += b + c;
  a *= _test[_i][_qp] * _well_sign[_qp];
  //////////////////////////////////////////////////////////////////////////////

  return (j + a) * _well_dir[_qp];
}

Real
MoskitoTransport::computeQpOffDiagJacobian(unsigned int jvar)
{

  RealVectorValue j = 0.0;
  RealVectorValue k = 0.0;
  if (jvar == _m_var_number)
  {
    j += _aux1_p[_qp] * _grad_p[_qp];
    j += _aux1_h[_qp] * _grad_h[_qp];
    j += _aux1_c[_qp] * _grad_u[_qp];
    j *= _u_g_m[_qp] * _phi[_j][_qp] * _u[_qp];

    j += _aux1[_qp] * _u_g_m[_qp] * _grad_phi[_j][_qp] * _u[_qp];

    j += _aux1[_qp] * _u_g_m[_qp] * _phi[_j][_qp] * _grad_u[_qp];

    j *= _test[_i][_qp] * _well_sign[_qp];
    ////////////////////////////////////////////////////////////////////////////
    k += _aux2_p[_qp] * _grad_p[_qp];
    k += _aux2_h[_qp] * _grad_h[_qp];
    k += _aux2_c[_qp] * _grad_u[_qp];
    k *= _u_l_m[_qp] * _phi[_j][_qp] * _u[_qp];

    k += _aux2[_qp] * _u_l_m[_qp] * _grad_phi[_j][_qp] * _u[_qp];

    k += _aux2[_qp] * _u_l_m[_qp] * _phi[_j][_qp] * _grad_u[_qp];

    k *= _test[_i][_qp] * _well_sign[_qp];
  }

  if (jvar == _h_var_number)
  {
    j += _aux1_h[_qp] * _grad_phi[_j][_qp] * _u_g[_qp];

    j += _aux1_h[_qp] * _phi[_j][_qp] * (_u_g_p[_qp]*_grad_p[_qp] + _u_g_h[_qp]*_grad_h[_qp] + _u_g_c[_qp]*_grad_u[_qp] + _u_g_m[_qp]*_grad_m[_qp]);
    j *= _u[_qp];

    j += _aux1_h[_qp] * _phi[_j][_qp] * _u_g[_qp] * _grad_u[_qp];

    j += (_aux1_p[_qp]*_grad_p[_qp] + _aux1_h[_qp]*_grad_h[_qp] + _aux1_c[_qp]*_grad_u[_qp])*_u_g_h[_qp]*_phi[_j][_qp]*_u[_qp];

    j += _aux1_h[_qp]*_u_g_h[_qp]*_grad_phi[_j][_qp]*_u[_qp];

    j += _aux1_h[_qp]*_u_g_h[_qp]*_phi[_j][_qp] * _grad_u[_qp];

    j *= _test[_i][_qp] * _well_sign[_qp];
    ///////////////////////////////////////////////////////////////////////////
    k += _aux2_h[_qp] * _grad_phi[_j][_qp] * _u_l[_qp];

    k += _aux2_h[_qp] * _phi[_j][_qp] * (_u_l_p[_qp]*_grad_p[_qp] + _u_l_h[_qp]*_grad_h[_qp] + _u_l_c[_qp]*_grad_u[_qp] + _u_l_m[_qp]*_grad_m[_qp]);
    k *= _u[_qp];

    k += _aux2_h[_qp] * _phi[_j][_qp] * _u_l[_qp] * _grad_u[_qp];

    k += (_aux2_p[_qp]*_grad_p[_qp] + _aux2_h[_qp]*_grad_h[_qp] + _aux2_c[_qp]*_grad_u[_qp])*_u_l_h[_qp]*_phi[_j][_qp]*_u[_qp];

    k += _aux2_h[_qp]*_u_l_h[_qp]*_grad_phi[_j][_qp]*_u[_qp];

    k += _aux2_h[_qp]*_u_l_h[_qp]*_phi[_j][_qp] * _grad_u[_qp];

    k *= _test[_i][_qp] * _well_sign[_qp];
  }

  if (jvar == _p_var_number)
  {
    j += _aux1_p[_qp] * _grad_phi[_j][_qp] * _u_g[_qp];

    j += _aux1_p[_qp] * _phi[_j][_qp] * (_u_g_p[_qp]*_grad_p[_qp] + _u_g_h[_qp]*_grad_h[_qp] + _u_g_c[_qp]*_grad_u[_qp] + _u_g_m[_qp]*_grad_m[_qp]);
    j *= _u[_qp];

    j += _aux1_p[_qp] * _phi[_j][_qp] * _u_g[_qp] * _grad_u[_qp];

    j += (_aux1_p[_qp]*_grad_p[_qp] + _aux1_h[_qp]*_grad_h[_qp] + _aux1_c[_qp]*_grad_u[_qp])*_u_g_p[_qp]*_phi[_j][_qp]*_u[_qp];

    j += _aux1_p[_qp]*_u_g_p[_qp]*_grad_phi[_j][_qp]*_u[_qp];

    j += _aux1_p[_qp]*_u_g_p[_qp]*_phi[_j][_qp] * _grad_u[_qp];

    j *= _test[_i][_qp] * _well_sign[_qp];
    ////////////////////////////////////////////////////////////////////////////
    k += _aux2_p[_qp] * _grad_phi[_j][_qp] * _u_l[_qp];

    k += _aux2_p[_qp] * _phi[_j][_qp] * (_u_l_p[_qp]*_grad_p[_qp] + _u_l_h[_qp]*_grad_h[_qp] + _u_l_c[_qp]*_grad_u[_qp] + _u_l_m[_qp]*_grad_m[_qp]);
    k *= _u[_qp];

    k += _aux2_p[_qp] * _phi[_j][_qp] * _u_l[_qp] * _grad_u[_qp];

    k += (_aux2_p[_qp]*_grad_p[_qp] + _aux2_h[_qp]*_grad_h[_qp] + _aux2_c[_qp]*_grad_u[_qp])*_u_l_p[_qp]*_phi[_j][_qp]*_u[_qp];

    k += _aux2_p[_qp]*_u_l_p[_qp]*_grad_phi[_j][_qp]*_u[_qp];

    k += _aux2_p[_qp]*_u_l_p[_qp]*_phi[_j][_qp] * _grad_u[_qp];

    k *= _test[_i][_qp] * _well_sign[_qp];
  }

  return (j + k) * _well_dir[_qp];
}
