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

#include "MoskitoGrad.h"

registerMooseObject("MoskitoApp", MoskitoGrad);

// template <>
// InputParameters
// validParams<MoskitoGrad>()
// {
//     InputParameters params = validParams<Material>();

InputParameters
MoskitoGrad::validParams()
{
  InputParameters params = Material::validParams();

    params.addRequiredCoupledVar("concentration", "concentration");
    params.addCoupledVar("c2", "Coupled value c2");
    params.addCoupledVar("c3", "Coupled value c3");
    params.addCoupledVar("c4", "Coupled value c4");
    params.addCoupledVar("c5", "Coupled value c5");
    params.addCoupledVar("c6", "Coupled value c6");
    params.addCoupledVar("c7", "Coupled value c7");
    params.addCoupledVar("c8", "Coupled value c8");

    //params.addRequiredParam<std::string>("u_g_c_name", "The name of u_g_c");


    return params;
}

MoskitoGrad::MoskitoGrad(const InputParameters & parameters)
  : Material(parameters),
  _drho_m_dc(declareProperty<Real>("drho_dc")),
  _drho_m_dcvect(declareProperty<RealVectorValue>("drho_dcvect")),
  _dgamma_dc(declareProperty<Real>("dgamma_dc")),
  _dgamma_dcvect(declareProperty<RealVectorValue>("dgamma_dcvect")),
  _dkappa_dc(declareProperty<Real>("dkappa_dc")),
  _dkappa_dcvect(declareProperty<RealVectorValue>("dkappa_dcvect")),
  _domega_dc(declareProperty<Real>("domega_dc")),
  _domega_dcvect(declareProperty<RealVectorValue>("domega_dcvect")),
  _T_c(declareProperty<Real>("T_c")),
  _T_cvect(declareProperty<RealVectorValue>("T_cvect")),
  _drho_cdot(declareProperty<Real>("drho_cdot")),
  _dgamma_dcdot(declareProperty<Real>("dgamma_dcdot")),

  _drho_dc1(getMaterialProperty<Real>("m5")),
  //_drho_dc2(getMaterialProperty<Real>("s5")),
  _drho_dc2(isCoupled("c2") ? getMaterialProperty<Real>("s5") : getMaterialProperty<Real>("zero_0")),
  _drho_dc3(isCoupled("c3") ? getMaterialProperty<Real>("a5") : getMaterialProperty<Real>("zero_0")),
  _drho_dc4(isCoupled("c4") ? getMaterialProperty<Real>("b5") : getMaterialProperty<Real>("zero_0")),
  _drho_dc5(isCoupled("c5") ? getMaterialProperty<Real>("c5") : getMaterialProperty<Real>("zero_0")),
  _drho_dc6(isCoupled("c6") ? getMaterialProperty<Real>("d5") : getMaterialProperty<Real>("zero_0")),
  _drho_dc7(isCoupled("c7") ? getMaterialProperty<Real>("e5") : getMaterialProperty<Real>("zero_0")),
  _drho_dc8(isCoupled("c8") ? getMaterialProperty<Real>("f5") : getMaterialProperty<Real>("zero_0")),

  _dgamma_dc1(getMaterialProperty<Real>("m6")),
  //_dgamma_dc2(getMaterialProperty<Real>("s6")),
  _dgamma_dc2(isCoupled("c2") ? getMaterialProperty<Real>("s6") : getMaterialProperty<Real>("zero_0")),
  _dgamma_dc3(isCoupled("c3") ? getMaterialProperty<Real>("a6") : getMaterialProperty<Real>("zero_0")),
  _dgamma_dc4(isCoupled("c4") ? getMaterialProperty<Real>("b6") : getMaterialProperty<Real>("zero_0")),
  _dgamma_dc5(isCoupled("c5") ? getMaterialProperty<Real>("c6") : getMaterialProperty<Real>("zero_0")),
  _dgamma_dc6(isCoupled("c6") ? getMaterialProperty<Real>("d6") : getMaterialProperty<Real>("zero_0")),
  _dgamma_dc7(isCoupled("c7") ? getMaterialProperty<Real>("e6") : getMaterialProperty<Real>("zero_0")),
  _dgamma_dc8(isCoupled("c8") ? getMaterialProperty<Real>("f6") : getMaterialProperty<Real>("zero_0")),


  _dkappa_dc1(getMaterialProperty<Real>("m7")),
  //_dkappa_dc2(getMaterialProperty<Real>("s7")),
  _dkappa_dc2(isCoupled("c2") ? getMaterialProperty<Real>("s7") : getMaterialProperty<Real>("zero_0")),
  _dkappa_dc3(isCoupled("c3") ? getMaterialProperty<Real>("a7") : getMaterialProperty<Real>("zero_0")),
  _dkappa_dc4(isCoupled("c4") ? getMaterialProperty<Real>("b7") : getMaterialProperty<Real>("zero_0")),
  _dkappa_dc5(isCoupled("c5") ? getMaterialProperty<Real>("c7") : getMaterialProperty<Real>("zero_0")),
  _dkappa_dc6(isCoupled("c6") ? getMaterialProperty<Real>("d7") : getMaterialProperty<Real>("zero_0")),
  _dkappa_dc7(isCoupled("c7") ? getMaterialProperty<Real>("e7") : getMaterialProperty<Real>("zero_0")),
  _dkappa_dc8(isCoupled("c8") ? getMaterialProperty<Real>("f7") : getMaterialProperty<Real>("zero_0")),


  _domega_dc1(getMaterialProperty<Real>("m8")),
  //_domega_dc2(getMaterialProperty<Real>("s8")),
  _domega_dc2(isCoupled("c2") ? getMaterialProperty<Real>("s8") : getMaterialProperty<Real>("zero_0")),
  _domega_dc3(isCoupled("c3") ? getMaterialProperty<Real>("a8") : getMaterialProperty<Real>("zero_0")),
  _domega_dc4(isCoupled("c4") ? getMaterialProperty<Real>("b8") : getMaterialProperty<Real>("zero_0")),
  _domega_dc5(isCoupled("c5") ? getMaterialProperty<Real>("c8") : getMaterialProperty<Real>("zero_0")),
  _domega_dc6(isCoupled("c6") ? getMaterialProperty<Real>("d8") : getMaterialProperty<Real>("zero_0")),
  _domega_dc7(isCoupled("c7") ? getMaterialProperty<Real>("e8") : getMaterialProperty<Real>("zero_0")),
  _domega_dc8(isCoupled("c8") ? getMaterialProperty<Real>("f8") : getMaterialProperty<Real>("zero_0")),


  _T_c1(getMaterialProperty<Real>("m9")),
  //_T_c2(getMaterialProperty<Real>("s9")),
  _T_c2(isCoupled("c2") ? getMaterialProperty<Real>("s9") : getMaterialProperty<Real>("zero_0")),
  _T_c3(isCoupled("c3") ? getMaterialProperty<Real>("a9") : getMaterialProperty<Real>("zero_0")),
  _T_c4(isCoupled("c4") ? getMaterialProperty<Real>("b9") : getMaterialProperty<Real>("zero_0")),
  _T_c5(isCoupled("c5") ? getMaterialProperty<Real>("c9") : getMaterialProperty<Real>("zero_0")),
  _T_c6(isCoupled("c6") ? getMaterialProperty<Real>("d9") : getMaterialProperty<Real>("zero_0")),
  _T_c7(isCoupled("c7") ? getMaterialProperty<Real>("e9") : getMaterialProperty<Real>("zero_0")),
  _T_c8(isCoupled("c8") ? getMaterialProperty<Real>("f9") : getMaterialProperty<Real>("zero_0")),


  _c(coupledValue("concentration")),
  _c2(isCoupled("c2") ? coupledValue("c2") : _zero),
  _c3(isCoupled("c3") ? coupledValue("c3") : _zero),
  _c4(isCoupled("c4") ? coupledValue("c4") : _zero),
  _c5(isCoupled("c5") ? coupledValue("c5") : _zero),
  _c6(isCoupled("c6") ? coupledValue("c6") : _zero),
  _c7(isCoupled("c7") ? coupledValue("c7") : _zero),
  _c8(isCoupled("c8") ? coupledValue("c8") : _zero),

  _c_dot(coupledDot("concentration")),
  //_c2_dot(coupledDot("c2")),
  _c2_dot(isCoupled("c2") ? coupledDot("c2") : _zero),
  _c3_dot(isCoupled("c3") ? coupledDot("c3") : _zero),
  _c4_dot(isCoupled("c4") ? coupledDot("c4") : _zero),
  _c5_dot(isCoupled("c5") ? coupledDot("c5") : _zero),
  _c6_dot(isCoupled("c6") ? coupledDot("c6") : _zero),
  _c7_dot(isCoupled("c7") ? coupledDot("c7") : _zero),
  _c8_dot(isCoupled("c8") ? coupledDot("c8") : _zero),

  //_grad_c2(isCoupled("c2") ? coupledGradient("c2") : _zero),
  _grad_c(coupledGradient("concentration")),
  _grad_c2(coupledGradient("c2")),
  _grad_c3(coupledGradient("c3")),
  _grad_c4(coupledGradient("c4")),
  _grad_c5(coupledGradient("c5")),
  _grad_c6(coupledGradient("c6")),
  _grad_c7(coupledGradient("c7")),
  _grad_c8(coupledGradient("c8"))
{
}

void
MoskitoGrad::computeQpProperties()
{

  //_drho_m_dc[_qp] = _drho_dc3[_qp] ;//+ _drho_dc2[_qp];
  _drho_m_dcvect[_qp] = _drho_dc1[_qp]*_grad_c[_qp] + _drho_dc2[_qp]*_grad_c2[_qp] + _drho_dc3[_qp]*_grad_c3[_qp] + _drho_dc4[_qp]*_grad_c4[_qp] +
   _drho_dc5[_qp]*_grad_c5[_qp] + _drho_dc6[_qp]*_grad_c6[_qp] + _drho_dc7[_qp]*_grad_c7[_qp] + _drho_dc8[_qp]*_grad_c8[_qp];

  _dgamma_dcvect[_qp] = _dgamma_dc1[_qp]*_grad_c[_qp] + _dgamma_dc2[_qp]*_grad_c2[_qp] + _dgamma_dc3[_qp]*_grad_c3[_qp] + _dgamma_dc4[_qp]*_grad_c4[_qp] +
   _dgamma_dc5[_qp]*_grad_c5[_qp] + _dgamma_dc6[_qp]*_grad_c6[_qp] + _dgamma_dc7[_qp]*_grad_c7[_qp] + _dgamma_dc8[_qp]*_grad_c8[_qp];

  _dkappa_dcvect[_qp] = _dkappa_dc1[_qp]*_grad_c[_qp] + _dkappa_dc2[_qp]*_grad_c2[_qp] + _dkappa_dc3[_qp]*_grad_c3[_qp] + _dkappa_dc4[_qp]*_grad_c4[_qp] +
   + _dkappa_dc5[_qp]*_grad_c5[_qp] + _dkappa_dc6[_qp]*_grad_c6[_qp] + _dkappa_dc7[_qp]*_grad_c7[_qp] + _dkappa_dc8[_qp]*_grad_c8[_qp];

  _domega_dcvect[_qp] = _domega_dc1[_qp]*_grad_c[_qp] + _domega_dc2[_qp]*_grad_c2[_qp] + _domega_dc3[_qp]*_grad_c3[_qp] + _domega_dc4[_qp]*_grad_c4[_qp] +
   + _domega_dc5[_qp]*_grad_c5[_qp] + _domega_dc6[_qp]*_grad_c6[_qp] + _domega_dc7[_qp]*_grad_c7[_qp] + _domega_dc8[_qp]*_grad_c8[_qp];

  _T_c[_qp] = _T_c1[_qp] ;//+ _T_c2[_qp];
  _T_cvect[_qp] = _T_c1[_qp]*_grad_c[_qp] + _T_c2[_qp]*_grad_c2[_qp];

  _drho_cdot[_qp] = _drho_dc1[_qp]*_c_dot[_qp] + _drho_dc2[_qp]*_c2_dot[_qp] + _drho_dc3[_qp]*_c3_dot[_qp] + _drho_dc4[_qp]*_c4_dot[_qp] +
   + _drho_dc5[_qp]*_c5_dot[_qp] + _drho_dc6[_qp]*_c6_dot[_qp] + _drho_dc7[_qp]*_c7_dot[_qp] + _drho_dc8[_qp]*_c8_dot[_qp];

  _dgamma_dcdot[_qp] = _dgamma_dc1[_qp]*_c_dot[_qp] + _dgamma_dc2[_qp]*_c2_dot[_qp] + _dgamma_dc3[_qp]*_c3_dot[_qp] + _dgamma_dc4[_qp]*_c4_dot[_qp] +
   + _dgamma_dc5[_qp]*_c5_dot[_qp] + _dgamma_dc6[_qp]*_c6_dot[_qp] + _dgamma_dc7[_qp]*_c7_dot[_qp] + _dgamma_dc8[_qp]*_c8_dot[_qp];

}
