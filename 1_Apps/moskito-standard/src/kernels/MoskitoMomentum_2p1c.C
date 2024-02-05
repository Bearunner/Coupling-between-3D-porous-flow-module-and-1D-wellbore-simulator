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

#include "MoskitoMomentum_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoMomentum_2p1c);

// template <>
// InputParameters
// validParams<MoskitoMomentum_2p1c>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoMomentum_2p1c::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable");
  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable");
  params.addRequiredCoupledVar("molarity", "gas concentration");
  params.addClassDescription("Momentum conservation equation for 2 phase "
        "pipe flow and it returns pressure.");
  return params;
}

MoskitoMomentum_2p1c::MoskitoMomentum_2p1c(const InputParameters & parameters)
  : Kernel(parameters),
    _m(coupledValue("massrate")),
    _grad_m(coupledGradient("massrate")),
    _m_var_number(coupled("massrate")),
    _grad_h(coupledGradient("enthalpy")),
    _h_var_number(coupled("enthalpy")),
    _grad_c(coupledGradient("molarity")),
    _c_var_number(coupled("molarity")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dp2(getMaterialProperty<Real>("drho_dp2")),
    _drho_dh(getMaterialProperty<Real>("drho_dh")),
    _drho_dh2(getMaterialProperty<Real>("drho_dh2")),
    _drho_dph(getMaterialProperty<Real>("drho_dph")),
    _drho_dc(getMaterialProperty<Real>("drho_dc")),
    _drho_m_dcvect(getMaterialProperty<RealVectorValue>("drho_dcvect")),
    _dgamma_dcvect(getMaterialProperty<RealVectorValue>("dgamma_dcvect")),
    _f(getMaterialProperty<Real>("well_moody_friction")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _area(getMaterialProperty<Real>("well_area")),
    _perimeter(getMaterialProperty<Real>("well_perimeter")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
    _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
    _dgamma_dp(getMaterialProperty<Real>("dgamma_dp")),
    _dgamma_dh(getMaterialProperty<Real>("dgamma_dh")),
    _dgamma_dm(getMaterialProperty<Real>("dgamma_dm")),
    _dgamma_dc(getMaterialProperty<Real>("dgamma_dc"))
{
}

Real
MoskitoMomentum_2p1c::computeQpResidual()
{
  RealVectorValue r = 0.0;
  Real a2 = _area[_qp] * _area[_qp];

  //r -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_dc[_qp] * _grad_c[_qp])
  r -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_m_dcvect[_qp])
        * _m[_qp] * _m[_qp] / (_rho[_qp] * _rho[_qp] * a2);
  r += 2.0 * _m[_qp] * _grad_m[_qp] / (_rho[_qp] * a2);
  r += _well_sign[_qp] * _well_dir[_qp] * _f[_qp] * _m[_qp] * _m[_qp]
        * _perimeter[_qp] / (8.0 * _area[_qp] * _rho[_qp] * a2);
  r += _grad_u[_qp];
  r -= _rho[_qp] * _gravity[_qp];
  r += _dgamma_dh[_qp] * _grad_h[_qp];
  r += _dgamma_dp[_qp] * _grad_u[_qp];
  r += _dgamma_dm[_qp] * _grad_m[_qp];
  //r += _dgamma_dc[_qp] * _grad_c[_qp];
  r += _dgamma_dcvect[_qp];

  return r * _well_dir[_qp] * _test[_i][_qp];
}

Real
MoskitoMomentum_2p1c::computeQpJacobian()
{
  RealVectorValue j = 0.0;
  Real a2 = _area[_qp] * _area[_qp];

  j -= (_drho_dp2[_qp] *_phi[_j][_qp] * _grad_u[_qp] + _drho_dp[_qp]
        * _grad_phi[_j][_qp] + _drho_dph[_qp] *_phi[_j][_qp] * _grad_h[_qp])
        * _m[_qp] * _m[_qp] / (_rho[_qp] * _rho[_qp] * a2);
  //j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_dc[_qp] * _grad_c[_qp]) * _m[_qp]
  j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_m_dcvect[_qp]) * _m[_qp]
        * _m[_qp] * -2.0 * _drho_dp[_qp] * _phi[_j][_qp] / (_rho[_qp] * _rho[_qp]
        * _rho[_qp] * a2);
  j += 2.0 * _m[_qp] * _grad_m[_qp] * -_drho_dp[_qp] *_phi[_j][_qp]
        / (_rho[_qp] * _rho[_qp] * a2);
  j += _well_sign[_qp] * _well_dir[_qp] * _f[_qp] * _m[_qp] * _m[_qp]
        * _perimeter[_qp] * -_drho_dp[_qp] * _phi[_j][_qp] / (8.0 * _area[_qp]
        * _rho[_qp] * _rho[_qp] * a2);
  j += _grad_phi[_j][_qp];
  j -= _drho_dp[_qp] * _phi[_j][_qp] * _gravity[_qp];
  j += _dgamma_dp[_qp] * _grad_phi[_j][_qp];

  return j * _well_dir[_qp] * _test[_i][_qp];
}

Real
MoskitoMomentum_2p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;
  Real a2 = _area[_qp] * _area[_qp];

  if (jvar == _m_var_number)
  {
    //j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_dc[_qp] * _grad_c[_qp])
    j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_m_dcvect[_qp])
          * 2.0 * _m[_qp] * _phi[_j][_qp] / (_rho[_qp] * _rho[_qp] * a2);
    j += 2.0 * (_phi[_j][_qp] * _grad_m[_qp] + _m[_qp] * _grad_phi[_j][_qp])
          / (_rho[_qp] * a2);
    j += _well_sign[_qp] * _well_dir[_qp] * _f[_qp] * 2.0 * _m[_qp]
          * _phi[_j][_qp] * _perimeter[_qp] / (8.0 * _area[_qp] * _rho[_qp] * a2);
    j += _dgamma_dm[_qp] * _grad_phi[_j][_qp];
  }

  if (jvar == _h_var_number)
  {
    j -= (_drho_dph[_qp] *_phi[_j][_qp] * _grad_u[_qp] + _drho_dh[_qp]
          * _grad_phi[_j][_qp] + _drho_dh2[_qp] *_phi[_j][_qp] * _grad_h[_qp])
          * _m[_qp] * _m[_qp] / (_rho[_qp] * _rho[_qp] * a2);
    //j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_dc[_qp] * _grad_c[_qp]) * _m[_qp]
    j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_m_dcvect[_qp]) * _m[_qp]
          * _m[_qp] * -2.0 * _drho_dh[_qp] *_phi[_j][_qp] / (_rho[_qp] * _rho[_qp]
          * _rho[_qp] * a2);
    j += 2.0 * _m[_qp] * _grad_m[_qp] * -_drho_dh[_qp] *_phi[_j][_qp]
          / (_rho[_qp] * _rho[_qp] * a2);
    j += _well_sign[_qp] * _well_dir[_qp] * _f[_qp] * _m[_qp] * _m[_qp]
          * _perimeter[_qp] * -_drho_dh[_qp] *_phi[_j][_qp] / (8.0 * _area[_qp]
          * _rho[_qp] * _rho[_qp] * a2);
    j -= _drho_dh[_qp] * _phi[_j][_qp] * _gravity[_qp];
    j += _dgamma_dh[_qp] * _grad_phi[_j][_qp];
  }
/*
  if (jvar == _c_var_number)
  {
    j -= ( _drho_dc[_qp]  * _grad_phi[_j][_qp] )
          * _m[_qp] * _m[_qp] / (_rho[_qp] * _rho[_qp] * a2);
    j -= (_drho_dp[_qp] * _grad_u[_qp] + _drho_dh[_qp] * _grad_h[_qp] + _drho_dc[_qp] * _grad_c[_qp]) * _m[_qp]
          * _m[_qp] * -2.0 * _drho_dc[_qp] *_phi[_j][_qp] / (_rho[_qp] * _rho[_qp]
          * _rho[_qp] * a2);
    j += 2.0 * _m[_qp] * _grad_m[_qp] * -_drho_dc[_qp] *_phi[_j][_qp]
          / (_rho[_qp] * _rho[_qp] * a2);
    j += _well_sign[_qp] * _well_dir[_qp] * _f[_qp] * _m[_qp] * _m[_qp]
          * _perimeter[_qp] * -_drho_dc[_qp] *_phi[_j][_qp] / (8.0 * _area[_qp]
          * _rho[_qp] * _rho[_qp] * a2);
    j -= _drho_dc[_qp] * _phi[_j][_qp] * _gravity[_qp];
    j += _dgamma_dc[_qp] * _grad_phi[_j][_qp];
  }
*/
  return j * _well_dir[_qp] * _test[_i][_qp];
}
