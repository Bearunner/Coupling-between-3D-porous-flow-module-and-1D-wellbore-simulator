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

#include "MoskitoEnergy_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoEnergy_2p1c);

// template <>
// InputParameters
// validParams<MoskitoEnergy_2p1c>()
// {
//   InputParameters params = validParams<Kernel>();

InputParameters
MoskitoEnergy_2p1c::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable");
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("molarity", "gas concentration");
  params.addParam<bool>("gravity_energy", false, "Consider potential energy "
          "caused by gravity acceleration");
  params.addClassDescription("Energy conservation equation for 2 phase "
          "pipe flow and it returns enthalpy");

  return params;
}

MoskitoEnergy_2p1c::MoskitoEnergy_2p1c(const InputParameters & parameters)
  : Kernel(parameters),
  _m(coupledValue("massrate")),
  _grad_m(coupledGradient("massrate")),
  _m_var_number(coupled("massrate")),
  _grad_p(coupledGradient("pressure")),
  _p_var_number(coupled("pressure")),
  _grad_c(coupledGradient("molarity")),
  _c_var_number(coupled("molarity")),
  _area(getMaterialProperty<Real>("well_area")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  _drho_m_dcvect(getMaterialProperty<RealVectorValue>("drho_dcvect")),
  _dgamma_dcvect(getMaterialProperty<RealVectorValue>("dgamma_dcvect")),
  _dkappa_dcvect(getMaterialProperty<RealVectorValue>("dkappa_dcvect")),
  _domega_dcvect(getMaterialProperty<RealVectorValue>("domega_dcvect")),
  _T_cvect(getMaterialProperty<RealVectorValue>("T_cvect")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  _lambda(getMaterialProperty<Real>("thermal_conductivity")),
  _cp(getMaterialProperty<Real>("specific_heat")),
  _rho(getMaterialProperty<Real>("density")),
  _drho_dp(getMaterialProperty<Real>("drho_dp")),
  _drho_dh(getMaterialProperty<Real>("drho_dh")),
  _drho_dc(getMaterialProperty<Real>("drho_dc")),
  _drho_dp2(getMaterialProperty<Real>("drho_dp2")),
  _drho_dh2(getMaterialProperty<Real>("drho_dh2")),
  _drho_dph(getMaterialProperty<Real>("drho_dph")),
  _add_g(getParam<bool>("gravity_energy")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity")),
  _dkappa_dp(getMaterialProperty<Real>("dkappa_dp")),
  _dkappa_dh(getMaterialProperty<Real>("dkappa_dh")),
  _dkappa_dm(getMaterialProperty<Real>("dkappa_dm")),
  _dkappa_dc(getMaterialProperty<Real>("dkappa_dc")),
  _dkappa_dph(getMaterialProperty<Real>("dkappa_dph")),
  _dkappa_dpm(getMaterialProperty<Real>("dkappa_dpm")),
  _dkappa_dhm(getMaterialProperty<Real>("dkappa_dhm")),
  _dkappa_dp2(getMaterialProperty<Real>("dkappa_dp2")),
  _dkappa_dh2(getMaterialProperty<Real>("dkappa_dh2")),
  _dkappa_dm2(getMaterialProperty<Real>("dkappa_dm2")),
  _domega_dp(getMaterialProperty<Real>("domega_dp")),
  _domega_dh(getMaterialProperty<Real>("domega_dh")),
  _domega_dm(getMaterialProperty<Real>("domega_dm")),
  _domega_dc(getMaterialProperty<Real>("domega_dc"))
{
  if(_add_g)
    _gfac = 1.0;
  else
    _gfac = 0.0;
}

Real
MoskitoEnergy_2p1c::computeQpResidual()
{
  RealVectorValue r = 0.0;

  r += (_u[_qp] + 1.5 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0))
        * _grad_m[_qp];
  r -= _gfac * _m[_qp] * _gravity[_qp];
  r /= _area[_qp];
  r += (_m[_qp] / _area[_qp] - std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0)
        * _drho_dh[_qp]) * _grad_u[_qp];
  r -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * _drho_dp[_qp]
        * _grad_p[_qp];
  //r -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * _drho_dc[_qp]
  //      * _grad_c[_qp];
  r -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * _drho_m_dcvect[_qp];
  r += _dkappa_dh[_qp] * _grad_u[_qp];
  r += _dkappa_dp[_qp] * _grad_p[_qp];
  //r += _dkappa_dc[_qp] * _grad_c[_qp];
  r += _dkappa_dcvect[_qp];
  r += _dkappa_dm[_qp] * _grad_m[_qp];
  //r += _domega_dh[_qp] * _grad_u[_qp] + _domega_dp[_qp] * _grad_p[_qp] + _domega_dm[_qp] * _grad_m[_qp] + _domega_dc[_qp] * _grad_c[_qp];
  r += _domega_dh[_qp] * _grad_u[_qp] + _domega_dp[_qp] * _grad_p[_qp] + _domega_dm[_qp] * _grad_m[_qp] + _domega_dcvect[_qp];

  return r * _test[_i][_qp] * _well_dir[_qp] * _well_sign[_qp];
}

Real
MoskitoEnergy_2p1c::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  j += (_phi[_j][_qp] - 3.0 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0)
        * _drho_dh[_qp] * _phi[_j][_qp] / _rho[_qp]) * _grad_m[_qp] / _area[_qp];
  j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (_drho_dh2[_qp]
        - 3.0 * _drho_dh[_qp] * _drho_dh[_qp] / _rho[_qp]) * _phi[_j][_qp]
             * _grad_u[_qp];
  j += (_m[_qp] / _area[_qp] - std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0)
              * _drho_dh[_qp]) * _grad_phi[_j][_qp];
  j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (_drho_dph[_qp]
        - 3.0 * _drho_dp[_qp] * _drho_dh[_qp] / _rho[_qp]) * _phi[_j][_qp]
        * _grad_p[_qp];
  //j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (0.0
  //      - 3.0 * _drho_dc[_qp] * _drho_dh[_qp] / _rho[_qp]) * _phi[_j][_qp]
  //      * _grad_c[_qp];
  j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (0.0
        - 3.0 * _drho_dh[_qp] / _rho[_qp]) * _phi[_j][_qp]
        * _drho_m_dcvect[_qp];
  j += _dkappa_dh[_qp] * _grad_phi[_j][_qp];
  // j += _dkappa_dph[_qp] * _phi[_j][_qp] * _grad_p[_qp];
  j += _domega_dh[_qp] * _grad_phi[_j][_qp];

  return j * _test[_i][_qp] * _well_dir[_qp] * _well_sign[_qp];
}

Real
MoskitoEnergy_2p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;

  if (jvar == _m_var_number)
  {
    j += 3.0 * std::pow(1.0 / (_rho[_qp] * _area[_qp]), 2.0) * _m[_qp]
          * _phi[_j][_qp] * _grad_m[_qp];
    j += (_u[_qp] + 1.5 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0))
          * _grad_phi[_j][_qp];
    j -= _gfac * _phi[_j][_qp] * _gravity[_qp];
    j /= _area[_qp];
    j += (1.0 / _area[_qp] - 3.0 * std::pow(_m[_qp] / (_rho[_qp]
          * _area[_qp]), 2.0) * _drho_dh[_qp] / (_rho[_qp] * _area[_qp]))
          * _phi[_j][_qp] * _grad_u[_qp];
    j -= 3.0 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0) * _drho_dp[_qp]
          / (_rho[_qp] * _area[_qp]) * _phi[_j][_qp] * _grad_p[_qp];
    //j -= 3.0 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0) * _drho_dc[_qp]
    //      / (_rho[_qp] * _area[_qp]) * _phi[_j][_qp] * _grad_c[_qp];
    j -= 3.0 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0)
          / (_rho[_qp] * _area[_qp]) * _phi[_j][_qp] *_drho_m_dcvect[_qp];
    j += _dkappa_dm[_qp] * _grad_phi[_j][_qp];
    j += _domega_dm[_qp] * _grad_phi[_j][_qp];
  }

  if (jvar == _p_var_number)
  {
    j -= 3.0 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0) * _drho_dp[_qp]
          * _phi[_j][_qp] / _rho[_qp] * _grad_m[_qp] / _area[_qp];
    j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (_drho_dph[_qp]
          - 3.0 * _drho_dh[_qp] * _drho_dp[_qp] / _rho[_qp]) * _phi[_j][_qp]
          * _grad_u[_qp];
    j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (_drho_dp2[_qp]
          - 3.0 * _drho_dp[_qp] * _drho_dp[_qp] / _rho[_qp]) * _phi[_j][_qp]
          * _grad_p[_qp];
    j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * _drho_dp[_qp]
          * _grad_phi[_j][_qp];
    j += _dkappa_dp[_qp] * _grad_phi[_j][_qp];
    // j += _dkappa_dp2[_qp] * _phi[_j][_qp] * _grad_p[_qp];
    j += _domega_dp[_qp] * _grad_phi[_j][_qp];
  }
/*
  if (jvar == _c_var_number)
  {
    j -= 3.0 * std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 2.0) * _drho_dc[_qp]
          * _phi[_j][_qp] / _rho[_qp] * _grad_m[_qp] / _area[_qp];
    j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (0.0
          - 3.0 * _drho_dh[_qp] * _drho_dc[_qp] / _rho[_qp]) * _phi[_j][_qp]
          * _grad_u[_qp];
    j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * (0.0
          - 3.0 * _drho_dc[_qp] * _drho_dc[_qp] / _rho[_qp]) * _phi[_j][_qp]
          * _grad_c[_qp];
    j -= std::pow(_m[_qp] / (_rho[_qp] * _area[_qp]), 3.0) * _drho_dc[_qp]
          * _grad_phi[_j][_qp];
    j += _dkappa_dp[_qp] * _grad_phi[_j][_qp];
    // j += _dkappa_dp2[_qp] * _phi[_j][_qp] * _grad_p[_qp];
    j += _domega_dc[_qp] * _grad_phi[_j][_qp];
  }
*/
  return j * _test[_i][_qp] * _well_dir[_qp] * _well_sign[_qp];
}
