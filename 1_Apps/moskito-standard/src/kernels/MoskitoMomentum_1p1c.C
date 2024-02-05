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

#include "MoskitoMomentum_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoMomentum_1p1c);

InputParameters
MoskitoMomentum_1p1c::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("temperature", "temperature");
  params.addClassDescription("Momentum conservation equation for 2 phase "
        "(liquid and gas) pipe flow and it returns flowrate. A positive "
        "value in BCs is the production scenario and the flow is in the opposite"
        " way of the well direction as defined in the material, and vice versa");
  return params;
}

MoskitoMomentum_1p1c::MoskitoMomentum_1p1c(const InputParameters & parameters)
  : Kernel(parameters),
    _p(coupledValue("pressure")),
    // _grad_p(coupledGradient("pressure")),
    // _grad_T(coupledGradient("temperature")),
    _p_var_number(coupled("pressure")),
    _T_var_number(coupled("temperature")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dT(getMaterialProperty<Real>("drho_dT")),
    _f(getMaterialProperty<Real>("well_moody_friction")),
    _gravity(getMaterialProperty<RealVectorValue>("gravity")),
    _area(getMaterialProperty<Real>("well_area")),
    _perimeter(getMaterialProperty<Real>("well_perimeter")),
    _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector"))
    // _well_sign(getMaterialProperty<Real>("flow_direction_sign"))
{
}

Real
MoskitoMomentum_1p1c::computeQpResidual()
{
  RealVectorValue r = 0.0;

  // r += _drho_dp[_qp] * _grad_p[_qp];
  // r += _drho_dT[_qp] * _grad_T[_qp];
  // r *= _u[_qp] * _u[_qp];
  // r += 2.0 * _rho[_qp] * _u[_qp] * _grad_u[_qp];
  // r += _well_sign[_qp] * _f[_qp] * _rho[_qp] * _u[_qp] * _u[_qp] * _perimeter[_qp]
  //       * _well_dir[_qp] / (8.0 * _area[_qp]);
  // r /= _area[_qp] * _area[_qp];
  // r += _grad_p[_qp];
  // r -= _rho[_qp] * _gravity[_qp];
  // r *= _test[_i][_qp];

  r += -(_p[_qp] + _rho[_qp] * _u[_qp] * _u[_qp] / _area[_qp] / _area[_qp]) * _grad_test[_i][_qp];

  r += -_rho[_qp] * _gravity[_qp] * _test[_i][_qp] + _well_dir[_qp] * _f[_qp] * _rho[_qp] * (_u[_qp] * _u[_qp] / _area[_qp] / _area[_qp]) * _perimeter[_qp] / (8.0 * _area[_qp]) * _test[_i][_qp];

  return r * _well_dir[_qp] ;
}

Real
MoskitoMomentum_1p1c::computeQpJacobian()
{
  RealVectorValue j = 0.0;

  // j += _drho_dp[_qp] * _grad_p[_qp];
  // j += _drho_dT[_qp] * _grad_T[_qp];
  // j *= 2.0 * _phi[_j][_qp] * _u[_qp];
  // j += 2.0 * _rho[_qp] * (_phi[_j][_qp]  * _grad_u[_qp]
  //       + _u[_qp] * _grad_phi[_j][_qp]);
  // j += _well_sign[_qp] * _f[_qp] * _rho[_qp] * _phi[_j][_qp] * _u[_qp] * _perimeter[_qp]
  //       * _well_dir[_qp] / (4.0 * _area[_qp]);
  // j /= _area[_qp] * _area[_qp];
  // j *= _test[_i][_qp];

  j += -(_rho[_qp] * 2 * _u[_qp] * _phi[_j][_qp] / _area[_qp] / _area[_qp]) * _grad_test[_i][_qp] + (_f[_qp] * _rho[_qp] * (2 * _u[_qp] * _phi[_j][_qp] / _area[_qp] / _area[_qp]) * _perimeter[_qp]
        / (8.0 * _area[_qp]) ) * _test[_i][_qp] * _well_dir[_qp];

  // j +=  (_f[_qp] * _rho[_qp] * (2 * _u[_qp] * _phi[_j][_qp] / _area[_qp] / _area[_qp]) * _perimeter[_qp]
  //       / (8.0 * _area[_qp]) ) * _test[_i][_qp];


  return j * _well_dir[_qp] ;
}

Real
MoskitoMomentum_1p1c::computeQpOffDiagJacobian(unsigned int jvar)
{
  RealVectorValue j = 0.0;

  if (jvar == _p_var_number)
  {
    // j += _drho_dp[_qp] * _grad_phi[_j][_qp];
    // j *= _u[_qp] * _u[_qp];
    // j += 2.0 * _drho_dp[_qp] * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp];
    // j += _well_sign[_qp] * _f[_qp] * _drho_dp[_qp] * _phi[_j][_qp] * _u[_qp]
    //       * _u[_qp] * _perimeter[_qp] * _well_dir[_qp] / (8.0 * _area[_qp]);
    // j /= _area[_qp] * _area[_qp];
    // j += _grad_phi[_j][_qp];
    // j -= _drho_dp[_qp] * _phi[_j][_qp] * _gravity[_qp];
    // j *= _test[_i][_qp];

    j += -(_phi[_j][_qp] + _drho_dp[_qp] * _phi[_j][_qp] * _u[_qp] * _u[_qp] / _area[_qp] / _area[_qp]) * _grad_test[_i][_qp];

    j += -_drho_dp[_qp] * _phi[_j][_qp] * _gravity[_qp] * _test[_i][_qp] + _well_dir[_qp] * _f[_qp] * _drho_dp[_qp] * _phi[_j][_qp] * (_u[_qp] * _u[_qp] / _area[_qp] / _area[_qp]) * _perimeter[_qp]
          / (8.0 * _area[_qp]) * _test[_i][_qp];
  }


  if (jvar == _T_var_number)
  {
    // j += _drho_dT[_qp] * _grad_phi[_j][_qp];
    // j *= _u[_qp] * _u[_qp];
    // j += 2.0 * _drho_dT[_qp] * _phi[_j][_qp] * _u[_qp] * _grad_u[_qp];
    // j += _well_sign[_qp] * _f[_qp] * _drho_dT[_qp] * _phi[_j][_qp] * _u[_qp]
    //       * _u[_qp] * _perimeter[_qp] * _well_dir[_qp] / (8.0 * _area[_qp]);
    // j /= _area[_qp] * _area[_qp];
    // j -= _drho_dT[_qp] * _phi[_j][_qp] * _gravity[_qp];
    // j *= _test[_i][_qp];
    j += -(_drho_dT[_qp] * _phi[_j][_qp] * _u[_qp] * _u[_qp] / _area[_qp] / _area[_qp]) * _grad_test[_i][_qp];

    j += -_drho_dT[_qp] * _phi[_j][_qp] * _gravity[_qp] * _test[_i][_qp] + _well_dir[_qp] * _f[_qp] * _drho_dT[_qp] * _phi[_j][_qp] * (_u[_qp] * _u[_qp] / _area[_qp] / _area[_qp]) * _perimeter[_qp]
          / (8.0 * _area[_qp]) * _test[_i][_qp];
    }
    return j * _well_dir[_qp];
}
