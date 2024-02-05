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

#include "MoskitoEOS1P_NaturalGas.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_NaturalGas);

// template <>
// InputParameters
// validParams<MoskitoEOS1P_NaturalGas>()
// {
//   InputParameters params = validParams<MoskitoEOS1P>();

InputParameters
MoskitoEOS1P_NaturalGas::validParams()
{
  InputParameters params = MoskitoEOS1P::validParams();

  params.addRequiredParam<Real>("molar_mass", "Molar mass of the gas (kg/mol)");
  params.addParam<Real>("specific_gravity", 1.0,
        "Specific gravity (air = 1.0)");
  params.addParam<Real>("specific_heat", 1.0e3,
        "Constant specific heat capacity at constant pressure (J/kg/K)");

  return params;
}

MoskitoEOS1P_NaturalGas::MoskitoEOS1P_NaturalGas(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _molar_mass(getParam<Real>("molar_mass")),
    _gamma_g(getParam<Real>("specific_gravity")),
    _R(8.3144598),
    _cp(getParam<Real>("specific_heat")),
    _lambda(0.0)
{
  Pseudo_Critical_Calc(_gamma_g);
  // should be checked later for a nonlinear cp
  _T_ref = 0.0;
  _h_ref = 0.0;
}

Real
MoskitoEOS1P_NaturalGas::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  Real z = z_factor(pressure, temperature);

  return pressure * _molar_mass / (z * _R * temperature);
}

void
MoskitoEOS1P_NaturalGas::rho_from_p_T(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
  Real z = z_factor(pressure, temperature);
  rho = this->rho_from_p_T(pressure, temperature);

  Real dz_dp, dz_dT, h;

  h = 0.0001 * pressure;
  dz_dp = (z_factor(pressure + h, temperature) - z_factor(pressure - h, temperature)) / (2.0 * h);
  drho_dp  = _molar_mass / (z * _R * temperature);
  drho_dp *= 1.0 - dz_dp * pressure / z;

  h = 0.0001 * temperature;
  dz_dT = (z_factor(pressure, temperature + h) - z_factor(pressure, temperature - h)) / (2.0 * h);
  drho_dT  = - rho * (dz_dT / z + 1.0 / temperature);
}

Real
MoskitoEOS1P_NaturalGas::cp(const Real & pressure, const Real & temperature) const
{
  return _cp;
}

Real
MoskitoEOS1P_NaturalGas::lambda(const Real & pressure, const Real & temperature) const
{
  return _lambda;
}

void
MoskitoEOS1P_NaturalGas::Pseudo_Critical_Calc(const Real & g)
{
  // Rankine scale
  _T_pc  = 120.1 + 425.0 * g - 62.9 * g * g;
  // Kelvin
  _T_pc *= 5.0 / 9.0;
  // PSI
  _P_pc  = 671.1 + 14.0 * g - 34.3 * g * g;
  // pascal
  _P_pc *= 6894.7572931783;
}

Real
MoskitoEOS1P_NaturalGas::z_factor(const Real & pressure, const Real & temperature) const
{
  Real T_pr, P_pr, t, z, y, A, B, C, D, E, F, G;

  T_pr = temperature / _T_pc;
  P_pr = pressure / _P_pc;
  t = 1.0 / T_pr;

  A = a[1] * t * exp(a[2] * pow(1.0 - t, 2.0)) * P_pr;
  B = a[3] * t + a[4] * t * t + a[5] * pow(t * P_pr, 6.0);
  C = a[9] + a[8] * t * P_pr + a[7] * pow(t * P_pr, 2.0) + a[6] * pow(t * P_pr, 3.0);
  D = a[10] * t * exp(a[11] * pow(1.0 - t, 2.0));
  E = a[12] * t + a[13] * t * t + a[14] * t * t * t;
  F = a[15] * t + a[16] * t * t + a[17] * t * t * t;
  G = a[18] + a[19] * t;

  y = D * P_pr / ((1.0 + A * A) / C - (A * A * B) / pow(C, 3.0));

  z  = D * P_pr * (1.0 + y + y * y - y * y * y);
  z /= D * P_pr + E * y * y - F * pow(y, G);
  z /= pow(1.0 - y, 3.0);

  return z;
}
