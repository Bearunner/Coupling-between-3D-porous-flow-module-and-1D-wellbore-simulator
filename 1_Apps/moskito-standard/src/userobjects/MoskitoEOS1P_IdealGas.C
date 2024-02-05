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

#include "MoskitoEOS1P_IdealGas.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_IdealGas);

// template <>
// InputParameters
// validParams<MoskitoEOS1P_IdealGas>()
// {
//   InputParameters params = validParams<MoskitoEOS1P>();

InputParameters
MoskitoEOS1P_IdealGas::validParams()
{
  InputParameters params = MoskitoEOS1P::validParams();

  params.addRequiredParam<Real>("molar_mass", "Molar mass of the gas (kg/mol)");
  params.addRequiredParam<Real>("specific_heat",
        "Specific heat capacity at constant volume (J/mol/K)");

  return params;
}

MoskitoEOS1P_IdealGas::MoskitoEOS1P_IdealGas(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _cp(getParam<Real>("specific_heat")),
    _lambda(0.0),
    _molar_mass(getParam<Real>("molar_mass")),
    _R(8.3144598)
{
  _cp += _R;
  _cp /= _molar_mass;
  _T_ref = 0.0;
  _h_ref = 0.0;
}

Real
MoskitoEOS1P_IdealGas::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  return pressure * _molar_mass / (_R * temperature);
}

void
MoskitoEOS1P_IdealGas::rho_from_p_T(const Real & pressure, const Real & temperature,
                            Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho_from_p_T(pressure, temperature);
  drho_dp = _molar_mass / (_R * temperature);
  drho_dT = -pressure * _molar_mass / (_R * temperature * temperature);
}

Real
MoskitoEOS1P_IdealGas::cp(const Real & pressure, const Real & temperature) const
{
  return _cp;
}

Real
MoskitoEOS1P_IdealGas::lambda(const Real & pressure, const Real & temperature) const
{
  return _lambda;
}
