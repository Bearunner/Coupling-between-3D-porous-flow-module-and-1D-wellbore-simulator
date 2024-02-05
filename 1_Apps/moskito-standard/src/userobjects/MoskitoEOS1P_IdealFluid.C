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

#include "MoskitoEOS1P_IdealFluid.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_IdealFluid);

// template <>
// InputParameters
// validParams<MoskitoEOS1P_IdealFluid>()
// {
//   InputParameters params = validParams<MoskitoEOS1P>();

InputParameters
MoskitoEOS1P_IdealFluid::validParams()
{
  InputParameters params = MoskitoEOS1P::validParams();

  params.addParam<Real>("thermal_expansion_0", 4.0E-4,
        "Constant coefficient of thermal expansion (1/K)");
  params.addParam<Real>("thermal_expansion_1", 0.0,
        "Constant coefficient of thermal expansion (1/K^2)");
  params.addParam<Real>("reference_density", 998.29,
        "Density at the reference pressure and temperature (kg/m^3)");
  params.addParam<Real>("reference_temperature", 293.15,
        "Reference temperature (K)");
  params.addParam<Real>("reference_pressure", 101325,
        "Reference pressure (Pa)");
  params.addParam<Real>("reference_enthalpy", 83950,
        "Specific enthalpy (J/kg)");
  params.addParam<Real>("specific_heat", 4200,
        "Constant specific heat at constant pressure (J/kg.K)");
  params.addParam<Real>("thermal_conductivity", 0.6,
        "Constant thermal conductivity (W/m/K)");
  params.addRangeCheckedParam<Real>("bulk_modulus", 2.15E9, "bulk_modulus>0",
        "Constant bulk modulus (Pa)");

  return params;
}

MoskitoEOS1P_IdealFluid::MoskitoEOS1P_IdealFluid(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _rho_ref(getParam<Real>("reference_density")),
    _P_ref(getParam<Real>("reference_pressure")),
    _cp(getParam<Real>("specific_heat")),
    _lambda(getParam<Real>("thermal_conductivity")),
    _thermal_expansion_0(getParam<Real>("thermal_expansion_0")),
    _thermal_expansion_1(getParam<Real>("thermal_expansion_1")),
    _bulk_modulus(getParam<Real>("bulk_modulus"))
{
  _T_ref = getParam<Real>("reference_temperature");
  _h_ref = getParam<Real>("reference_enthalpy");
}

Real
MoskitoEOS1P_IdealFluid::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  return _rho_ref * std::exp((pressure-_P_ref) / _bulk_modulus - _thermal_expansion_0 *
          (temperature - _T_ref) - 0.5 * _thermal_expansion_1 * (temperature * temperature - _T_ref * _T_ref));
}

void
MoskitoEOS1P_IdealFluid::rho_from_p_T(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho_from_p_T(pressure, temperature);
  drho_dp = rho / _bulk_modulus;
  drho_dT = (-_thermal_expansion_0 - _thermal_expansion_1 * temperature) * rho;
}

Real
MoskitoEOS1P_IdealFluid::cp(const Real & pressure, const Real & temperature) const
{
  return _cp;
}

Real
MoskitoEOS1P_IdealFluid::lambda(const Real & pressure, const Real & temperature) const
{
  return _lambda;
}
