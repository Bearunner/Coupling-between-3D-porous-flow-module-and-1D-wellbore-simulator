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

#include "MoskitoEOS1P_SimpleFluidProperties.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_SimpleFluidProperties);

// template <>
// InputParameters
// validParams<MoskitoEOS1P_SimpleFluidProperties>()
// {
//   InputParameters params = validParams<MoskitoEOS1P>();

InputParameters
MoskitoEOS1P_SimpleFluidProperties::validParams()
{
  InputParameters params = MoskitoEOS1P::validParams();

  params.addParam<Real>("thermal_expansion", 0,
        "Constant coefficient of thermal expansion (1/K)");
  params.addParam<Real>("reference_density", 1000,
        "Density at the reference pressure and temperature (kg/m^3)");
  params.addParam<Real>("specific_heat", 4186,
        "Constant specific heat at constant pressure (J/kg.K)");
  params.addParam<Real>("thermal_conductivity", 0.6,
        "Constant thermal conductivity (W/m/K)");
  params.addRangeCheckedParam<Real>("bulk_modulus", 2E9, "bulk_modulus>0",
        "Constant bulk modulus (Pa)");

  return params;
}

MoskitoEOS1P_SimpleFluidProperties::MoskitoEOS1P_SimpleFluidProperties(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _rho_ref(getParam<Real>("reference_density")),
    _cp(getParam<Real>("specific_heat")),
    _lambda(getParam<Real>("thermal_conductivity")),
    _bulk_modulus(getParam<Real>("bulk_modulus"))
{
}

Real
MoskitoEOS1P_SimpleFluidProperties::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  return _rho_ref * std::exp(pressure / _bulk_modulus - _thermal_expansion * temperature);
}

void
MoskitoEOS1P_SimpleFluidProperties::rho_from_p_T(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho_from_p_T(pressure, temperature);
  drho_dp = rho / _bulk_modulus;
  drho_dT = -_thermal_expansion * rho;
}

Real
MoskitoEOS1P_SimpleFluidProperties::cp(const Real & pressure, const Real & temperature) const
{
  return _cp;
}

Real
MoskitoEOS1P_SimpleFluidProperties::lambda(const Real & pressure, const Real & temperature) const
{
  return _lambda;
}
