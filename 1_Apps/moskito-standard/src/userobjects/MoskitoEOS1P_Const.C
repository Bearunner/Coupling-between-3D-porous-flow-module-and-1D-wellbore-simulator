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

#include "MoskitoEOS1P_Const.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_Const);

// template <>
// InputParameters
// validParams<MoskitoEOS1P_Const>()
// {
//   InputParameters params = validParams<MoskitoEOS1P>();
InputParameters
MoskitoEOS1P_Const::validParams()
{
  InputParameters params = MoskitoEOS1P::validParams();
  
  params.addParam<Real>("specific_heat", 4200,
        "Constant specific heat at constant pressure (J/kg.K)");
  params.addParam<Real>("thermal_conductivity", 0.6,
        "Constant thermal conductivity (W/m/K)");
  params.addParam<Real>("density", 1000.0, "Constant density");

  return params;
}

MoskitoEOS1P_Const::MoskitoEOS1P_Const(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _cp(getParam<Real>("specific_heat")),
    _lambda(getParam<Real>("thermal_conductivity")),
    _density(getParam<Real>("density"))
{
}

Real
MoskitoEOS1P_Const::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  return _density;
}

void
MoskitoEOS1P_Const::rho_from_p_T(const Real & pressure, const Real & temperature,
                              Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho_from_p_T(pressure, temperature);
  drho_dp = 0.0;
  drho_dT = 0.0;
}

Real
MoskitoEOS1P_Const::cp(const Real & pressure, const Real & temperature) const
{
  return _cp;
}

Real
MoskitoEOS1P_Const::lambda(const Real & pressure, const Real & temperature) const
{
  return _lambda;
}
