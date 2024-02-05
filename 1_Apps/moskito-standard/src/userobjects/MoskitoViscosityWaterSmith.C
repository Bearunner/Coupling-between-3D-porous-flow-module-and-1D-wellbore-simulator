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

#include "MoskitoViscosityWaterSmith.h"

registerMooseObject("MoskitoApp", MoskitoViscosityWaterSmith);

// template <>
// InputParameters
// validParams<MoskitoViscosityWaterSmith>()
// {
//   InputParameters params = validParams<MoskitoViscosity1P>();

InputParameters
MoskitoViscosityWaterSmith::validParams()
{
  InputParameters params = MoskitoViscosity1P::validParams();

  /* The viscosity is based on Smith and Chapmann 1983 */

  return params;
}

MoskitoViscosityWaterSmith::MoskitoViscosityWaterSmith(const InputParameters & parameters)
  : MoskitoViscosity1P(parameters)
{
}

Real
MoskitoViscosityWaterSmith::mu(Real /* pressure */, Real temperature) const
{
  return 2.4e-5 * std::pow(10.0, 248.37 / (temperature - 140));
}

void
MoskitoViscosityWaterSmith::dmu_dpT(
    Real pressure, Real temperature, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  mu = this->mu(pressure, temperature);
  dmu_dp = 0.0;
  dmu_dT = mu * -571.893 / std::pow((temperature - 140), 2.0);
}
