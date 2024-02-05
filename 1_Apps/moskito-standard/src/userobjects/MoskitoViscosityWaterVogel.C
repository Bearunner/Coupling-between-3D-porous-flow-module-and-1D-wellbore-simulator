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

#include "MoskitoViscosityWaterVogel.h"

registerMooseObject("MoskitoApp", MoskitoViscosityWaterVogel);

// template <>
// InputParameters
// validParams<MoskitoViscosityWaterVogel>()
// {
//   InputParameters params = validParams<MoskitoViscosity1P>();

InputParameters
MoskitoViscosityWaterVogel::validParams()
{
  InputParameters params = MoskitoViscosity1P::validParams();

  /* The viscosity is based on Vogel viscosity equation */

  return params;
}

MoskitoViscosityWaterVogel::MoskitoViscosityWaterVogel(const InputParameters & parameters)
  : MoskitoViscosity1P(parameters)
{
}

Real
MoskitoViscosityWaterVogel::mu(Real /* pressure */, Real temperature) const
{
  return 1e-3*exp(-3.7188+578.919/(-137.546+temperature));
}

void
MoskitoViscosityWaterVogel::dmu_dpT(
    Real pressure, Real temperature, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  mu = this->mu(pressure, temperature);
  dmu_dp = 0.0;
  dmu_dT = mu*-578.919/std::pow(-137.546+temperature,2.0);
}
