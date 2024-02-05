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

#include "MoskitoViscosityConst.h"

registerMooseObject("MoskitoApp", MoskitoViscosityConst);

// template <>
// InputParameters
// validParams<MoskitoViscosityConst>()
// {
//   InputParameters params = validParams<MoskitoViscosity1P>();

InputParameters
MoskitoViscosityConst::validParams()
{
  InputParameters params = MoskitoViscosity1P::validParams();

  params.addParam<Real>("viscosity", 0.001, "Constant viscosity value (Pa.s)");

  return params;
}

MoskitoViscosityConst::MoskitoViscosityConst(const InputParameters & parameters)
  : MoskitoViscosity1P(parameters), _mu0(getParam<Real>("viscosity"))
{
}

Real MoskitoViscosityConst::mu(Real /* pressure */, Real /* temperature */) const { return _mu0; }

void
MoskitoViscosityConst::dmu_dpT(
    Real pressure, Real temperature, Real & mu, Real & dmu_dp, Real & dmu_dT) const
{
  mu = this->mu(pressure, temperature);
  dmu_dp = 0.0;
  dmu_dT = 0.0;
}
