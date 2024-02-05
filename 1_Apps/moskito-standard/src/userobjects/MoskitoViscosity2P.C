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

#include "MoskitoViscosity2P.h"

registerMooseObject("MoskitoApp", MoskitoViscosity2P);

// template <>
// InputParameters
// validParams<MoskitoViscosity2P>()
// {
//   InputParameters params = validParams<GeneralUserObject>();

InputParameters
MoskitoViscosity2P::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<UserObjectName>("ve_uo_gas",
        "The name of the viscosity equation userobject for gas");
  params.addRequiredParam<UserObjectName>("ve_uo_liquid",
        "The name of the viscosity equation userobject for liquid");
  MooseEnum Mixture
        ("Series Parallel ME1 ME2 EMT Mean_ME12",
        "Mean_ME12");
  params.addParam<MooseEnum>("mixing_type", Mixture,
        "Type of mixing to calculate the mixture viscosity [Series, Parallel, "
        "Maxwell Eucken 1, Maxwell Eucken 2, Effective Medium Theory and "
        "Mean Maxwell Eucken 1 & 2");
  return params;
}

MoskitoViscosity2P::MoskitoViscosity2P(const InputParameters & parameters)
  : GeneralUserObject(parameters),
  gas(getUserObject<MoskitoViscosity1P>("ve_uo_gas")),
  liquid(getUserObject<MoskitoViscosity1P>("ve_uo_liquid")),
  _mt(getParam<MooseEnum>("mixing_type"))
{
}

MoskitoViscosity2P::~MoskitoViscosity2P() {}

Real
MoskitoViscosity2P::mixture_mu(Real pressure, Real temperature, Real mass_fraction) const
{
  Real mu = 0.0;
  Real & x = mass_fraction;
  Real mu_l = liquid.mu(pressure, temperature);
  Real mu_g = gas.mu(pressure, temperature);

  switch (_mt)
  {
    case MT::Series:
      mu = (1.0 - x) / mu_l + x / mu_g;
      mu = 1.0 / mu;
    break;
    case MT::Parallel:
      mu = (1.0 - x) * mu_l + x * mu_g;
      break;
    case MT::ME1:
      mu = ME1_calc(mu_l, mu_g, x);
      break;
    case MT::ME2:
      mu = ME2_calc(mu_l, mu_g, x);
      break;
    case MT::EMT:
      mu = EMT_calc(mu_l, mu_g, x);
      break;
    case MT::Mean_ME12:
      mu = 0.5 * (ME1_calc(mu_l, mu_g, x) + ME2_calc(mu_l, mu_g, x));
      break;
  }

  return mu;
}

Real
MoskitoViscosity2P::ME1_calc(Real mu_l, Real mu_g, Real x) const
{
  Real mu;
  mu  = 2.0 * mu_l + mu_g - 2.0 * (mu_l - mu_g) * x;
  mu /= 2.0 * mu_l + mu_g + (mu_l - mu_g) * x;
  mu *= mu_l;
  return mu;
}

Real
MoskitoViscosity2P::ME2_calc(Real mu_l, Real mu_g, Real x) const
{
  Real mu;
  mu  = 2.0 * mu_g + mu_l - 2.0 * (mu_g - mu_l) * (1.0 - x);
  mu /= 2.0 * mu_g + mu_l + (mu_g - mu_l) * (1.0 - x);
  mu *= mu_g;
  return mu;
}

Real
MoskitoViscosity2P::EMT_calc(Real mu_l, Real mu_g, Real x) const
{
  Real mu;
  mu  = (3.0 * x - 1.0) * mu_g;
  mu += mu_l * (3.0 * (1.0 - x) - 1.0);
  mu += std::sqrt(std::pow((3.0 * x - 1.0) * mu_g + mu_l *
        (3.0 * (1.0 - x) - 1.0),2.0) + 8.0 * mu_l * mu_g);
  mu *= 0.25;
  return mu;
}
