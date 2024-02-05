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

#include "MoskitoEOS1P.h"

// template <>
// InputParameters
// validParams<MoskitoEOS1P>()
// {
//   InputParameters params = validParams<GeneralUserObject>();

InputParameters
MoskitoEOS1P::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  // MooseEnum Linearity("yes=1 no=0");
  // params.addParam<MooseEnum>(
  //     "cp_linearity", Linearity="yes", "if cp is constant or linear, yes. "
  //     "Otherwise, no.");
  // params.addParam<Real>("deltaT", 20.0, "if cp is nonlinear, an approximation "
  //     "based on , deltaT value, in Celcius is considered for integrations "
  //     "in the case of the enthalpy and temperature conversions");

  return params;
}

MoskitoEOS1P::MoskitoEOS1P(const InputParameters & parameters)
  : GeneralUserObject(parameters)
  // _cp_linear(getParam<MooseEnum>("cp_linearity")),
  // _deltaT(getParam<Real>("deltaT"))
  {}

MoskitoEOS1P::~MoskitoEOS1P() {}

Real
MoskitoEOS1P::h_from_p_T(const Real & pressure, const Real & temperature) const
{
  // if (temperature < _T_ref)
  //   mooseError(name(), ": Temperature should not be less than the reference temperature");

  Real h = 0.0;

  if (_cp_linear)
    h += 0.5 * (cp(pressure, temperature) + cp(pressure, _T_ref)) * (temperature - _T_ref);
  else
  {
    Real n = std::ceil((temperature - _T_ref) / _deltaT);
    Real dT = (temperature - _T_ref) / n;

    for(int i=1;i<n;i++)
      h += cp(pressure, _T_ref + i * dT);

    h += 0.5 * (cp(pressure, _T_ref) + cp(pressure, temperature));
    h *= dT;
  }

  h += _h_ref;

  return h;
}
