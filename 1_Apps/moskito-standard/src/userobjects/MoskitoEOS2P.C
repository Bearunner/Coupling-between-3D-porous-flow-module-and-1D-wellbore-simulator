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

#include "MoskitoEOS2P.h"

// template <>
// InputParameters
// validParams<MoskitoEOS2P>()
// {
//   InputParameters params = validParams<GeneralUserObject>();

InputParameters
MoskitoEOS2P::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addParam<Real>("derivative_tolerance", 0.001,
        "Tolerance to calculate derivatives based on numerical differentiation");

  return params;
}

MoskitoEOS2P::MoskitoEOS2P(const InputParameters & parameters)
  : GeneralUserObject(parameters),
  _tol(getParam<Real>("derivative_tolerance"))
{
}

void
MoskitoEOS2P::rho_m_by_p(const Real & pressure, const Real & enthalpy, Real & rho, Real & drho_dp, Real & drho_dp_2) const
{
  Real rho_plus_tol, rho_minus_tol, tol_p;
  tol_p = _tol * pressure;
  rho_plus_tol = rho_m_from_p_h(pressure + tol_p, enthalpy);
  rho_minus_tol = rho_m_from_p_h(pressure - tol_p, enthalpy);
  rho = rho_m_from_p_h(pressure, enthalpy);

  drho_dp   = (rho_plus_tol - rho_minus_tol) / 2.0 / tol_p;
  drho_dp_2 = (rho_plus_tol - 2.0 * rho + rho_minus_tol) / tol_p / tol_p;
  if(fabs(drho_dp_2)<1.e-11)
    drho_dp_2 = 0.0;
}

void
MoskitoEOS2P::rho_m_by_h(const Real & pressure, const Real & enthalpy, Real & rho, Real & drho_dh, Real & drho_dh_2) const
{
  Real rho_plus_tol, rho_minus_tol, tol_h;
  tol_h = _tol * enthalpy;
  rho_plus_tol = rho_m_from_p_h(pressure, enthalpy + tol_h);
  rho_minus_tol = rho_m_from_p_h(pressure, enthalpy - tol_h);
  rho = rho_m_from_p_h(pressure, enthalpy);

  drho_dh   = (rho_plus_tol - rho_minus_tol) / 2.0 / tol_h;
  drho_dh_2 = (rho_plus_tol - 2.0 * rho + rho_minus_tol) / tol_h / tol_h;
  if(fabs(drho_dh_2)<1.e-11)
    drho_dh_2 = 0.0;
}

void
MoskitoEOS2P::rho_m_by_ph(const Real & pressure, const Real & enthalpy,Real & drho_dph) const
{
  Real dh = _tol * enthalpy;
  Real dp = _tol * pressure;

  drho_dph  = rho_m_from_p_h(pressure + dp,enthalpy + dh) + rho_m_from_p_h(pressure - dp,enthalpy - dh);
  drho_dph -= rho_m_from_p_h(pressure + dp,enthalpy - dh) + rho_m_from_p_h(pressure - dp,enthalpy + dh);
  drho_dph /= 4.0 * dh * dp;
  if(fabs(drho_dph)<1.e-11)
    drho_dph = 0.0;

}
