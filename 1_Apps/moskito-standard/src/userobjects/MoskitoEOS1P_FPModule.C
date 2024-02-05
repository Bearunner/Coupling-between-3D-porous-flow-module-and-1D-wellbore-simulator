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

#include "MoskitoEOS1P_FPModule.h"

registerMooseObject("MoskitoApp", MoskitoEOS1P_FPModule);

// template <>
// InputParameters
// validParams<MoskitoEOS1P_FPModule>()
// {
//   InputParameters params = validParams<MoskitoEOS1P>();

InputParameters
MoskitoEOS1P_FPModule::validParams()
{
  InputParameters params = MoskitoEOS1P::validParams();

  params.addRequiredParam<UserObjectName>("SinglePhase_fp",
          "The name of the FluidProperties UserObject");

  return params;
}

MoskitoEOS1P_FPModule::MoskitoEOS1P_FPModule(const InputParameters & parameters)
  : MoskitoEOS1P(parameters),
    _fp_eos(getUserObject<SinglePhaseFluidProperties>("SinglePhase_fp"))
{
}

Real
MoskitoEOS1P_FPModule::h_from_p_T(const Real & pressure, const Real & temperature) const
{
  return _fp_eos.h_from_p_T(pressure, temperature);
}

void
MoskitoEOS1P_FPModule::h_from_p_T(const Real & pressure, const Real & temperature,
                            Real & h, Real & dh_dp, Real & dh_dT) const
{
  _fp_eos.h_from_p_T(pressure, temperature, h, dh_dp, dh_dT);
}

Real
MoskitoEOS1P_FPModule::rho_from_p_T(const Real & pressure, const Real & temperature) const
{
  return _fp_eos.rho_from_p_T(pressure, temperature);
}

void
MoskitoEOS1P_FPModule::rho_from_p_T(const Real & pressure, const Real & temperature,
                            Real & rho, Real & drho_dp, Real & drho_dT) const
{
  _fp_eos.rho_from_p_T(pressure, temperature, rho, drho_dp, drho_dT);
}

// Real
// MoskitoEOS1P_FPModule::mu_from_p_T(const Real & pressure, const Real & temperature) const
// {
//   return _fp_eos.mu_from_p_T(pressure, temperature);
// }
//
// void
// MoskitoEOS1P_FPModule::mu_from_p_T(const Real & pressure, const Real & temperature,
//                             Real & mu, Real & dmu_dp, Real & dmu_dT) const
// {
//   _fp_eos.mu_from_p_T(pressure, temperature, mu, dmu_dp, dmu_dT);
// }

Real
MoskitoEOS1P_FPModule::cp(const Real & pressure, const Real & temperature) const
{
  return _fp_eos.cp_from_p_T(pressure, temperature);
}

Real
MoskitoEOS1P_FPModule::lambda(const Real & pressure, const Real & temperature) const
{
  return _fp_eos.k_from_p_T(pressure, temperature);
}
