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

#pragma once

#include "MoskitoEOS1P.h"
#include "SinglePhaseFluidProperties.h"

// class MoskitoEOS1P_FPModule;
//
// template <>
// InputParameters validParams<MoskitoEOS1P_FPModule>();

class MoskitoEOS1P_FPModule : public MoskitoEOS1P
{
public:
  static InputParameters validParams();
  MoskitoEOS1P_FPModule(const InputParameters & parameters);

  virtual Real h_from_p_T(const Real & pressure, const Real & temperature) const override;
  virtual void h_from_p_T(const Real & pressure, const Real & temperature, Real & h, Real & dh_dp, Real & dh_dT) const;
  virtual Real rho_from_p_T(const Real & pressure, const Real & temperature) const override;
  virtual void rho_from_p_T(const Real & pressure, const Real & temperature,
                        Real & rho, Real & drho_dp, Real & drho_dT) const override;
  // virtual Real mu_from_p_T(const Real & pressure, const Real & temperature) const;
  // virtual void mu_from_p_T(const Real & pressure, const Real & temperature,
  //                       Real & mu, Real & dmu_dp, Real & dmu_dT) const;
  virtual Real cp(const Real & pressure, const Real & temperature) const override;
  virtual Real lambda(const Real & pressure, const Real & temperature) const override;

protected:
  const SinglePhaseFluidProperties & _fp_eos;
};
