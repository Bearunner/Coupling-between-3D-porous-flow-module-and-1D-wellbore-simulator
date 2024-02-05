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

// class MoskitoEOS1P_IdealGas;
//
// template <>
// InputParameters validParams<MoskitoEOS1P_IdealGas>();

class MoskitoEOS1P_IdealGas : public MoskitoEOS1P
{
public:
  static InputParameters validParams();
  MoskitoEOS1P_IdealGas(const InputParameters & parameters);

  virtual Real rho_from_p_T(const Real & pressure, const Real & temperature) const override;
  virtual void rho_from_p_T(const Real & pressure, const Real & temperature,
                        Real & rho, Real & drho_dp, Real & drho_dT) const override;
  virtual Real cp(const Real & pressure, const Real & temperature) const override;
  virtual Real lambda(const Real & pressure, const Real & temperature) const override;

protected:
  Real _cp;
  const Real _lambda;
  const Real _molar_mass;
  // Universal Gas constant (J/mol.K)
  const Real _R;
};
