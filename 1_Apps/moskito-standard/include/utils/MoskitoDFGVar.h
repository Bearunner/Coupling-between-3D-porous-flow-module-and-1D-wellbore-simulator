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

#ifndef MOSKITODFGVAR_H
#define MOSKITODFGVAR_H

#include "Material.h"

// it is for storing global variables for drift flux uo
class MoskitoDFGVar
{
public:
  MoskitoDFGVar(Real v_m, Real rho_g, Real rho_l, const Real & mfrac, const Real & vfrac, Real dia,
    const Real & dir, const Real & friction, const RealVectorValue & gravity,
    const RealVectorValue & well_dir);

  void DFMOutput(Real & FlowPat, Real & C0, Real & vd);

  // Flow pattern 0 = nothing, 1 = bubbly, 2 = dispersed_bubbly, 3 = slug, 4 = churn, 5 = annular
  Real _FlowPat;
  // superficial velocities
  Real _v_sg;
  Real _v_sl;
  // Drift Flux parameters
  Real _C0;
  Real _vd;

  // mixture velocity of 2P-system
  Real _v_m;
  // Gas density
  Real _rho_g;
  // Liquid density
  Real _rho_l;
  // Mass fraction of void phase
  const Real _mfrac;
  // Volume fraction of void phase
  const Real _vfrac;
  // Well diameter
  Real _dia;
  // Flow direction
  Real _dir;
  // Well friction
  const Real _friction;
  // The gravity acceleration as a vector
  const RealVectorValue _gravity;
  // unit vector along well
  const RealVectorValue _well_dir;
  // gravity acceleration value
  Real _grav;
  // Angle between gravity vector and well_unity_vector
  Real _angle;
};

#endif /* MOSKITODFGVAR_H */
