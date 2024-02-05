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

#include "MoskitoDFGVar.h"

MoskitoDFGVar::MoskitoDFGVar(Real v_m, Real rho_g, Real rho_l, const Real & mfrac, const Real & vfrac,
  Real dia, const Real & dir, const Real & friction, const RealVectorValue & gravity,
  const RealVectorValue & well_dir)
  : _FlowPat(0),
    _v_sg(0.0),
    _v_sl(0.0),
    _C0(0.0),
    _vd(0.0),
    _v_m(fabs(v_m)),
    _rho_g(rho_g),
    _rho_l(rho_l),
    _mfrac(mfrac),
    _vfrac(vfrac),
    _dia(dia),
    _dir(dir),
    _friction(friction),
    _gravity(gravity),
    _well_dir(well_dir),
    _grav(_gravity.norm()),
    _angle(acos(fabs(_gravity * _well_dir / _grav)))
{
}

void
MoskitoDFGVar::DFMOutput(Real & FlowPat, Real & C0, Real & vd)
{
  FlowPat = _FlowPat;
  C0 = _C0;
  vd = _vd;
}
