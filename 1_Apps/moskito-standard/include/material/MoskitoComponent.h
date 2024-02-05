/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*  Co-developed by Sebastian Held                                        */
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

#include "Material.h"
#include "MoskitoEOS2P.h"
#include "MoskitoDriftFlux.h"
#include "MoskitoViscosity2P.h"


// class MoskitoComponent;
//
// template <>
// InputParameters validParams<MoskitoComponent>();

class MoskitoComponent : public Material
{
public:
  static InputParameters validParams();
  MoskitoComponent(const InputParameters & parameters);
  virtual void computeQpProperties() override;

  void Iteration(const Real & h, const Real & p, const Real & m,
     const Real c_vect[], Real aux1_arr[], Real aux2_arr[], Real & g_ro, Real & l_rho, Real & lhx, Real & gas_h, Real & mass_frac, Real & non_aqueous_mf, Real & rho_mix,
   Real & u_g, Real & u_l, Real & gamma, Real & kappa, Real & omega, Real & temper);

   void Derivatives_c();

  //Real ResistivityNoAnnulus(const int & begin, const int & end, const bool & hf);
  //Real nonDtimefunction();

  //virtual Real computeReferenceResidual(const Real trial_value, const Real scalar) override;



protected:

  const MoskitoEOS2P & eos_uo;
  const MoskitoDriftFlux & dfm_uo;
  const MoskitoViscosity2P & viscosity_uo;

  const MaterialProperty<Real> & _area;
  const MaterialProperty<Real> & _friction;
  const MaterialProperty<Real> & _dia;
  //const MaterialProperty<Real> & _u_f;
  const MaterialProperty<Real> & _well_sign;
  const MaterialProperty<RealVectorValue> & _gravity;
  const MaterialProperty<RealVectorValue> & _well_dir;

  //MaterialProperty<Real> & _aux1_p;
  //MaterialProperty<Real> & _aux2_p;
  //MaterialProperty<Real> & _u_l_p;
  //MaterialProperty<Real> & _u_g_p;
  //MaterialProperty<Real> & _aux1_h;
  //MaterialProperty<Real> & _aux2_h;
  //MaterialProperty<Real> & _u_l_h;
  //MaterialProperty<Real> & _u_g_h;
  MaterialProperty<Real> & _aux1_c;
  MaterialProperty<Real> & _aux2_c;
  MaterialProperty<Real> & _u_l_c;
  MaterialProperty<Real> & _u_g_c;
  MaterialProperty<Real> & _drho_m_dcc;
  MaterialProperty<Real> & _dgamma_dcc;
  MaterialProperty<Real> & _dkappa_dcc;
  MaterialProperty<Real> & _domega_dcc;
  MaterialProperty<Real> & _T_cc;

  const int _component;
  //const int _ico2;

  const VariableValue & _P;
  const VariableValue & _h;
  const VariableValue & _m;
  const VariableValue & _c;
  const VariableValue & _c2;
  const VariableValue & _c3;
  const VariableValue & _c4;
  const VariableValue & _c5;
  const VariableValue & _c6;
  const VariableValue & _c7;
  const VariableValue & _c8;

  // imported from other materials
  //const MaterialProperty<Real> & _Dti;

  Real _aux1_new[10], _aux2_new[10], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7;
  Real u_g, u_l, gamma, kappa, omega, temper;
  Real c_vect[8];

};
