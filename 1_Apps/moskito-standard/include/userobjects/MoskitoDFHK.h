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

#ifndef MOSKITODFHK_H
#define MOSKITODFHK_H

#include "MoskitoDriftFlux.h"

// class MoskitoDFHK;
//
// template <>
// InputParameters validParams<MoskitoDFHK>();

// forward declaration
 class MoskitoHKLVar;

class MoskitoDFHK : public MoskitoDriftFlux
{
public:
  static InputParameters validParams();
  MoskitoDFHK(const InputParameters & parameters);

  virtual void DFMCalculator(MoskitoDFGVar & input) const override;

  // Preparation of HK individual parameters
  void HKinitialisation(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  // final calculator
  void HKcalculator(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  // Calculation of Void fraction
  // void HKvfrac(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;

protected:
  // Main code determination of flow pattern, C0, vd and void fraction
  void cal_v_s(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;

  // Equations for calcualation of bubbly / Taylor rise velocities
  void cal_vd_b(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void cal_vd_tb(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void cal_vd_mix(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;

  // Calculation of thresholds of transitions
  void cal_v_gb(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void cal_v_gc(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void cal_v_ms(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;

  // Determine drift flow parameters for each flow pattern
  void Det_db_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void Det_bubbly_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void Det_churn_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void Det_annular_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;
  void Det_slug_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const;

  // Interpolation between different C0
  Real interpol(const Real & C0_1, const Real & C0_2, const Real & v_denom, const Real & v_num) const;

private:
  Real _surf_ten;
  const Real _C0b = 1.2;
  const Real _C0db = 1.2;
  const Real _C0s_u = 1.2;
  const Real _C0s_d = 1.12;
  const Real _C0c_u = 1.15;
  const Real _C0c_d = 1.12;
  const Real _C0a = 1.0;
};

class MoskitoHKLVar
{
public:
  // Superficial gas velocity in m/s
  Real v_sg= 0.0;
  // Superficial fluid velocity in m/s
  Real v_sl= 0.0;
  // velocity of dispersed bubbles
  Real vd_b= 0.0;
  // velocity of Taylor bubbles
  Real vd_tb= 0.0;
  // mixture velocity between Taylor bubbles velocity and dispersed bubble velocity for slug and churn flow
  Real vd_mix= 0.0;
  //  Thereshold for transition from bubbly to slug flow
  Real v_gb= 0.0;
  //  Thereshold for transition from churn to annular flow
  Real v_gc= 0.0;
  // Thereshold for transition from bubbly / slug flow to d_dubbly and slug to churn flow
  Real v_ms= 0.0;
};

#endif /* MOSKITODFHK_H */
