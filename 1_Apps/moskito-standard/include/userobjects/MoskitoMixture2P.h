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

#ifndef MOSKITOMIXTURE2P_H
#define MOSKITOMIXTURE2P_H

#include "MoskitoEOS2P.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

// class MoskitoMixture2P;
//
// template <>
// InputParameters validParams<MoskitoMixture2P>();

class MoskitoMixture2P : public MoskitoEOS2P
{
public:
  static InputParameters validParams();
  MoskitoMixture2P(const InputParameters & parameters);

  virtual void VMFrac_T_from_p_h(
      const Real & press, const Real & enthalpy, Real & vmfrac, Real & temp, Real & phase) const override;

  virtual void rho_g_from_p_T(
      const Real & press, const Real & temp, Real & rho, Real & drho_dp, Real & drho_dT, const Real & phase) const override;

  virtual void rho_l_from_p_T(
      const Real & press, const Real & temp, Real & rho, Real & drho_dp, Real & drho_dT, const Real & phase) const override;

  virtual Real rho_g_from_p_T(const Real & press, const Real & temp, const Real & phase) const override;

  virtual Real rho_l_from_p_T(const Real & press, const Real & temp, const Real & phase) const override;

  virtual Real cp_m_from_p_T(
      const Real & press, const Real & temp, const Real & vmfrac, const Real & phase) const override;

  virtual Real rho_m_from_p_h(const Real & press, const Real & enthalpy) const override;

  virtual void PENG_PURE(const Real & press, const Real & temp, Real tc[], Real pc[], Real m[], Real a[], Real b[]) const override;

protected:
  virtual void h_lat(const Real & press, Real & hlat, Real & hsatl, Real & hsatg) const override;

  virtual void PENG_MIX(const Real & press, const Real & temp, Real y[], const Real a[], const Real b[],
   Real & capital_a, Real & capital_b, Real & sum_b, Real & sum_a, Real binary_a[][6], Real kappa[][6]) const override;

  virtual void ROOT_FINDER(const Real & capital_a, const Real & capital_b, Real & z) const override;

  virtual void CONCENTRATION_CALCULATOR(const Real k[], const Real tmf[], Real & non_aqueous_mf, Real x[], Real y[]) const override;

  virtual void FRAC_CONVERTER(const Real omf[], const Real tmf[], const Real & non_aqueous_mf, const Real x[],
     const Real y[], const Real m_mas[], Real s_molal[], Real s_molef[], Real s_massf[], Real molef[], Real massf[],
      Real & non_aqueous_massf, Real s_molality[], Real s_molefrac[], Real & mf_h2o) const override;

  virtual void ACTIVITY(const Real & press, const Real & t, const Real x[],
     const Real tmf[], const Real & non_aqueous_mf, Real gama[], const Real s_molality[]) const override;

  virtual void HENRY(const Real & press, const Real & temp, Real hen[]) const override;

  virtual void FUG(const Real & capital_a, const Real & capital_b, const Real & z, const Real y[], const Real b[],
     const Real & sum_a, const Real & sum_b, const Real binary_a[][6], Real phi[]) const override;

  virtual void K_CONVERGENCE(const Real & press, const Real & temp, const Real hen[], const Real gama[],
     const Real phi[], Real k[], Real k_old[]) const override;

  virtual void TH_DRIESNER(const Real & press, const Real & temp, const Real s_molef[], const Real molef[], Real & th) const override;

  virtual Real SATURATED_H_WATER(const Real & th) const override;

  virtual Real SPECIFIC_V_WATER(const Real & press, const Real & th) const override;

  virtual Real T_EXPANSION_WATER(const Real & th) const override;

  virtual Real P_SAT_WATER(const Real & th) const override;

  virtual void SALT_ENTHALPY(const Real & temp, const Real s_molality[], Real sh[]) const override;

  virtual void LIQUID_ENTHALPY(const Real s_massf[], const Real massf[], const Real sh[],
     Real & liq_h, const Real & water_h, const Real gas_ha[], const Real molarv_deriv[]) const override;

  virtual void GAS_H(const Real & t, const Real y[], const Real & z, const Real & sum_a, const Real & sum_b,
     const Real & capital_b, const Real a[], const Real kappa[][6], Real & gas_h, const Real m[], const Real tc[],
     const Real pc[], const Real & mf_h2o, const Real & hh2o) const override;

  virtual void GAS_RO(const Real & p, const Real & t, const Real & z, const Real y[],
     const Real tc[], const Real pc[], Real & g_ro, const Real & roh2o, const Real & mf_h2o) const override;

  virtual void TV_DRIESNER(const Real & press, const Real & temp, const Real tmf[], Real & tv) const override;

  virtual void LIQUID_RHO(const Real & p, const Real & t, Real & l_rho, const Real s_molality[],
     const Real s_massf[], const Real massf[], const Real molarv_vol[]) const override;

  //virtual void GAS_VIS(const Real & press, const Real & temp, const Real tc[], const Real pc[],
     //const Real y[], const Real & g_ro, Real & g_vis) const override;

  virtual void LIQUID_VIS(const Real & t, const Real s_massf[], const Real massf[], Real & l_vis) const override;

  virtual void VAR_SEQUENCE(const Real & var, vector<Real> & var_seq) const override;

  virtual void VAR_GRADIENT(const Real & press, const Real & temp, const Real tmf[],
     const vector<Real> & var_seq, Real & var_p, Real & var_t, Real & var_z) const override;

  virtual void Transport_PROP(const Real & g_ro, const Real & l_ro, const Real & non_aqueous_mf,
     const Real k[], Real & rog_k_n, Real & rol_n) const override;

  virtual void mix_PROP(const Real x[], const Real y[], const Real & non_aqueous_mf,
     const Real & liq_prop, const Real & gas_prop, const Real m_mas[], Real & mix_prop, Real & mass_frac) const override;

  virtual void Itera(const Real & h, const Real & p, const Real & m,const Real c_vect[], Real & g_ro, Real & l_rho,
     Real & lhx, Real & gas_h, Real & mass_frac, Real & non_aqueous_mf, Real & rho_mix, Real & temper, Real x[],
      Real aux1_arr[], Real aux2_arr[], Real & viscosity, Real & conductivity) const override;
  virtual void water_v(const Real y[], const Real & p, const Real & t, Real & roh2o,
     Real & hh2o, Real & rov, Real & rol, Real & usv, Real & usl, Real & xv) const override;
  virtual void GAS_VIS(const Real & p, const Real & t, const Real tc[], const Real pc[],
     const Real y[], const Real ro[], Real & g_vis, Real vis[], const Real & xv) const override;
  virtual void GAS_IN_LIQ(const Real & t, const Real & p, Real gas_ha[], Real ro[]) const override;
  virtual void CONDUCTIVITY_GAS(const Real & t, const Real & p, Real & g_con, const Real & g_vis,
     Real con[], const Real vis[], const Real y[]) const override;
  virtual void GAS_MV(const Real & t, const Real & p, Real molarv[]) const override;
  virtual void GAS_MV_DERIV(const Real & t, const Real & p, Real molarv_deriv[]) const override;
  virtual void GAS_LD(const Real & t, Real molarv_vol[]) const override;
  virtual void enti(const Real & non_aqueous_massf, const Real &  gas_prop, const Real &  liq_prop, Real & ent) const override;
  virtual void LIQ_CONDUCTIVITY(const Real & t, Real & l_con, const Real s_massf[], const Real massf[]) const override;
  //virtual void LIQUID_VIS(const Real & temp, const Real y[], Real & l_vis) const override;
  //virtual void LIQUID_VIS(const Real & temp, const Real y[], Real & l_vis) const override;
  //virtual void LIQUID_VIS(const Real & temp, const Real y[], Real & l_vis) const override;

  Real t;

};

#endif /* MoskitoMixture2P_H */
