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

#include "MoskitoFluidWellGeneral.h"
#include "MoskitoEOS2P.h"
#include "MoskitoViscosity2P.h"
#include "MoskitoDriftFlux.h"

// class MoskitoFluidWell_2p1c;
//
// template <>
// InputParameters validParams<MoskitoFluidWell_2p1c>();

class MoskitoFluidWell_2p1c : public MoskitoFluidWellGeneral
{
public:
  static InputParameters validParams();
  MoskitoFluidWell_2p1c(const InputParameters & parameters);
  virtual void computeQpProperties() override;
  void PhaseVelocities();
  void GammaDerivatives();
  void KappaDerivatives();
  void OmegaDerivatives();
  Real gamma(const Real & m, const Real & g_ro, const Real & l_rho, const Real & lhx,
     const Real & gas_h, const Real & mass_frac, const Real & non_aqueous_mf, const Real & rho_mix);
  Real kappa(const Real & m, const Real & g_ro, const Real & l_rho, const Real & lhx,
     const Real & gas_h, const Real & mass_frac, const Real & non_aqueous_mf, const Real & rho_mix);
  Real omega(const Real & m, const Real & g_ro, const Real & l_rho, const Real & lhx,
     const Real & gas_h, const Real & mass_frac, const Real & non_aqueous_mf, const Real & rho_mix);
  void Iteration(const Real & h, const Real & p, const Real & m,
     const Real c_vect[], Real aux1_arr[], Real aux2_arr[], Real & g_ro, Real & l_rho, Real & lhx, Real & gas_h, Real & mass_frac, Real & non_aqueous_mf, Real & rho_mix);
  void Derivatives_mat();

protected:
  // Userobject to equation of state
  const MoskitoEOS2P & eos_uo;
  // Userobject to Viscosity Eq
  const MoskitoViscosity2P & viscosity_uo;
  // Userobject to Drift Fluc model
  const MoskitoDriftFlux & dfm_uo;

  // temperature
  MaterialProperty<Real> & _T;
  // The specific heat of mixture at constant pressure
  MaterialProperty<Real> & _cp_m;
  // Density of gas
  MaterialProperty<Real> & _rho_g;
  // Density of liquid
  MaterialProperty<Real> & _rho_l;
  // Density of mixture
  MaterialProperty<Real> & _rho_m;
  // Profile-adjusted density of mixture
  MaterialProperty<Real> & _rho_pam;
  // The first derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_m_dp;
  // The second derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_m_dp2;
  // The first derivative of mixture density wrt enthalpy
  MaterialProperty<Real> & _drho_m_dh;
  // The second derivative of mixture density wrt enthalpy
  MaterialProperty<Real> & _drho_m_dh2;
  // The second derivative of mixture density wrt enthalpy and pressure
  MaterialProperty<Real> & _drho_m_dph;

  //MaterialProperty<Real> & _drho_m_dc;
  // mass_fraction
  MaterialProperty<Real> & _vmfrac;
  // Gas velocity
  MaterialProperty<Real> & _u_g;
  MaterialProperty<Real> & _u_g_p_plus;
  MaterialProperty<Real> & _u_g_p_minus;
  MaterialProperty<Real> & _u_g_h_plus;
  MaterialProperty<Real> & _u_g_h_minus;
  //MaterialProperty<Real> & _u_g_c_plus;
  //MaterialProperty<Real> & _u_g_c_minus;
  MaterialProperty<Real> & _u_g_m_plus;
  MaterialProperty<Real> & _u_g_m_minus;
  // Liquid velocity
  MaterialProperty<Real> & _u_l;
  MaterialProperty<Real> & _u_l_p_plus;
  MaterialProperty<Real> & _u_l_p_minus;
  MaterialProperty<Real> & _u_l_h_plus;
  MaterialProperty<Real> & _u_l_h_minus;
  //MaterialProperty<Real> & _u_l_c_plus;
  //MaterialProperty<Real> & _u_l_c_minus;
  MaterialProperty<Real> & _u_l_m_plus;
  MaterialProperty<Real> & _u_l_m_minus;
  MaterialProperty<Real> & _T_h_plus;
  MaterialProperty<Real> & _T_h_minus;
  MaterialProperty<Real> & _T_h;
  MaterialProperty<Real> & _T_p_plus;
  MaterialProperty<Real> & _T_p_minus;
  MaterialProperty<Real> & _T_p;
  //MaterialProperty<Real> & _T_c_plus;
  //MaterialProperty<Real> & _T_c_minus;
  //MaterialProperty<Real> & _T_c;

  // void_fraction
  MaterialProperty<Real> & _vfrac;
  // current phase
  MaterialProperty<Real> & _phase;
  // drift velocity
  MaterialProperty<Real> & _u_d;
  // flow type parameter
  MaterialProperty<Real> & _c0;
  // flow pattern
  MaterialProperty<Real> & _flow_pat;

  // The gamma derivatives
  MaterialProperty<Real> & _dgamma_dp;
  MaterialProperty<Real> & _dgamma_dh;
  MaterialProperty<Real> & _dgamma_dm;
  //MaterialProperty<Real> & _dgamma_dc;

  // The kappa first derivatives
  MaterialProperty<Real> & _dkappa_dp;
  // The kappa first derivatives
  MaterialProperty<Real> & _dkappa_dh;
  // The kappa first derivatives
  MaterialProperty<Real> & _dkappa_dm;
  //MaterialProperty<Real> & _dkappa_dc;
  // The kappa second derivatives
  MaterialProperty<Real> & _dkappa_dph;
  // The kappa second derivatives
  MaterialProperty<Real> & _dkappa_dpm;
  // The kappa second derivatives
  MaterialProperty<Real> & _dkappa_dhm;
  // The kappa second derivatives
  MaterialProperty<Real> & _dkappa_dp2;
  // The kappa second derivatives
  MaterialProperty<Real> & _dkappa_dh2;
  // The kappa second derivatives
  MaterialProperty<Real> & _dkappa_dm2;

  // The omega derivatives
  MaterialProperty<Real> & _domega_dp;
  MaterialProperty<Real> & _domega_dh;
  MaterialProperty<Real> & _domega_dm;
  //MaterialProperty<Real> & _domega_dc;
  MaterialProperty<Real> & _aux1;
  MaterialProperty<Real> & _aux2;
  MaterialProperty<Real> & _aux1_CH4;
  MaterialProperty<Real> & _aux2_CH4;
  MaterialProperty<Real> & _aux1_N2;
  MaterialProperty<Real> & _aux2_N2;
  MaterialProperty<Real> & _aux1_H2S;
  MaterialProperty<Real> & _aux2_H2S;
  MaterialProperty<Real> & _aux1_NaCl;
  MaterialProperty<Real> & _aux2_NaCl;
  MaterialProperty<Real> & _aux1_KCl;
  MaterialProperty<Real> & _aux2_KCl;
  MaterialProperty<Real> & _aux1_CaCl2;
  MaterialProperty<Real> & _aux2_CaCl2;
  MaterialProperty<Real> & _aux1_MgCl2;
  MaterialProperty<Real> & _aux2_MgCl2;
  MaterialProperty<Real> & _aux1_p;
  MaterialProperty<Real> & _aux2_p;
  MaterialProperty<Real> & _aux1_p_CH4;
  MaterialProperty<Real> & _aux2_p_CH4;
  MaterialProperty<Real> & _aux1_p_N2;
  MaterialProperty<Real> & _aux2_p_N2;
  MaterialProperty<Real> & _aux1_p_H2S;
  MaterialProperty<Real> & _aux2_p_H2S;
  MaterialProperty<Real> & _aux1_p_NaCl;
  MaterialProperty<Real> & _aux2_p_NaCl;
  MaterialProperty<Real> & _aux1_p_KCl;
  MaterialProperty<Real> & _aux2_p_KCl;
  MaterialProperty<Real> & _aux1_p_CaCl2;
  MaterialProperty<Real> & _aux2_p_CaCl2;
  MaterialProperty<Real> & _aux1_p_MgCl2;
  MaterialProperty<Real> & _aux2_p_MgCl2;
  MaterialProperty<Real> & _u_l_p;
  MaterialProperty<Real> & _u_g_p;
  MaterialProperty<Real> & _aux1_h;
  MaterialProperty<Real> & _aux2_h;
  MaterialProperty<Real> & _aux1_h_CH4;
  MaterialProperty<Real> & _aux2_h_CH4;
  MaterialProperty<Real> & _aux1_h_N2;
  MaterialProperty<Real> & _aux2_h_N2;
  MaterialProperty<Real> & _aux1_h_H2S;
  MaterialProperty<Real> & _aux2_h_H2S;
  MaterialProperty<Real> & _aux1_h_NaCl;
  MaterialProperty<Real> & _aux2_h_NaCl;
  MaterialProperty<Real> & _aux1_h_KCl;
  MaterialProperty<Real> & _aux2_h_KCl;
  MaterialProperty<Real> & _aux1_h_CaCl2;
  MaterialProperty<Real> & _aux2_h_CaCl2;
  MaterialProperty<Real> & _aux1_h_MgCl2;
  MaterialProperty<Real> & _aux2_h_MgCl2;
  MaterialProperty<Real> & _u_l_h;
  MaterialProperty<Real> & _u_g_h;
  //MaterialProperty<Real> & _aux1_c;
  //MaterialProperty<Real> & _aux2_c;
  //MaterialProperty<Real> & _u_l_c;
  //MaterialProperty<Real> & _u_g_c;
  MaterialProperty<Real> & _u_l_m;
  MaterialProperty<Real> & _u_g_m;
  MaterialProperty<Real> & _x3;
  MaterialProperty<Real> & _x2;
  MaterialProperty<Real> & _zero_0;
  MaterialProperty<Real> & _viscosity;
  MaterialProperty<Real> & _lambda;

  //const MaterialProperty<Real> & _drho_dc1;
  //const MaterialProperty<Real> & _drho_dc2;
  //const MaterialProperty<Real> & _dgamma_dc1;
  //const MaterialProperty<Real> & _dgamma_dc2;
  //const MaterialProperty<Real> & _dkappa_dc1;
  //const MaterialProperty<Real> & _dkappa_dc2;
  //const MaterialProperty<Real> & _domega_dc1;
  //const MaterialProperty<Real> & _domega_dc2;
  //const MaterialProperty<Real> & _T_c1;
  //const MaterialProperty<Real> & _T_c2;

  //MaterialProperty<Real> & _drho_m_dc;
  //MaterialProperty<Real> & _dgamma_dc;
  //MaterialProperty<Real> & _dkappa_dc;
  //MaterialProperty<Real> & _domega_dc;
  //MaterialProperty<Real> & _T_c;


  // The coupled enthalpy
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

  // The gradient of the coupled variables
  const VariableGradient & _grad_c;
  const VariableGradient & _grad_c2;
  const VariableGradient & _grad_m;
  const VariableGradient & _grad_h;
  const VariableGradient & _grad_p;

  //Real aux1;
  //Real aux2;

  Real shimba;
  Real dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9;
  Real c_vect[8], aux1_arr[10], aux2_arr[10];
};
