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

#include "MoskitoFluidWell_2p1c.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell_2p1c);

// template <>
// InputParameters
// validParams<MoskitoFluidWell_2p1c>()
// {
//   InputParameters params = validParams<MoskitoFluidWellGeneral>();

InputParameters
MoskitoFluidWell_2p1c::validParams()
{
  InputParameters params = MoskitoFluidWellGeneral::validParams();

  params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable (J/kg)");
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for 2 phase EOS");
  params.addRequiredParam<UserObjectName>("viscosity_uo",
        "The name of the userobject for 2 phase viscosity Eq");
  params.addRequiredParam<UserObjectName>("drift_flux_uo",
        "The name of the userobject for drift flux model");
  params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable (m^3/s)");
  params.addRequiredCoupledVar("concentration", "concentration");
  params.addCoupledVar("c2", "Coupled value c2");
  params.addCoupledVar("c3", "Coupled value c3");
  params.addCoupledVar("c4", "Coupled value c4");
  params.addCoupledVar("c5", "Coupled value c5");
  params.addCoupledVar("c6", "Coupled value c6");
  params.addCoupledVar("c7", "Coupled value c7");
  params.addCoupledVar("c8", "Coupled value c8");
  return params;
}

MoskitoFluidWell_2p1c::MoskitoFluidWell_2p1c(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS2P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity2P>("viscosity_uo")),
    dfm_uo(getUserObject<MoskitoDriftFlux>("drift_flux_uo")),
    _T(declareProperty<Real>("temperature")),
    _cp_m(declareProperty<Real>("specific_heat")),
    _rho_g(declareProperty<Real>("gas_density")),
    _rho_l(declareProperty<Real>("liquid_density")),
    _rho_m(declareProperty<Real>("density")),
    _rho_pam(declareProperty<Real>("profile_mixture_density")),
    _drho_m_dp(declareProperty<Real>("drho_dp")),
    _drho_m_dp2(declareProperty<Real>("drho_dp2")),
    _drho_m_dh(declareProperty<Real>("drho_dh")),
    _drho_m_dh2(declareProperty<Real>("drho_dh2")),
    _drho_m_dph(declareProperty<Real>("drho_dph")),
    //_drho_m_dc(declareProperty<Real>("drho_dc")),
    _vmfrac(declareProperty<Real>("mass_fraction")),
    _u_g(declareProperty<Real>("gas_velocity")),
    _u_g_p_plus(declareProperty<Real>("u_g_p_plus")),
    _u_g_p_minus(declareProperty<Real>("u_g_p_minus")),
    _u_g_h_plus(declareProperty<Real>("u_g_h_plus")),
    _u_g_h_minus(declareProperty<Real>("u_g_h_minus")),
    //_u_g_c_plus(declareProperty<Real>("u_g_c_plus")),
    //_u_g_c_minus(declareProperty<Real>("u_g_c_minus")),
    _u_g_m_plus(declareProperty<Real>("u_g_m_plus")),
    _u_g_m_minus(declareProperty<Real>("u_g_m_minus")),
    _u_l(declareProperty<Real>("liquid_velocity")),
    _u_l_p_plus(declareProperty<Real>("u_l_p_plus")),
    _u_l_p_minus(declareProperty<Real>("u_l_p_minus")),
    _u_l_h_plus(declareProperty<Real>("u_l_h_plus")),
    _u_l_h_minus(declareProperty<Real>("u_l_h_minus")),
    //_u_l_c_plus(declareProperty<Real>("u_l_c_plus")),
    //_u_l_c_minus(declareProperty<Real>("u_l_c_minus")),
    _u_l_m_plus(declareProperty<Real>("u_l_m_plus")),
    _u_l_m_minus(declareProperty<Real>("u_l_m_minus")),
    _T_h_plus(declareProperty<Real>("T_h_plus")),
    _T_h_minus(declareProperty<Real>("T_h_minus")),
    _T_h(declareProperty<Real>("T_h")),
    _T_p_plus(declareProperty<Real>("T_p_plus")),
    _T_p_minus(declareProperty<Real>("T_p_minus")),
    _T_p(declareProperty<Real>("T_p")),
    //_T_c_plus(declareProperty<Real>("T_c_plus")),
    //_T_c_minus(declareProperty<Real>("T_c_minus")),
    //_T_c(declareProperty<Real>("T_c")),
    _vfrac(declareProperty<Real>("void_fraction")),
    _phase(declareProperty<Real>("current_phase")),
    _u_d(declareProperty<Real>("drift_velocity")),
    _c0(declareProperty<Real>("flow_type_c0")),
    _flow_pat(declareProperty<Real>("flow_pattern")),
    _dgamma_dp(declareProperty<Real>("dgamma_dp")),
    _dgamma_dh(declareProperty<Real>("dgamma_dh")),
    _dgamma_dm(declareProperty<Real>("dgamma_dm")),
    //_dgamma_dc(declareProperty<Real>("dgamma_dc")),
    _dkappa_dp(declareProperty<Real>("dkappa_dp")),
    _dkappa_dh(declareProperty<Real>("dkappa_dh")),
    _dkappa_dm(declareProperty<Real>("dkappa_dm")),
    //_dkappa_dc(declareProperty<Real>("dkappa_dc")),
    _dkappa_dph(declareProperty<Real>("dkappa_dph")),
    _dkappa_dpm(declareProperty<Real>("dkappa_dpm")),
    _dkappa_dhm(declareProperty<Real>("dkappa_dhm")),
    _dkappa_dp2(declareProperty<Real>("dkappa_dp2")),
    _dkappa_dh2(declareProperty<Real>("dkappa_dh2")),
    _dkappa_dm2(declareProperty<Real>("dkappa_dm2")),
    _domega_dp(declareProperty<Real>("domega_dp")),
    _domega_dh(declareProperty<Real>("domega_dh")),
    _domega_dm(declareProperty<Real>("domega_dm")),
    //_domega_dc(declareProperty<Real>("domega_dc")),
    _aux1(declareProperty<Real>("aux1")),
    _aux2(declareProperty<Real>("aux2")),
    _aux1_CH4(declareProperty<Real>("aux1_CH4")),
    _aux2_CH4(declareProperty<Real>("aux2_CH4")),
    _aux1_N2(declareProperty<Real>("aux1_N2")),
    _aux2_N2(declareProperty<Real>("aux2_N2")),
    _aux1_H2S(declareProperty<Real>("aux1_H2S")),
    _aux2_H2S(declareProperty<Real>("aux2_H2S")),
    _aux1_NaCl(declareProperty<Real>("aux1_NaCl")),
    _aux2_NaCl(declareProperty<Real>("aux2_NaCl")),
    _aux1_KCl(declareProperty<Real>("aux1_KCl")),
    _aux2_KCl(declareProperty<Real>("aux2_KCl")),
    _aux1_CaCl2(declareProperty<Real>("aux1_CaCl2")),
    _aux2_CaCl2(declareProperty<Real>("aux2_CaCl2")),
    _aux1_MgCl2(declareProperty<Real>("aux1_MgCl2")),
    _aux2_MgCl2(declareProperty<Real>("aux2_MgCl2")),
    _aux1_p(declareProperty<Real>("aux1_p")),
    _aux2_p(declareProperty<Real>("aux2_p")),
    _aux1_p_CH4(declareProperty<Real>("aux1_p_CH4")),
    _aux2_p_CH4(declareProperty<Real>("aux2_p_CH4")),
    _aux1_p_N2(declareProperty<Real>("aux1_p_N2")),
    _aux2_p_N2(declareProperty<Real>("aux2_p_N2")),
    _aux1_p_H2S(declareProperty<Real>("aux1_p_H2S")),
    _aux2_p_H2S(declareProperty<Real>("aux2_p_H2S")),
    _aux1_p_NaCl(declareProperty<Real>("aux1_p_NaCl")),
    _aux2_p_NaCl(declareProperty<Real>("aux2_p_NaCl")),
    _aux1_p_KCl(declareProperty<Real>("aux1_p_KCl")),
    _aux2_p_KCl(declareProperty<Real>("aux2_p_KCl")),
    _aux1_p_CaCl2(declareProperty<Real>("aux1_p_CaCl2")),
    _aux2_p_CaCl2(declareProperty<Real>("aux2_p_CaCl2")),
    _aux1_p_MgCl2(declareProperty<Real>("aux1_p_MgCl2")),
    _aux2_p_MgCl2(declareProperty<Real>("aux2_p_MgCl2")),
    _u_l_p(declareProperty<Real>("u_l_p")),
    _u_g_p(declareProperty<Real>("u_g_p")),
    _aux1_h(declareProperty<Real>("aux1_h")),
    _aux2_h(declareProperty<Real>("aux2_h")),
    _aux1_h_CH4(declareProperty<Real>("aux1_h_CH4")),
    _aux2_h_CH4(declareProperty<Real>("aux2_h_CH4")),
    _aux1_h_N2(declareProperty<Real>("aux1_h_N2")),
    _aux2_h_N2(declareProperty<Real>("aux2_h_N2")),
    _aux1_h_H2S(declareProperty<Real>("aux1_h_H2S")),
    _aux2_h_H2S(declareProperty<Real>("aux2_h_H2S")),
    _aux1_h_NaCl(declareProperty<Real>("aux1_h_NaCl")),
    _aux2_h_NaCl(declareProperty<Real>("aux2_h_NaCl")),
    _aux1_h_KCl(declareProperty<Real>("aux1_h_KCl")),
    _aux2_h_KCl(declareProperty<Real>("aux2_h_KCl")),
    _aux1_h_CaCl2(declareProperty<Real>("aux1_h_CaCl2")),
    _aux2_h_CaCl2(declareProperty<Real>("aux2_h_CaCl2")),
    _aux1_h_MgCl2(declareProperty<Real>("aux1_h_MgCl2")),
    _aux2_h_MgCl2(declareProperty<Real>("aux2_h_MgCl2")),
    _u_l_h(declareProperty<Real>("u_l_h")),
    _u_g_h(declareProperty<Real>("u_g_h")),
    //_aux1_c(declareProperty<Real>("aux1_c")),
    //_aux2_c(declareProperty<Real>("aux2_c")),
    //_u_l_c(declareProperty<Real>("u_l_c")),
    //_u_g_c(declareProperty<Real>("u_g_c")),
    _u_l_m(declareProperty<Real>("u_l_m")),
    _u_g_m(declareProperty<Real>("u_g_m")),
    _x3(declareProperty<Real>("x3")),
    _x2(declareProperty<Real>("x2")),
    _zero_0(declareProperty<Real>("zero_0")),
    _viscosity(declareProperty<Real>("viscosity")),
    _lambda(declareProperty<Real>("thermal_conductivity")),
    //_drho_dc1(getMaterialProperty<Real>("m5")),
    //_drho_dc2(getMaterialProperty<Real>("s5")),
    //_dgamma_dc1(getMaterialProperty<Real>("m6")),
    //_dgamma_dc2(getMaterialProperty<Real>("s6")),
    //_dkappa_dc1(getMaterialProperty<Real>("m7")),
    //_dkappa_dc2(getMaterialProperty<Real>("s7")),
    //_domega_dc1(getMaterialProperty<Real>("m8")),
    //_domega_dc2(getMaterialProperty<Real>("s8")),
    //_T_c1(getMaterialProperty<Real>("m9")),
    //_T_c2(getMaterialProperty<Real>("s9")),

    //_drho_m_dc(declareProperty<Real>("drho_dc")),
    //_dgamma_dc(declareProperty<Real>("dgamma_dc")),
    //_dkappa_dc(declareProperty<Real>("dkappa_dc")),
    //_domega_dc(declareProperty<Real>("domega_dc")),
    //_T_c(declareProperty<Real>("T_c")),

    _h(coupledValue("enthalpy")),
    _m(coupledValue("massrate")),
    _c(coupledValue("concentration")),
    //_c2(coupledValue("c2")),
    _c2(isCoupled("c2") ? coupledValue("c2") : _zero),
    _c3(isCoupled("c3") ? coupledValue("c3") : _zero),
    _c4(isCoupled("c4") ? coupledValue("c4") : _zero),
    _c5(isCoupled("c5") ? coupledValue("c5") : _zero),
    _c6(isCoupled("c6") ? coupledValue("c6") : _zero),
    _c7(isCoupled("c7") ? coupledValue("c7") : _zero),
    _c8(isCoupled("c8") ? coupledValue("c8") : _zero),
    _grad_c(coupledGradient("concentration")),
    _grad_c2(coupledGradient("c2")),
    _grad_m(coupledGradient("massrate")),
    _grad_h(coupledGradient("enthalpy")),
    _grad_p(coupledGradient("pressure"))
{
}

void
MoskitoFluidWell_2p1c::computeQpProperties()
{
  MoskitoFluidWellGeneral::computeQpProperties();

  _zero_0[_qp] = 0.0;

  c_vect[0] = _c[_qp];
  c_vect[1] = _c2[_qp];
  c_vect[2] = _c3[_qp];
  c_vect[3] = _c4[_qp];
  c_vect[4] = _c5[_qp];
  c_vect[5] = _c6[_qp];
  c_vect[6] = _c7[_qp];
  c_vect[7] = _c8[_qp];

  //cout<<"well2p c_vect[0] = "<<c_vect[0]<<" c_vect[1] = "<<c_vect[1]<<" c_vect[2] = "<<c_vect[2]<<endl;
/*
  for(int kk=0;kk<=7;kk++)
  {
    if(c_vect[kk] <= 1e-100)
    {
      c_vect[kk] = 0.0;
    }
  }
*/
  //cout<<"well2p"<<c_vect[1]<<endl;
  //cout<<"c_vect[2] in material = "<<c_vect[2]<<endl;
  //cout<<"_c3[_qp] = "<<_c3[_qp]<<endl;

  _phase[_qp] = 2.0;

  //_drho_m_dp[_qp] = 0.00000000001;
  _drho_m_dp2[_qp] = 0.00000000001;
  //_drho_m_dh[_qp] = -0.00000000001;
  _drho_m_dh2[_qp] = -0.00000000001;
  _drho_m_dph[_qp] = -0.00000000001;

  //Iteration(-8344.44, 10000000.0, 10, 0.3, _aux1[_qp], _aux2[_qp]);

  Derivatives_mat();

  Iteration(_h[_qp], _P[_qp], _m[_qp], c_vect, aux1_arr, aux2_arr, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);

  //std::cout<<_u_g[_qp]<<','<<_u_l[_qp]<<std::endl;

  //GammaDerivatives();
  //KappaDerivatives();
  //OmegaDerivatives();

}
void
MoskitoFluidWell_2p1c::PhaseVelocities()
{
  // based on mass weighted flow rate
  // momentum eq is valid only by mass mixing flow rate
  if (_u[_qp] != 0.0)
  {
    _u_g[_qp]  = _c0[_qp] * _rho_m[_qp] * _u[_qp] + _rho_l[_qp] * _u_d[_qp];
    _u_g[_qp] /= _rho_pam[_qp];
    _u_l[_qp]  = (1.0 - _vfrac[_qp] * _c0[_qp]) * _rho_m[_qp]  * _u[_qp] - _rho_g[_qp] * _vfrac[_qp] * _u_d[_qp];
    _u_l[_qp] /= (1.0 - _vfrac[_qp]) * _rho_pam[_qp];
  }

}

Real
MoskitoFluidWell_2p1c::gamma(const Real & m, const Real & g_ro, const Real & l_rho, const Real & lhx,
   const Real & gas_h, const Real & mass_frac, const Real & non_aqueous_mf, const Real & rho_mix)
{
  Real vmfrac, gamma = 0.0;

  if(m != 0.0)
  {
    Real vfrac, rho_l, rho_g, rho_m, rho_pam;
    rho_l = l_rho;
    rho_g = g_ro;
    vmfrac = mass_frac;
    rho_m = rho_mix;
    vfrac = non_aqueous_mf;
    rho_pam = rho_g * _c0[_qp]  * vfrac + (1.0 - vfrac * _c0[_qp]) * rho_l;

    gamma  = vfrac / (1.0 - vfrac);
    gamma *= rho_g * rho_l * rho_m / (rho_pam * rho_pam);
    gamma *= std::pow((_c0[_qp] - 1.0) * (m / rho_m / _area[_qp]) + _u_d[_qp] , 2.0);
  }

  return gamma;
}

Real
MoskitoFluidWell_2p1c::kappa(const Real & m, const Real & g_ro, const Real & l_rho, const Real & lhx,
   const Real & gas_h, const Real & mass_frac, const Real & non_aqueous_mf, const Real & rho_mix)
{
  Real vmfrac, kappa = 0.0;

  if(m != 0.0)
  {
    Real vfrac, rho_l, rho_g, rho_m, rho_pam, h_g, h_l;
    rho_l = l_rho;
    rho_g = g_ro;
    vmfrac = mass_frac;
    rho_m = rho_l * rho_g / (vmfrac * (rho_l - rho_g) + rho_g);
    vfrac = non_aqueous_mf;

    Real dummy, c0, u_d;
    MoskitoDFGVar DFinp(m / rho_m / _area[_qp], rho_g, rho_l, vmfrac, vfrac,
        _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
    dfm_uo.DFMCalculator(DFinp);
    DFinp.DFMOutput(dummy, c0, u_d);

    rho_pam = rho_g * c0 * vfrac + (1.0 - vfrac * c0) * rho_l;
    h_g = gas_h;
    h_l = lhx;

    kappa  = vfrac * rho_g * rho_l / rho_pam * (h_g - h_l);
    kappa *= (c0 - 1.0) * (m / rho_m / _area[_qp]) + u_d;
  }

  return kappa;
}

Real
MoskitoFluidWell_2p1c::omega(const Real & m, const Real & g_ro, const Real & l_rho, const Real & lhx,
   const Real & gas_h, const Real & mass_frac, const Real & non_aqueous_mf, const Real & rho_mix)
{
  Real vmfrac, omega = 0.0;

  if(m != 0.0)
  {
    Real vfrac, rho_l, rho_g, rho_m, u_g, u_l, rho_pam, dummy, c0, u_d;
    rho_l = l_rho;
    rho_g = g_ro;
    vmfrac = mass_frac;
    rho_m = rho_mix;
    vfrac = non_aqueous_mf;
    rho_pam = rho_g * c0  * vfrac + (1.0 - vfrac * c0) * rho_l;

    Real v = m / rho_m / _area[_qp] ;

    u_g  = (c0 * rho_m * v + rho_l * u_d) / rho_pam;
    u_l  = (1.0 - vfrac * c0) * rho_m * v - rho_g * vfrac * u_d;
    u_l /= (1.0 - vfrac) * rho_pam;

    omega -= 3.0 * u_g * u_l * v;
    omega += std::pow(u_g,3.0) * (1.0 + rho_g * vfrac /rho_m);
    omega += std::pow(u_l,3.0) * (1.0 + rho_l * (1.0 - vfrac) / rho_m);
    omega *= 0.5 * vfrac * (1.0 - vfrac) * rho_g * rho_l / rho_m;
  }

  return omega;
}
/*
void
MoskitoFluidWell_2p1c::GammaDerivatives()
{
  _dgamma_dp[_qp]  = 0.0; _dgamma_dh[_qp]  = 0.0; _dgamma_dm[_qp]  = 0.0;

  if (_phase[_qp] == 2.0)
  {
    Real dh, dm, dp;
    Real a,b,c;
    Real tol = 1.0e-3;
    dh = tol * _h[_qp]; dp = tol * _P[_qp]; dm = tol * _m[_qp];

    if (dh != 0.0)
    {
      _dgamma_dh[_qp]  = gamma(_h[_qp] + dh, _P[_qp], _m[_qp]) - gamma(_h[_qp] - dh, _P[_qp], _m[_qp]);
      _dgamma_dh[_qp] /= 2.0 * dh;
    }

    if (dp != 0.0)
    {
    _dgamma_dp[_qp]  = gamma(_h[_qp], _P[_qp] + dp, _m[_qp]) - gamma(_h[_qp], _P[_qp] - dp, _m[_qp]);
    _dgamma_dp[_qp] /= 2.0 * dp;
    }

    if (dm != 0.0)
    {
    _dgamma_dm[_qp]  = gamma(_h[_qp], _P[_qp], _m[_qp] + dm) - gamma(_h[_qp], _P[_qp], _m[_qp] - dm);
    _dgamma_dm[_qp] /= 2.0 * dm;
    }
  }
}

void
MoskitoFluidWell_2p1c::KappaDerivatives()
{
  _dkappa_dp[_qp] = 0.0; _dkappa_dh[_qp] = 0.0; _dkappa_dm[_qp] = 0.0;
  _dkappa_dph[_qp] = 0.0; _dkappa_dpm[_qp] = 0.0; _dkappa_dhm[_qp] = 0.0;
  _dkappa_dp2[_qp] = 0.0; _dkappa_dh2[_qp] = 0.0; _dkappa_dm2[_qp] = 0.0;

  if (_phase[_qp] == 2.0)
  {
    Real dh, dm, dp;
    Real tol = 1.0e-3;
    dh = tol * _h[_qp]; dp = tol * _P[_qp]; dm = tol * _m[_qp];

    if (dh != 0.0)
    {
      _dkappa_dh[_qp]  = kappa(_h[_qp] + dh, _P[_qp], _m[_qp]) - kappa(_h[_qp] - dh, _P[_qp], _m[_qp]);
      _dkappa_dh[_qp] /= 2.0 * dh;
    }

    if (dp != 0.0)
    {
    _dkappa_dp[_qp]  = kappa(_h[_qp], _P[_qp] + dp, _m[_qp]) - kappa(_h[_qp], _P[_qp] - dp, _m[_qp]);
    _dkappa_dp[_qp] /= 2.0 * dp;
    }

    if (dm != 0.0)
    {
    _dkappa_dm[_qp]  = kappa(_h[_qp], _P[_qp], _m[_qp] + dm) - kappa(_h[_qp], _P[_qp], _m[_qp] - dm);
    _dkappa_dm[_qp] /= 2.0 * dm;
    }
    if (dp * dh != 0.0)
    {
    _dkappa_dph[_qp]  = kappa(_h[_qp] + dh, _P[_qp] + dp, _m[_qp]) + kappa(_h[_qp] - dh, _P[_qp] - dp, _m[_qp]);
    _dkappa_dph[_qp] -= kappa(_h[_qp] + dh, _P[_qp] - dp, _m[_qp]) + kappa(_h[_qp] - dh, _P[_qp] + dp, _m[_qp]);
    _dkappa_dph[_qp] /= 4.0 * dh * dp;
    }

    if (dh * dm != 0.0)
    {
    _dkappa_dhm[_qp]  = kappa(_h[_qp] + dh, _P[_qp], _m[_qp] + dm) + kappa(_h[_qp] - dh, _P[_qp], _m[_qp] - dm);
    _dkappa_dhm[_qp] -= kappa(_h[_qp] + dh, _P[_qp], _m[_qp] - dm) + kappa(_h[_qp] - dh, _P[_qp], _m[_qp] + dm);
    _dkappa_dhm[_qp] /= 4.0 * dh * dm;
    }

    if (dp * dm != 0.0)
    {
    _dkappa_dpm[_qp]  = kappa(_h[_qp], _P[_qp] + dp, _m[_qp] + dm) + kappa(_h[_qp], _P[_qp] - dp, _m[_qp] - dm);
    _dkappa_dpm[_qp] -= kappa(_h[_qp], _P[_qp] + dp, _m[_qp] - dm) + kappa(_h[_qp], _P[_qp] - dp, _m[_qp] + dm);
    _dkappa_dpm[_qp] /= 4.0 * dp * dm;
    }

    if (dp != 0.0)
    {
    _dkappa_dp2[_qp]  = kappa(_h[_qp], _P[_qp] + dp, _m[_qp]) + kappa(_h[_qp], _P[_qp] - dp, _m[_qp]);
    _dkappa_dp2[_qp] -= 2.0 * kappa(_h[_qp], _P[_qp], _m[_qp]);
    _dkappa_dp2[_qp] /=  dp * dp;
    }
    if (dh != 0.0)
    {
    _dkappa_dh2[_qp]  = kappa(_h[_qp] + dh, _P[_qp], _m[_qp]) + kappa(_h[_qp] - dh, _P[_qp], _m[_qp]);
    _dkappa_dh2[_qp] -= 2.0 * kappa(_h[_qp], _P[_qp], _m[_qp]);
    _dkappa_dh2[_qp] /=  dh * dh;
    }
    if (dm != 0.0)
    {
    _dkappa_dm2[_qp]  = kappa(_h[_qp], _P[_qp], _m[_qp] + dm) + kappa(_h[_qp], _P[_qp], _m[_qp] - dm);
    _dkappa_dm2[_qp] -= 2.0 * kappa(_h[_qp], _P[_qp], _m[_qp]);
    _dkappa_dm2[_qp] /=  dm * dm;
    }
  }
}

void
MoskitoFluidWell_2p1c::OmegaDerivatives()
{
  _domega_dp[_qp] = 0.0; _domega_dh[_qp] = 0.0; _domega_dm[_qp] = 0.0;

  if (_phase[_qp] == 2.0)
  {
    Real dh, dm, dp;
    Real tol = 1.0e-3;
    dh = tol * _h[_qp]; dp = tol * _P[_qp]; dm = tol * _m[_qp];

    if (dh != 0.0)
    {
      _domega_dh[_qp]  = omega(_h[_qp] + dh, _P[_qp], _m[_qp]) - omega(_h[_qp] - dh, _P[_qp], _m[_qp]);
      _domega_dh[_qp] /= 2.0 * dh;
    }

    if (dp != 0.0)
    {
    _domega_dp[_qp]  = omega(_h[_qp], _P[_qp] + dp, _m[_qp]) - omega(_h[_qp], _P[_qp] - dp, _m[_qp]);
    _domega_dp[_qp] /= 2.0 * dp;
    }

    if (dm != 0.0)
    {
    _domega_dm[_qp]  = omega(_h[_qp], _P[_qp], _m[_qp] + dm) - omega(_h[_qp], _P[_qp], _m[_qp] - dm);
    _domega_dm[_qp] /= 2.0 * dm;
    }
  }
}
*/
void
MoskitoFluidWell_2p1c::Iteration(const Real & h, const Real & p, const Real & m,
   const Real c_vect[], Real aux1_arr[], Real aux2_arr[], Real & g_ro, Real & l_rho, Real & lhx, Real & gas_h, Real & mass_frac, Real & non_aqueous_mf, Real & rho_mix)
{
  //std::cout<<_P[_qp]<<','<<_h[_qp]<<','<<_m[_qp]<<','<<_c[_qp]<<','<<std::endl;
Real temper, x[6], vis, con;

    eos_uo.Itera(h, p, m, c_vect, g_ro, l_rho, lhx, gas_h, mass_frac, non_aqueous_mf, rho_mix, temper, x, aux1_arr, aux2_arr, vis, con);

    _aux1[_qp] = aux1_arr[3];
    _aux2[_qp] = aux2_arr[3];

    _aux1_CH4[_qp] = aux1_arr[2];
    _aux2_CH4[_qp] = aux2_arr[2];

    _aux1_N2[_qp] = aux1_arr[4];
    _aux2_N2[_qp] = aux2_arr[4];

    _aux1_H2S[_qp] = aux1_arr[5];
    _aux2_H2S[_qp] = aux2_arr[5];

    _aux1_NaCl[_qp] = aux1_arr[6];
    _aux2_NaCl[_qp] = aux2_arr[6];

    _aux1_KCl[_qp] = aux1_arr[7];
    _aux2_KCl[_qp] = aux2_arr[7];

    _aux1_CaCl2[_qp] = aux1_arr[8];
    _aux2_CaCl2[_qp] = aux2_arr[8];

    _aux1_MgCl2[_qp] = aux1_arr[9];
    _aux2_MgCl2[_qp] = aux2_arr[9];


/*
  std::cout<<non_aqueous_mf<<','<<x[1]<<','<<x[3]<<','<<y[1]<<','<<y[3]<<std::endl;
  std::cout<<mm[1]<<','<<mm[3]<<','<<a[1]<<','<<a[3]<<','<<b[1]<<','<<b[3]<<std::endl;
  std::cout<<capital_a<<','<<capital_b<<','<<sum_a<<','<<sum_b<<std::endl;
  std::cout<<z<<std::endl;
  std::cout<<hen[1]<<','<<hen[2]<<','<<hen[3]<<','<<hen[4]<<std::endl;
  std::cout<<phi[1]<<','<<phi[2]<<','<<phi[3]<<','<<phi[4]<<std::endl;
  std::cout<<k[1]<<','<<k[2]<<','<<k[3]<<','<<k[4]<<','<<k[5]<<std::endl;
  std::cout<<th<<std::endl;
  std::cout<<water_h<<std::endl;
  std::cout<<sh[1]<<','<<sh[2]<<','<<sh[3]<<std::endl;
  std::cout<<liq_h<<std::endl;
  std::cout<<h_depart<<std::endl;
  std::cout<<g_ro<<std::endl;
  std::cout<<tv<<std::endl;
  std::cout<<l_rho<<std::endl;
  std::cout<<g_vis<<std::endl;
  std::cout<<l_vis<<std::endl;
  std::cout<<aux1<<','<<aux2<<std::endl;
*/
  _x3[_qp] = x[3];
  _x2[_qp] = x[2];
  _vmfrac[_qp] = mass_frac;
  _T[_qp] = temper;
  _phase[_qp] = 2;
  _rho_m[_qp] = rho_mix;
  _rho_l[_qp] = l_rho;
  _rho_g[_qp] = g_ro;
  _vfrac[_qp]  = non_aqueous_mf;
  _cp_m[_qp]  = 4000.0;

  _u[_qp] = fabs(m) / _rho_m[_qp] / _area[_qp];
  //_Re[_qp] = _rho_m[_qp] * _dia[_qp] * _u[_qp] / viscosity_uo.mixture_mu(_P[_qp], _T[_qp], _vmfrac[_qp]);
  _Re[_qp] = _rho_m[_qp] * _dia[_qp] * _u[_qp] / vis;
  _viscosity[_qp] = vis;
  _lambda[_qp] = con;
  //std::cout<<m<<','<<_rho_m[_qp]<<','<<_area[_qp]<<','<<_u[_qp]<<std::endl;
  //std::cout<<_dia[_qp]<<','<<viscosity_uo.mixture_mu(_P[_qp], _T[_qp], _vmfrac[_qp])<<','<<_Re[_qp]<<std::endl<<std::endl;
  if (_f_defined)
    _friction[_qp] = _u_f;
  else
    MoodyFrictionFactor(_friction[_qp], _rel_roughness, _Re[_qp], _roughness_type);

  // drift-flux calculator section
    MoskitoDFGVar DFinp(_u[_qp], _rho_g[_qp], _rho_l[_qp], _vmfrac[_qp], _vfrac[_qp],
      _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
    dfm_uo.DFMCalculator(DFinp);
    DFinp.DFMOutput(_flow_pat[_qp], _c0[_qp], _u_d[_qp]);

  _rho_pam[_qp] = _rho_g[_qp] * _c0[_qp]  * _vfrac[_qp] + (1.0 - _vfrac[_qp] * _c0[_qp]) * _rho_l[_qp];

  //std::cout<<_friction[_qp]<<','<<_c0[_qp]<<','<<_u_d[_qp]<<std::endl;

  //std::cout<<_c0[_qp]<<','<<_u_d[_qp]<<std::endl;

  PhaseVelocities();


  //std::cout<<_u_g[_qp]<<','<<_u_l[_qp]<<','<<_u[_qp]<<std::endl;

  //std::cout<<_u_l[_qp]<<','<<_u_g[_qp]<<','<<_u_d[_qp]<<','<<_u[_qp]<<std::endl;

  //std::cout<<t<<','<<it<<std::endl;
  //std::cout<<aux1<<','<<aux2<<std::endl;
  //cout<<"well2p   "<<g_ro<<" "<<l_rho<<" "<<lhx<<" "<<gas_h<<" "<<mass_frac<<" "<<non_aqueous_mf<<" "<<rho_mix<<" "<<temper<<" "<<aux1_arr[4]<<" "<<aux2_arr[4]<<" "<<endl;

}

void
MoskitoFluidWell_2p1c::Derivatives_mat()
{

  Real dh, dm, dp, dc, aux1_p_plus[10], aux2_p_plus[10], aux1_p_minus[10], aux2_p_minus[10];
  Real aux1_h_plus[10], aux2_h_plus[10], aux1_h_minus[10], aux2_h_minus[10];
  Real aux1_c_plus[10], aux2_c_plus[10], aux1_c_minus[10], aux2_c_minus[10];
  Real aux1_m_plus[10], aux2_m_plus[10], aux1_m_minus[10], aux2_m_minus[10];
  Real dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8, dummy9;
  Real gamma_p_plus, gamma_p_minus, gamma_h_plus, gamma_h_minus, gamma_c_plus, gamma_c_minus, gamma_m_plus, gamma_m_minus;
  Real kappa_p_plus, kappa_p_minus, kappa_h_plus, kappa_h_minus, kappa_c_plus, kappa_c_minus, kappa_m_plus, kappa_m_minus;
  Real omega_p_plus, omega_p_minus, omega_h_plus, omega_h_minus, omega_c_plus, omega_c_minus, omega_m_plus, omega_m_minus;
  Real rho_p_plus, rho_p_minus, rho_h_plus, rho_h_minus, rho_c_plus, rho_c_minus;
  Real tol = 1.0e-1;
  dh = tol * _h[_qp]; dp = tol * _P[_qp]; dm = tol * _m[_qp]; dc = tol * _c[_qp];

  Iteration(_h[_qp], _P[_qp] + dp, _m[_qp], c_vect, aux1_p_plus, aux2_p_plus, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  gamma_p_plus = gamma(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  kappa_p_plus = kappa(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  omega_p_plus = omega(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  rho_p_plus = dummy7;
  _u_l_p_plus[_qp] = _u_l[_qp];
  _u_g_p_plus[_qp] = _u_g[_qp];
  _T_p_plus[_qp] = _T[_qp];
  //std::cout<<aux1_p_plus<<','<<aux2_p_plus<<','<<_u_l_p_plus[_qp]<<','<<_u_g_p_plus[_qp]<<std::endl;
  Iteration(_h[_qp], _P[_qp] - dp, _m[_qp], c_vect, aux1_p_minus, aux2_p_minus, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  gamma_p_minus = gamma(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  kappa_p_minus = kappa(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  omega_p_minus = omega(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  rho_p_minus = dummy7;
  _u_l_p_minus[_qp] = _u_l[_qp];
  _u_g_p_minus[_qp] = _u_g[_qp];
  _T_p_minus[_qp] = _T[_qp];
  //std::cout<<aux1_p_minus<<','<<aux2_p_minus<<','<<_u_l_p_minus[_qp]<<','<<_u_g_p_minus[_qp]<<std::endl;

  Iteration(_h[_qp] + dh, _P[_qp], _m[_qp], c_vect, aux1_h_plus, aux2_h_plus, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  gamma_h_plus = gamma(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  kappa_h_plus = kappa(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  omega_h_plus = omega(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  rho_h_plus = dummy7;
  _u_l_h_plus[_qp] = _u_l[_qp];
  _u_g_h_plus[_qp] = _u_g[_qp];
  _T_h_plus[_qp] = _T[_qp];
  //std::cout<<aux1_h_plus<<','<<aux2_h_plus<<','<<_u_l_h_plus[_qp]<<','<<_u_g_h_plus[_qp]<<std::endl;
  Iteration(_h[_qp] - dh, _P[_qp], _m[_qp], c_vect, aux1_h_minus, aux2_h_minus, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  gamma_h_minus = gamma(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  kappa_h_minus = kappa(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  omega_h_minus = omega(_m[_qp], dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  rho_h_minus = dummy7;
  _u_l_h_minus[_qp] = _u_l[_qp];
  _u_g_h_minus[_qp] = _u_g[_qp];
  _T_h_minus[_qp] = _T[_qp];
  //std::cout<<aux1_h_minus<<','<<aux2_h_minus<<','<<_u_l_h_minus[_qp]<<','<<_u_g_h_minus[_qp]<<std::endl;


  Iteration(_h[_qp], _P[_qp], _m[_qp] + dm, c_vect, aux1_m_plus, aux2_m_plus, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  gamma_m_plus = gamma(_m[_qp] + dm, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  kappa_m_plus = kappa(_m[_qp] + dm, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  omega_m_plus = omega(_m[_qp] + dm, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  _u_l_m_plus[_qp] = _u_l[_qp];
  _u_g_m_plus[_qp] = _u_g[_qp];
  //std::cout<<aux1_m_plus<<','<<aux2_m_plus<<','<<_u_l_m_plus[_qp]<<','<<_u_g_m_plus[_qp]<<std::endl;
  Iteration(_h[_qp], _P[_qp], _m[_qp] - dm, c_vect, aux1_m_minus, aux2_m_minus, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  gamma_m_minus = gamma(_m[_qp] - dm, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  kappa_m_minus = kappa(_m[_qp] - dm, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  omega_m_minus = omega(_m[_qp] - dm, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7);
  _u_l_m_minus[_qp] = _u_l[_qp];
  _u_g_m_minus[_qp] = _u_g[_qp];
  //std::cout<<aux1_m_minus<<','<<aux2_m_minus<<','<<_u_l_m_minus[_qp]<<','<<_u_g_m_minus[_qp]<<std::endl;

  _aux1_p[_qp] = (aux1_p_plus[3] - aux1_p_minus[3])/2.0/dp ;
  _aux2_p[_qp] = (aux2_p_plus[3] - aux2_p_minus[3])/2.0/dp ;
  _aux1_p_CH4[_qp] = (aux1_p_plus[2] - aux1_p_minus[2])/2.0/dp ;
  _aux2_p_CH4[_qp] = (aux2_p_plus[2] - aux2_p_minus[2])/2.0/dp ;
  _aux1_p_N2[_qp] = (aux1_p_plus[4] - aux1_p_minus[4])/2.0/dp ;
  _aux2_p_N2[_qp] = (aux2_p_plus[4] - aux2_p_minus[4])/2.0/dp ;
  _aux1_p_H2S[_qp] = (aux1_p_plus[5] - aux1_p_minus[5])/2.0/dp ;
  _aux2_p_H2S[_qp] = (aux2_p_plus[5] - aux2_p_minus[5])/2.0/dp ;
  _aux1_p_NaCl[_qp] = (aux1_p_plus[6] - aux1_p_minus[6])/2.0/dp ;
  _aux2_p_NaCl[_qp] = (aux2_p_plus[6] - aux2_p_minus[6])/2.0/dp ;
  _aux1_p_KCl[_qp] = (aux1_p_plus[7] - aux1_p_minus[7])/2.0/dp ;
  _aux2_p_KCl[_qp] = (aux2_p_plus[7] - aux2_p_minus[7])/2.0/dp ;
  _aux1_p_CaCl2[_qp] = (aux1_p_plus[8] - aux1_p_minus[8])/2.0/dp ;
  _aux2_p_CaCl2[_qp] = (aux2_p_plus[8] - aux2_p_minus[8])/2.0/dp ;
  _aux1_p_MgCl2[_qp] = (aux1_p_plus[9] - aux1_p_minus[9])/2.0/dp ;
  _aux2_p_MgCl2[_qp] = (aux2_p_plus[9] - aux2_p_minus[9])/2.0/dp ;
  _u_l_p[_qp] = (_u_l_p_plus[_qp] - _u_l_p_minus[_qp])/2.0/dp ;
  _u_g_p[_qp] = (_u_g_p_plus[_qp] - _u_g_p_minus[_qp])/2.0/dp ;
  _T_p[_qp] = (_T_p_plus[_qp] - _T_p_minus[_qp])/2.0/dp ;


  _aux1_h[_qp] = (aux1_h_plus[3] - aux1_h_minus[3])/2.0/dh ;
  _aux2_h[_qp] = (aux2_h_plus[3] - aux2_h_minus[3])/2.0/dh ;
  _aux1_h_CH4[_qp] = (aux1_h_plus[2] - aux1_h_minus[2])/2.0/dh ;
  _aux2_h_CH4[_qp] = (aux2_h_plus[2] - aux2_h_minus[2])/2.0/dh ;
  _aux1_h_N2[_qp] = (aux1_h_plus[4] - aux1_h_minus[4])/2.0/dh ;
  _aux2_h_N2[_qp] = (aux2_h_plus[4] - aux2_h_minus[4])/2.0/dh ;
  _aux1_h_H2S[_qp] = (aux1_h_plus[5] - aux1_h_minus[5])/2.0/dh ;
  _aux2_h_H2S[_qp] = (aux2_h_plus[5] - aux2_h_minus[5])/2.0/dh ;
  _aux1_h_NaCl[_qp] = (aux1_h_plus[6] - aux1_h_minus[6])/2.0/dh ;
  _aux2_h_NaCl[_qp] = (aux2_h_plus[6] - aux2_h_minus[6])/2.0/dh ;
  _aux1_h_KCl[_qp] = (aux1_h_plus[7] - aux1_h_minus[7])/2.0/dh ;
  _aux2_h_KCl[_qp] = (aux2_h_plus[7] - aux2_h_minus[7])/2.0/dh ;
  _aux1_h_CaCl2[_qp] = (aux1_h_plus[8] - aux1_h_minus[8])/2.0/dh ;
  _aux2_h_CaCl2[_qp] = (aux2_h_plus[8] - aux2_h_minus[8])/2.0/dh ;
  _aux1_h_MgCl2[_qp] = (aux1_h_plus[9] - aux1_h_minus[9])/2.0/dh ;
  _aux2_h_MgCl2[_qp] = (aux2_h_plus[9] - aux2_h_minus[9])/2.0/dh ;
  _u_l_h[_qp] = (_u_l_h_plus[_qp] - _u_l_h_minus[_qp])/2.0/dh ;
  _u_g_h[_qp] = (_u_g_h_plus[_qp] - _u_g_h_minus[_qp])/2.0/dh ;
  _T_h[_qp] = (_T_h_plus[_qp] - _T_h_minus[_qp])/2.0/dh ;

  //_aux1_c[_qp] = (aux1_c_plus - aux1_c_minus)/2.0/dc ;
  //_aux2_c[_qp] = (aux2_c_plus - aux2_c_minus)/2.0/dc ;
  //_u_l_c[_qp] = (_u_l_c_plus[_qp] - _u_l_c_minus[_qp])/2.0/dc ;
  //_u_g_c[_qp] = (_u_g_c_plus[_qp] - _u_g_c_minus[_qp])/2.0/dc ;
//  _T_c[_qp] = (_T_c_plus[_qp] - _T_c_minus[_qp])/2.0/dc ;

  _u_l_m[_qp] = (_u_l_m_plus[_qp] - _u_l_m_minus[_qp])/2.0/dm ;
  _u_g_m[_qp] = (_u_g_m_plus[_qp] - _u_g_m_minus[_qp])/2.0/dm ;

  _dgamma_dp[_qp]  = 0.0; _dgamma_dh[_qp]  = 0.0; _dgamma_dm[_qp]  = 0.0;
  _dkappa_dp[_qp] = 0.0; _dkappa_dh[_qp] = 0.0; _dkappa_dm[_qp] = 0.0;
  _dkappa_dph[_qp] = 0.0; _dkappa_dpm[_qp] = 0.0; _dkappa_dhm[_qp] = 0.0;
  _dkappa_dp2[_qp] = 0.0; _dkappa_dh2[_qp] = 0.0; _dkappa_dm2[_qp] = 0.0;
  _domega_dp[_qp] = 0.0; _domega_dh[_qp] = 0.0; _domega_dm[_qp] = 0.0;

  _dgamma_dp[_qp]  = (gamma_p_plus - gamma_p_minus)/2.0/dp;
  _dkappa_dp[_qp]  = (kappa_p_plus - kappa_p_minus)/2.0/dp;
  _domega_dp[_qp]  = (omega_p_plus - omega_p_minus)/2.0/dp;
  _drho_m_dp[_qp] =  (rho_p_plus - rho_p_minus)/2.0/dp;

  _dgamma_dh[_qp]  = (gamma_h_plus - gamma_h_minus)/2.0/dh;
  _dkappa_dh[_qp]  = (kappa_h_plus - kappa_h_minus)/2.0/dh;
  _domega_dh[_qp]  = (omega_h_plus - omega_h_minus)/2.0/dh;
  _drho_m_dh[_qp] = (rho_h_plus - rho_h_minus)/2.0/dh;

  _dgamma_dm[_qp]  = (gamma_m_plus - gamma_m_minus)/2.0/dm;
  _dkappa_dm[_qp]  = (kappa_m_plus - kappa_m_minus)/2.0/dm;
  _domega_dm[_qp]  = (omega_m_plus - omega_m_minus)/2.0/dm;

  //_drho_m_dc[_qp] = (rho_c_plus - rho_c_minus)/2.0/dc;
  //_dgamma_dc[_qp]  = (gamma_c_plus - gamma_c_minus)/2.0/dc;
  //_dkappa_dc[_qp]  = (kappa_c_plus - kappa_c_minus)/2.0/dc;
  //_domega_dc[_qp]  = (omega_c_plus - omega_c_minus)/2.0/dc;


}
