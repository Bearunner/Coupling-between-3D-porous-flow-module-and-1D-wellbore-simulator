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

#include "MoskitoComponent.h"

registerMooseObject("MoskitoApp", MoskitoComponent);

// template <>
// InputParameters
// validParams<MoskitoComponent>()
// {
//     InputParameters params = validParams<Material>();

InputParameters
MoskitoComponent::validParams()
{
  InputParameters params = Material::validParams();

    params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable (Pa)");
    params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable (J/kg)");
    params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable (m^3/s)");
    params.addRequiredCoupledVar("concentration", "concentration");
    //params.addCoupledVar("c2", 0.05, "Coupled value");
    params.addCoupledVar("c2", "Coupled value c2");
    params.addCoupledVar("c3", "Coupled value c3");
    params.addCoupledVar("c4", "Coupled value c4");
    params.addCoupledVar("c5", "Coupled value c5");
    params.addCoupledVar("c6", "Coupled value c6");
    params.addCoupledVar("c7", "Coupled value c7");
    params.addCoupledVar("c8", "Coupled value c8");

    params.addRequiredParam<UserObjectName>("eos_uo",
          "The name of the userobject for 2 phase EOS");
    params.addRequiredParam<UserObjectName>("drift_flux_uo",
          "The name of the userobject for drift flux model");
    params.addRequiredParam<UserObjectName>("viscosity_uo",
          "The name of the userobject for 2 phase viscosity Eq");
    params.addRequiredParam<MaterialPropertyName>("aux1_c_name", "The name of aux1_c");
    params.addRequiredParam<MaterialPropertyName>("aux2_c_name", "The name of aux2_c");
    params.addRequiredParam<MaterialPropertyName>("u_l_c_name", "The name of u_l_c");
    params.addRequiredParam<MaterialPropertyName>("u_g_c_name", "The name of u_g_c");
    params.addRequiredParam<MaterialPropertyName>("drho_m_dcc_name", "The name of aux1_c");
    params.addRequiredParam<MaterialPropertyName>("dgamma_dcc_name", "The name of aux2_c");
    params.addRequiredParam<MaterialPropertyName>("dkappa_dcc_name", "The name of aux2_c");
    params.addRequiredParam<MaterialPropertyName>("domega_dcc_name", "The name of u_l_c");
    params.addRequiredParam<MaterialPropertyName>("T_cc_name", "The name of u_g_c");
    params.addRequiredParam<int>("component", "The component");
    //params.addParam<int>("index_CO2", 100, "CO2 index in variable vector");
    //params.addRequiredParam<std::string>("u_g_c_name", "The name of u_g_c");


    return params;
}

MoskitoComponent::MoskitoComponent(const InputParameters & parameters)
  : Material(parameters),
  eos_uo(getUserObject<MoskitoEOS2P>("eos_uo")),
  dfm_uo(getUserObject<MoskitoDriftFlux>("drift_flux_uo")),
  viscosity_uo(getUserObject<MoskitoViscosity2P>("viscosity_uo")),
  _area(getMaterialProperty<Real>("well_area")),
  _friction(getMaterialProperty<Real>("well_moody_friction")),
  _dia(getMaterialProperty<Real>("well_diameter")),
  //_u_f(getParam<Real>("manual_friction_factor")),
  _well_sign(getMaterialProperty<Real>("flow_direction_sign")),
  _gravity(getMaterialProperty<RealVectorValue>("gravity")),
  _well_dir(getMaterialProperty<RealVectorValue>("well_direction_vector")),
  //_aux1_p(declareProperty<Real>("aux1_p")),
  //_aux2_p(declareProperty<Real>("aux2_p")),
  //_u_l_p(declareProperty<Real>("u_l_p")),
  //_u_g_p(declareProperty<Real>("u_g_p")),
  //_aux1_h(declareProperty<Real>("aux1_h")),
  //_aux2_h(declareProperty<Real>("aux2_h")),
  //_u_l_h(declareProperty<Real>("u_l_h")),
  //_u_g_h(declareProperty<Real>("u_g_h")),
  _aux1_c(declareProperty<Real>("aux1_c_name")),
  _aux2_c(declareProperty<Real>("aux2_c_name")),
  _u_l_c(declareProperty<Real>("u_l_c_name")),
  _u_g_c(declareProperty<Real>("u_g_c_name")),
  _drho_m_dcc(declareProperty<Real>("drho_m_dcc_name")),
  _dgamma_dcc(declareProperty<Real>("dgamma_dcc_name")),
  _dkappa_dcc(declareProperty<Real>("dkappa_dcc_name")),
  _domega_dcc(declareProperty<Real>("domega_dcc_name")),
  _T_cc(declareProperty<Real>("T_cc_name")),
  _component(getParam<int>("component")),
  //_ico2(getParam<int>("index_CO2")),
  _P(coupledValue("pressure")),
  _h(coupledValue("enthalpy")),
  _m(coupledValue("massrate")),
  _c(coupledValue("concentration")),
  //_c2(coupledValue("c2"))
  _c2(isCoupled("c2") ? coupledValue("c2") : _zero),
  _c3(isCoupled("c3") ? coupledValue("c3") : _zero),
  _c4(isCoupled("c4") ? coupledValue("c4") : _zero),
  _c5(isCoupled("c5") ? coupledValue("c5") : _zero),
  _c6(isCoupled("c6") ? coupledValue("c6") : _zero),
  _c7(isCoupled("c7") ? coupledValue("c7") : _zero),
  _c8(isCoupled("c8") ? coupledValue("c8") : _zero)
{
}

void
MoskitoComponent::computeQpProperties()
{

  c_vect[0] = _c[_qp];
  c_vect[1] = _c2[_qp];
  c_vect[2] = _c3[_qp];
  c_vect[3] = _c4[_qp];
  c_vect[4] = _c5[_qp];
  c_vect[5] = _c6[_qp];
  c_vect[6] = _c7[_qp];
  c_vect[7] = _c8[_qp];

  //cout<<"component"<<c_vect[1]<<endl;

  //cout<<"component c_vect[0] = "<<c_vect[0]<<" c_vect[1] = "<<c_vect[1]<<" c_vect[2] = "<<c_vect[2]<<endl;

  Iteration(_h[_qp], _P[_qp], _m[_qp], c_vect, _aux1_new, _aux2_new, dummy1, dummy2,
     dummy3, dummy4, dummy5, dummy6, dummy7, u_g, u_l, gamma, kappa, omega, temper);

  Derivatives_c();

     //std::cout<<u_g<<','<<u_l<<" sec"<<std::endl;

}

void
MoskitoComponent::Iteration(const Real & h, const Real & p, const Real & m,
   const Real c_vect[], Real aux1_arr[], Real aux2_arr[], Real & g_ro, Real & l_rho, Real & lhx, Real & gas_h, Real & mass_frac, Real & non_aqueous_mf, Real & rho_mix,
 Real & u_g, Real & u_l, Real & gamma, Real & kappa, Real & omega, Real & temper)
{

    Real x[6], u, _Rey, _rho_pa, _flow_path, _c00, _u_dd, vis, con;

    eos_uo.Itera(h, p, m, c_vect, g_ro, l_rho, lhx, gas_h, mass_frac, non_aqueous_mf, rho_mix, temper, x, aux1_arr, aux2_arr, vis, con);

    u = fabs(m) / rho_mix / _area[_qp];

    //_Rey = rho_mix * _dia[_qp] * u / viscosity_uo.mixture_mu(p, temper, non_aqueous_mf);
    _Rey = rho_mix * _dia[_qp] * u / vis;

    // drift-flux calculator section
      MoskitoDFGVar DFinp(u, g_ro, l_rho, mass_frac, non_aqueous_mf,
        _dia[_qp], _well_sign[_qp], _friction[_qp], _gravity[_qp], _well_dir[_qp]);
      dfm_uo.DFMCalculator(DFinp);
      DFinp.DFMOutput(_flow_path, _c00, _u_dd);

    _rho_pa = g_ro * _c00  * non_aqueous_mf + (1.0 - non_aqueous_mf * _c00) * l_rho;
////////////////////////////////////////////////////////////////////////////////
u_g, u_l;
    if (u != 0.0)
    {
      u_g  = _c00 * rho_mix * u + l_rho * _u_dd;
      u_g /= _rho_pa;
      u_l  = (1.0 - non_aqueous_mf * _c00) * rho_mix  * u - g_ro * non_aqueous_mf * _u_dd;
      u_l /= (1.0 - non_aqueous_mf) * _rho_pa;
    }
////////////////////////////////////////////////////////////////////////////////
gamma = 0.0;

if(m != 0.0)
{
  gamma  = non_aqueous_mf / (1.0 - non_aqueous_mf);
  gamma *= g_ro * l_rho * rho_mix / (_rho_pa * _rho_pa);
  gamma *= std::pow((_c00 - 1.0) * (m / rho_mix / _area[_qp]) + _u_dd , 2.0);
}
  ////////////////////////////////////////////////////////////////////////////////
kappa = 0.0;

  Real h_g, h_l, rho_m;

  if(m != 0.0)
  {
    rho_m = l_rho * g_ro / (mass_frac * (l_rho - g_ro) + g_ro);

    h_g = gas_h;
    h_l = lhx;

    kappa  = non_aqueous_mf * g_ro * l_rho / _rho_pa * (h_g - h_l);
    kappa *= (_c00 - 1.0) * (m / rho_m / _area[_qp]) + _u_dd;
  }
  ////////////////////////////////////////////////////////////////////////////////
   omega = 0.0;

  if(m != 0.0)
  {
    rho_m = rho_mix;

    Real v = m / rho_m / _area[_qp] ;

    u_g  = (_c00 * rho_m * v + l_rho * _u_dd) / _rho_pa;
    u_l  = (1.0 - non_aqueous_mf * _c00) * rho_m * v - g_ro * non_aqueous_mf * _u_dd;
    u_l /= (1.0 - non_aqueous_mf) * _rho_pa;

    omega -= 3.0 * u_g * u_l * v;
    omega += std::pow(u_g,3.0) * (1.0 + g_ro * non_aqueous_mf /rho_m);
    omega += std::pow(u_l,3.0) * (1.0 + l_rho * (1.0 - non_aqueous_mf) / rho_m);
    omega *= 0.5 * non_aqueous_mf * (1.0 - non_aqueous_mf) * g_ro * l_rho / rho_m;
  }
//Itera(h, p, m, c_vect, g_ro, l_rho, lhx, gas_h, mass_frac, non_aqueous_mf, rho_mix, temper, x, aux1_arr, aux2_arr);
//cout<<g_ro<<" "<<l_rho<<" "<<lhx<<" "<<gas_h<<" "<<mass_frac<<" "<<non_aqueous_mf<<" "<<rho_mix<<" "<<temper<<" "<<aux1_arr[4]<<" "<<aux2_arr[4]<<" "<<endl;
//cout<<"component  "<<g_ro<<" "<<l_rho<<" "<<lhx<<" "<<gas_h<<" "<<mass_frac<<" "<<non_aqueous_mf<<" "<<rho_mix<<" "<<temper<<" "<<aux1_arr[4]<<" "<<aux2_arr[4]<<" "<<endl;
//cout<<aux1_arr[1]<< " "<<aux1_arr[2]<< " "<<aux1_arr[3]<< " "<<aux1_arr[4]<< " "<<aux1_arr[5]<< " "<<endl;
//cout<<g_ro<<" "<<_component<<endl;
//cout<<l_rho<<endl;
}

void
MoskitoComponent::Derivatives_c()
{

  Real dh, dm, dp, aux1_p_plus, aux2_p_plus, aux1_p_minus, aux2_p_minus;
  Real aux1_h_plus, aux2_h_plus, aux1_h_minus, aux2_h_minus;
  Real aux1_c_plus[10], aux2_c_plus[10], aux1_c_minus[10], aux2_c_minus[10];
  Real dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
  Real gamma_pp_plus, gamma_pp_minus, gamma_hh_plus, gamma_hh_minus, gamma_cc_plus, gamma_cc_minus;
  Real kappa_pp_plus, kappa_pp_minus, kappa_hh_plus, kappa_hh_minus, kappa_cc_plus, kappa_cc_minus;
  Real omega_pp_plus, omega_pp_minus, omega_hh_plus, omega_hh_minus, omega_cc_plus, omega_cc_minus;
  Real _u_g_p_plus, _u_g_p_minus, _u_g_h_plus, _u_g_h_minus, _u_g_c_plus, _u_g_c_minus;
  Real _u_l_p_plus, _u_l_p_minus, _u_l_h_plus, _u_l_h_minus, _u_l_c_plus, _u_l_c_minus;
  Real rho_p_plus, rho_p_minus, rho_h_plus, rho_h_minus, rho_c_plus, rho_c_minus, temper_cc_plus, temper_cc_minus;
  Real tol = 1.0e-1;
  dh = tol * _h[_qp]; dp = tol * _P[_qp]; dm = tol * _m[_qp];

/*
  Real dc[2];
  Real c_vect_plus[2];
  int i = _component;
  dc[i] = c_vect[i]*0.1;
  c_vect_plus[0] = c_vect[0];
  c_vect_plus[1] = c_vect[1];
  c_vect_plus[i] = c_vect[i] + dc[i];
*/
  Real dc[8];
  Real c_vect_plus[8];
  int i = _component;
  dc[i] = c_vect[i]*0.1;
  if (dc[i] == 0.0)
  {
    dc[i] = 0.0001;
  }
  c_vect_plus[0] = c_vect[0];
  c_vect_plus[1] = c_vect[1];
  c_vect_plus[2] = c_vect[2];
  c_vect_plus[3] = c_vect[3];
  c_vect_plus[4] = c_vect[4];
  c_vect_plus[5] = c_vect[5];
  c_vect_plus[6] = c_vect[6];
  c_vect_plus[7] = c_vect[7];
  c_vect_plus[i] = c_vect[i] + dc[i];

  Iteration(_h[_qp], _P[_qp], _m[_qp], c_vect_plus, aux1_c_plus, aux2_c_plus,
     dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, rho_c_plus, _u_g_c_plus, _u_l_c_plus, gamma_cc_plus, kappa_cc_plus, omega_cc_plus, temper_cc_plus);

  //std::cout<<aux1_c_plus<<','<<aux2_c_plus<<','<<_u_l_c_plus[_qp]<<','<<_u_g_c_plus[_qp]<<std::endl;
/*
  Real c_vect_minus[2];
  dc[i] = c_vect[i]*0.1;
  if (dc[i] == 0.0)
  {
    dc[i] = 0.0001;
  }
  c_vect_minus[0] = c_vect[0];
  c_vect_minus[1] = c_vect[1];
  c_vect_minus[i] = c_vect[i] - dc[i];


  Real c_vect_minus[8];
  dc[i] = c_vect[i]*0.1;
  if (dc[i] == 0.0)
  {
    dc[i] = 0.0001;
  }
  c_vect_minus[0] = c_vect[0];
  c_vect_minus[1] = c_vect[1];
  c_vect_minus[2] = c_vect[2];
  c_vect_minus[3] = c_vect[3];
  c_vect_minus[4] = c_vect[4];
  c_vect_minus[5] = c_vect[5];
  c_vect_minus[6] = c_vect[6];
  c_vect_minus[7] = c_vect[7];
  c_vect_minus[i] = c_vect[i] + dc[i];
*/

  c_vect_plus[i] = c_vect[i] - dc[i];

  Iteration(_h[_qp], _P[_qp], _m[_qp], c_vect_plus, aux1_c_minus, aux2_c_minus,
     dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, rho_c_minus, _u_g_c_minus, _u_l_c_minus, gamma_cc_minus, kappa_cc_minus, omega_cc_minus, temper_cc_minus);

  //std::cout<<aux1_c_plus[i]<<','<<aux1_c_minus[i]<<','<<_u_l_c_minus<<','<<_u_g_c_minus<<std::endl;

  //_aux1_p[_qp] = (aux1_p_plus - aux1_p_minus)/2.0/dp ;
  //_aux2_p[_qp] = (aux2_p_plus - aux2_p_minus)/2.0/dp ;
  //_T_p[_qp] = (_T_p_plus[_qp] - _T_p_minus[_qp])/2.0/dp ;


  //_aux1_h[_qp] = (aux1_h_plus - aux1_h_minus)/2.0/dh ;
  //_aux2_h[_qp] = (aux2_h_plus - aux2_h_minus)/2.0/dh ;
  //_T_h[_qp] = (_T_h_plus[_qp] - _T_h_minus[_qp])/2.0/dh ;

  _aux1_c[_qp] = (aux1_c_plus[i+2] - aux1_c_minus[i+2])/2.0/dc[i] ;
  _aux2_c[_qp] = (aux2_c_plus[i+2] - aux2_c_minus[i+2])/2.0/dc[i] ;
  _u_l_c[_qp] = (_u_l_c_plus - _u_l_c_minus)/2.0/dc[i] ;
  _u_g_c[_qp] = (_u_g_c_plus - _u_g_c_minus)/2.0/dc[i] ;
  _T_cc[_qp] = (temper_cc_plus - temper_cc_minus)/2.0/dc[i] ;


  _drho_m_dcc[_qp] = (rho_c_plus - rho_c_minus)/2.0/dc[i];
  _dgamma_dcc[_qp]  = (gamma_cc_plus - gamma_cc_minus)/2.0/dc[i];
  _dkappa_dcc[_qp]  = (kappa_cc_plus - kappa_cc_minus)/2.0/dc[i];
  _domega_dcc[_qp]  = (omega_cc_plus - omega_cc_minus)/2.0/dc[i];



}
