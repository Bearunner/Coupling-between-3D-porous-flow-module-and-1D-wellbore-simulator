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

#include "MoskitoFluidWell_1p1c.h"

registerMooseObject("MoskitoApp", MoskitoFluidWell_1p1c);

// template <>
// InputParameters
// validParams<MoskitoFluidWell_1p1c>()
// {
//   InputParameters params = validParams<MoskitoFluidWellGeneral>();

InputParameters
MoskitoFluidWell_1p1c::validParams()
{
  InputParameters params = MoskitoFluidWellGeneral::validParams();

  params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable (K)");
  params.addRequiredCoupledVar("flowrate", "Mixture flow rate nonlinear variable (m^3/s)");
  params.addRequiredParam<UserObjectName>("eos_uo",
        "The name of the userobject for EOS");
  params.addRequiredParam<UserObjectName>("viscosity_uo",
        "The name of the userobject for viscosity Eq");

  return params;
}

MoskitoFluidWell_1p1c::MoskitoFluidWell_1p1c(const InputParameters & parameters)
  : MoskitoFluidWellGeneral(parameters),
    eos_uo(getUserObject<MoskitoEOS1P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity1P>("viscosity_uo")),
    _hf(declareProperty<Real>("convective_heat_factor")),
    _vis(declareProperty<Real>("viscosity")),
    _dmu_dp(declareProperty<Real>("dmu_dp")),
    _dmu_dT(declareProperty<Real>("dmu_dT")),
    _lambda(declareProperty<Real>("fluid_thermal_conductivity")),
    _cp(declareProperty<Real>("specific_heat")),
    _rho(declareProperty<Real>("density")),
    _drho_dp(declareProperty<Real>("drho_dp")),
    _drho_dT(declareProperty<Real>("drho_dT")),
    _h(declareProperty<Real>("h_from_p_T")),
    _T(coupledValue("temperature")),
    _flow(coupledValue("flowrate")),
    _wv(declareProperty<RealVectorValue>("well_velocity_vector")),
    _dwv_dT(declareProperty<RealVectorValue>("d_well_velocity_dT")),
    _dwv_dp_phi(declareProperty<RealVectorValue>("d_well_velocity_dp_phi")),
    _dwv_dp_gradphi(declareProperty<RankTwoTensor>("d_well_velocity_dp_gradphi"))
{
}

void
MoskitoFluidWell_1p1c::computeQpProperties()
{
  MoskitoFluidWellGeneral::computeQpProperties();

  _vis[_qp] = viscosity_uo.mu(_P[_qp], _T[_qp]);
  // _dmu_dp[_qp] = 0;
  // _dmu_dT[_qp] = viscosity_uo.dmu_dpT(_P[_qp], _T[_qp], _vis[_qp], _dmu_dp[_qp], _dmu_dT[_qp]);

  // _vis[_qp] = viscosity_uo.mu(_P[_qp], _T[_qp]);
  _lambda[_qp] = eos_uo.lambda(_P[_qp], _T[_qp]);
  _cp[_qp] = eos_uo.cp(_P[_qp], _T[_qp]);
  _h[_qp] = eos_uo.h_from_p_T(_P[_qp], _T[_qp]);
  eos_uo.rho_from_p_T(_P[_qp], _T[_qp], _rho[_qp], _drho_dp[_qp], _drho_dT[_qp]);

  _u[_qp] = _flow[_qp] / _area[_qp];
  _well_P[_qp] = _P[_qp];
  _well_T[_qp] = _T[_qp];
  _well_q[_qp] = _q[_qp];
  // _constant
  // RealVectorValue _wv(-_flow[_qp] / _area[_qp],  0.0,  0.0);// not working
  // _wv[_qp](-_flow[_qp] / _area[_qp],  0.0,  0.0);// not working
  RealVectorValue wv(-_flow[_qp] / _area[_qp],  0.0,  0.0);
  _wv[_qp] = wv;
  _dwv_dT[_qp] = 0.0;
  _dwv_dp_phi[_qp] = 0.0;
  _dwv_dp_gradphi[_qp] = 0.0;
  _Re[_qp] = _rho[_qp] * _dia[_qp] * fabs(_u[_qp]) / _vis[_qp];
  if (_f_defined)
    _friction[_qp] = _u_f;
  else
    MoodyFrictionFactor(_friction[_qp], _rel_roughness, _Re[_qp], _roughness_type);

  _hf[_qp] = Conv_coeff();
  // _lambda[_qp]  = (1.0 - (_d * _d) / std::pow(_d + _thickness , 2.0)) * _lambda0;
  // _lambda[_qp] += (_d * _d) / std::pow(_d + _thickness , 2.0) * eos_uo.lambda(_P[_qp], _T[_qp]);
}

Real
MoskitoFluidWell_1p1c::Conv_coeff()
{
  Real pr_b, nusselt = 0.0, gamma;
  if (_Re[_qp]>0.0)
  {
    pr_b = _vis[_qp] * _cp[_qp] / _lambda[_qp];
    if (_Re[_qp]<2300.0)
      nusselt = 4.364 ;
    else if (_Re[_qp]<10000.0)
    {
      gamma = (_Re[_qp] - 2300.0)/(10000.0 - 2300.0);
      nusselt = (1.0 - gamma) * 4.364 + gamma * 0.023 * pow(_Re[_qp], 0.8) * pow(pr_b, 0.3);
    }
    else
      nusselt = 0.023 * pow(_Re[_qp], 0.8) * pow(pr_b, 0.3);
  }

  return nusselt * _lambda[_qp] / _H_dia[_qp];
}
