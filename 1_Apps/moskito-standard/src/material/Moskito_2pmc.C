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

#include "Moskito_2pmc.h"

registerMooseObject("MoskitoApp", Moskito_2pmc);
//
// template <>
// InputParameters
// validParams<Moskito_2pmc>()
// {
//     InputParameters params = validParams<Material>();

InputParameters
Moskito_2pmc::validParams()
{
  InputParameters params = Material::validParams();

    //params += validParams<NewtonIteration>();
    params.addClassDescription("describe it!");
    params.addRequiredCoupledVar("enthalpy", "Specific enthalpy nonlinear variable (J/kg)");
    params.addRequiredParam<UserObjectName>("eos_2p",
          "The name of the userobject for 2 phase EOS");
    params.addRequiredCoupledVar("massrate", "Massrate nonlinear variable (m^3/s)");
    params.addRequiredCoupledVar("pressure", "pressure nonlinear variable (m^3/s)");
    return params;
}

Moskito_2pmc::Moskito_2pmc(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    eos_2p(getUserObject<MoskitoEOS2P>("eos_2p")),
    _rho_gas(declareProperty<Real>("gas_den")),
    _rho_liq(declareProperty<Real>("liquid_den")),
    _rho_mix(declareProperty<Real>("densitymix")),
    _drho_mix_dp(declareProperty<Real>("drho_mix_dp")),
    _drho_mix_dp2(declareProperty<Real>("drho_mix_dp2")),
    _h(coupledValue("enthalpy")),
    _m(coupledValue("massrate")),
    _p(coupledValue("pressure")),
    //_p_var(coupled("pressure")),
    _p_name(getVar("pressure", 0)->name()),
    _grad_m(coupledGradient("massrate")),
    _grad_h(coupledGradient("enthalpy")),
    _grad_p(coupledGradient("pressure")),
    _d_rho_mixd_p(declarePropertyDerivative<Real>("densitymix", _p_name)),
    _new_gradient(declareProperty<Real>("new_gradient"))
{
}

void
Moskito_2pmc::computeQpProperties()
{
  //eos_2p.rho_m_by_p(_p[_qp], _h[_qp], _rho_mix[_qp], _drho_mix_dp[_qp], _drho_mix_dp2[_qp]);
  _new_gradient[_qp] = 1.0;// d_rho_mixd_p[_qp];

  //std::cout<<_new_gradient[_qp];
  //break;
  //real i;




}
