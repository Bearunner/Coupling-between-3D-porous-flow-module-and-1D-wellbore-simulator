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

#include "MoskitoDFShi.h"

registerMooseObject("MoskitoApp", MoskitoDFShi);

// template <>
// InputParameters
// validParams<MoskitoDFShi>()
// {
//     InputParameters params = validParams<MoskitoDriftFlux>();

InputParameters
MoskitoDFShi::validParams()
{
  InputParameters params = MoskitoDriftFlux::validParams();

    params.addRequiredParam<Real>("surface_tension",
          "Surface tension of the liquid phase (kg/s²)");
    params.addParam<Real>("Shi_param_cKu",142,
          "Fitting parameter for Kutateladze parameter; dimesionless; 75 < cKu < 142; original value was 75 (Richter), Pan uses 142");
    params.addParam<Real>("Shi_param_a1",0.06,
          "Fitting parameter for transition curve; dimesionless; a1 = 0.06");
    params.addParam<Real>("Shi_param_a2",0.21,
          "Fitting parameter for transition curve; dimesionless; 0.12 < a2 < 0.21");
    params.addParam<Real>("Shi_param_Fv",1.0,
          "Fitting parameter to increase C0's Sensitivty towards gas-water fractions ; dimesionless ");
    params.addParam<Real>("Pan_param_cMax",1.2,
          "Fitting parameter for transition curve called A in Shi et al. 2005; dimesionless; 1.0 < cMax < 1.5; 1.2 is the cefault value in ECLIPSE; the larger values (=>1.2 match smaller diameter borehole, for larger diameters cMax should approach 1.0 ");
    // params.addRequiredParam<std::vector<Real>>("inclination_vector",
    //       "Parameterization of Shi inclination function; 1.27 < m0 < 1.85;"
    //        "0.21 < n1 0.24; 0.95 < n2 1.08");
    return params;
}

MoskitoDFShi::MoskitoDFShi(const InputParameters & parameters)
  : MoskitoDriftFlux(parameters),
    _surf_ten(getParam<Real>("surface_tension")),
    _cKu(getParam<Real>("Shi_param_cKu")),
    _tran_param_a1(getParam<Real>("Shi_param_a1")),
    _tran_param_a2(getParam<Real>("Shi_param_a2")),
    _C0_Fv(getParam<Real>("Shi_param_Fv")),
    _C0_cMax(getParam<Real>("Pan_param_cMax"))
    // _Shi_incl_triple(getParam<std::vector<Real>>("inclination_vector"))
    {
    }

void
MoskitoDFShi::DFMCalculator(MoskitoDFGVar & input) const
{
  MoskitoShiLVar tmp;
  if (input._mfrac > 0.0 && input._mfrac < 1.0 && input._v_m > 0.0)
  {
      //check constraints of Shi approach limiting the approach to a deviation of 70° from vertical
      if (input._angle > 0.3888 * PI)
        mooseError(name(), ": Angle > 70°, violating Shi limitations");

    Shiinitialisation(input, tmp);
    Shicalculator(input, tmp);

    // correction for ud sign for injection and production
    if (input._dir == 1.0)
      input._vd *= -1.0;
  }
  else
  {
    input._vd = 0.0;
    input._C0 = 1.0;
    input._FlowPat = 0.0;
  }
}

void
MoskitoDFShi::Shiinitialisation(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const
{
  cal_Kutateladze(input, LVar);
  cal_velocities(input, LVar);
}

void
MoskitoDFShi::cal_Kutateladze(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const
{
	//Bond number is needed to calculate Kutateladze number
  LVar.Bond = input._dia * input._dia;
  LVar.Bond *= (input._grav * (input._rho_l - input._rho_g)/_surf_ten);

   //Kutateladze number is needed as a parameter to compute v_sgf and k_fct, Used formula by Pan et al 2011
  LVar.Kutateladze = std::pow(_cKu / std::pow(LVar.Bond,0.5) * (std::pow(1.0 + LVar.Bond / (_cKu * _cKu * _cw) ,0.5) - 1.0),0.5);
}

void
MoskitoDFShi::cal_velocities(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const
{
	//v_c is the characteristic velocity, used to calculate v_sgf and in the final Vd-function
  LVar.v_c = std::pow(input._grav * _surf_ten * (input._rho_l - input._rho_g) / (input._rho_l * input._rho_l),0.25);
    //v_sgf is a C0-parameter used in C0_B
  LVar.v_sgf = LVar.Kutateladze * std::pow(input._rho_l / input._rho_g,0.5) * LVar.v_c;
}

void
MoskitoDFShi::Shicalculator(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const
{
	//C0 is the parameter quantifying the cross section (liquid-gas mixture) of the well; used in _vd and in well2pMaterial
  input._C0 = _C0_cMax / (1.0 + (_C0_cMax - 1.0) * cal_auxvar_nu(input, LVar) * cal_auxvar_nu(input, LVar));

  input._vd = (1.0 - input._C0 * input._vfrac) * LVar.v_c * cal_transition_fct(input,LVar);
  input._vd /= input._C0 * input._vfrac * std::pow(input._rho_g / input._rho_l ,0.5) + 1.0 - input._vfrac * input._C0;
  // inclinaciton correction after HK and Shi
  input._vd *= std::pow(std::cos(input._angle),0.5) * std::pow(1.0 + std::sin(input._angle),1.2);
  // Pan used correction for inclination different from HK and Shi (see next line). If activatd inclination vector has to be activated also
  // input._vd *= _Shi_incl_triple[0] * std::pow(std::cos(input._angle), _Shi_incl_triple[1]) * std::pow(1.0 + std::sin(input._angle),_Shi_incl_triple[2]);
}

Real
MoskitoDFShi::cal_auxvar_nu(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const
{
	//C0-parameter; for convergence should be B<(2-C0_cMax)/C0_cMax (Shi et al), here appraoch used after Pan et al 2011
  Real C0_B = 2.0 / _C0_cMax - 1.0667;

	//C0-parameter; 0<C0_beta<1
  // Shi et al used the following expression
  // Real C0_beta = std::max(input._vfrac, _C0_Fv * input._vfrac * input._v_m / LVar.v_sgf);
  // which is not used but simplified as Shis formual leads to discontinuous curve for u_d
  Real C0_beta = input._vfrac ;

  Real nu = (C0_beta - C0_B) / (1.0 - C0_B);
  if (nu < 0)
    nu = 0;
  //if beta is < 1, which is the case in the simplified way, we don#t need to constrain nu
  // if (nu > 1.0)
  //   nu = 1.0;

  return nu;
}

Real
MoskitoDFShi::cal_transition_fct(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const
{
	//transitionfunction K; interpolates continously between bubbly and flooding film conditions
  Real K_fct;
  if (input._vfrac <= _tran_param_a1)
    K_fct = 1.53;
  else
  {
    if (input._vfrac <= _tran_param_a2)
    {
      K_fct = input._C0 * LVar.Kutateladze - 1.53;
      K_fct /= 2.0;
      K_fct *= (1.0 - std::cos(PI * (input._vfrac - _tran_param_a1) / (_tran_param_a2 -_tran_param_a1)));
      K_fct += 1.53;
    }
    else
      K_fct = input._C0 * LVar.Kutateladze;

  }
  return K_fct;
}
