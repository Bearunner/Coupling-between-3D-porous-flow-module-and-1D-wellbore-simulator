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

#include "MoskitoDFHK.h"

registerMooseObject("MoskitoApp", MoskitoDFHK);

// template <>
// InputParameters
// validParams<MoskitoDFHK>()
// {
//   InputParameters params = validParams<MoskitoDriftFlux>();

InputParameters
MoskitoDFHK::validParams()
{
  InputParameters params = MoskitoDriftFlux::validParams();

  params.addRequiredParam<Real>("surface_tension",
        "Surface tension of the liquid phase (kg/s²)");
  return params;
}

MoskitoDFHK::MoskitoDFHK(const InputParameters & parameters)
  : MoskitoDriftFlux(parameters),
    _surf_ten(getParam<Real>("surface_tension"))
{
  _surf_ten *= kg_to_lbm;
}

void
MoskitoDFHK::DFMCalculator(MoskitoDFGVar & input) const
{
  // To avoid calculation for 1 phase flow and return correct value for DF method
  if (input._mfrac > 0.0 && input._mfrac < 1.0 && input._v_m > 0.0)
  {
    // conversion SI --> American petroleum units
    input._rho_g *= kg_to_lbm / m3_to_ft3;
    input._rho_l *= kg_to_lbm / m3_to_ft3;
    input._grav *= m_to_ft;
    input._dia *= m_to_ft;
    input._v_m *= m_to_ft;

    // To match new sign for flow direction
    input._dir *= -1.0;

    MoskitoHKLVar tmp;

    //check constraints of Hasan Kabir approach
    if (input._angle > 0.25 * PI)
      mooseError(name(), ": Angle > 45°, violating Hasan & Kabir limitations");

    HKinitialisation(input, tmp);
    HKcalculator(input, tmp);
    // HKvfrac(input, tmp);

    // conversion back to SI
    input._vd /= m_to_ft;


    // correction for ud sign for injection and production
    if (input._dir == -1.0)
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
MoskitoDFHK::HKinitialisation(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  cal_v_s(input, LVar);
  cal_vd_b(input, LVar);
  cal_vd_tb(input, LVar);
  cal_vd_mix(input, LVar);
  cal_v_gb(input, LVar);
  cal_v_ms(input, LVar);
  cal_v_gc(input, LVar);

  input._vd = LVar.vd_b;
  input._C0 = _C0b;
}

void
MoskitoDFHK::HKcalculator(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 0;

  // Check transition between slug and bubbly flow
  if (LVar.v_sg <= LVar.v_gb && LVar.vd_tb > LVar.vd_b)
  {
      // Check transition between dispersed bubbly and bubbly flow
      if (LVar.v_ms < input._v_m )
        Det_db_flow(input, LVar);
      else
        Det_bubbly_flow(input, LVar);
  }
  else
  {
    // Check transition from slug to churn flow and from d_bubbly to churn flow
    if (LVar.v_sg > 1.08 * LVar.v_sl && LVar.v_ms < input._v_m )
    {
      // Check transition between churn and annular flow
      if (LVar.v_sg < LVar.v_gc)
        Det_churn_flow(input, LVar);
      else
        Det_annular_flow(input, LVar);
    }
    else
    {
      // Check transition between slug and dispersed bubbly flow
      if (input._v_m < LVar.v_ms)
        Det_slug_flow(input, LVar);
      else
        Det_db_flow(input, LVar);
    }
  }
}

void
MoskitoDFHK::cal_v_s(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  if (input._dir == 1.0 || input._dir == -1.0) // cocurrent flow
    LVar.v_sg =  input._v_m * input._rho_l * input._mfrac  / ( - input._mfrac * input._rho_g + input._rho_l * input._mfrac + input._rho_g);
    input._v_sg = LVar.v_sg/m_to_ft;
  // else                 // countercurrent
  //   LVar.v_sg =  - input._v_m * input._rho_l * input._mfrac / ( - input._mfrac * input._rho_g  - input._rho_l * input._mfrac + input._rho_g);

  if (input._dir == 1.0 || input._dir == -1.0) // cocurrent flow
    LVar.v_sl = input._v_m - LVar.v_sg;
    input._v_sl = LVar.v_sl/m_to_ft;
  // else          // countercurrent
  //   LVar.v_sl = LVar.v_sg - input._v_m;
}

void
MoskitoDFHK::cal_vd_b(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.vd_b = 1.53 * std::pow(input._grav * _surf_ten  * (input._rho_l - input._rho_g) / std::pow(input._rho_l,2.0),1.0 / 4.0);
}

void
MoskitoDFHK::cal_vd_tb(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.vd_tb = 0.35 * std::pow(input._grav * input._dia  * (input._rho_l - input._rho_g) /
        input._rho_l ,0.5) * std::pow(std::cos(input._angle),0.5) *
        std::pow(1.0 + std::sin(input._angle),1.2);
}

void
MoskitoDFHK::cal_vd_mix(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.vd_mix = LVar.vd_b * (1.0 - std::exp(-0.1 * LVar.v_gb / (LVar.v_sg - LVar.v_gb))) +  LVar.vd_tb * std::exp(-0.1 * LVar.v_gb / (LVar.v_sg - LVar.v_gb));
}

void
MoskitoDFHK::cal_v_gb(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.v_gb = ((_C0b / (4.0 - _C0b)) * LVar.v_sl + (1.0 / (4.0 - _C0b)) * LVar.vd_b * input._dir) * cos(input._angle);
}

void
MoskitoDFHK::cal_v_gc(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.v_gc = 3.1 * std::pow(input._grav * _surf_ten * (input._rho_l - input._rho_g) / (std::pow(input._rho_g,2.0)),1.0 / 4.0);
}

void
MoskitoDFHK::cal_v_ms(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  LVar.v_ms = std::pow((0.725 + 4.15 * std::pow(LVar.v_sg / input._v_m, 0.5)) /
          (2.0 * std::pow(0.4 * _surf_ten  / ((input._rho_l - input._rho_g ) *
          input._grav),0.5) * std::pow(input._rho_l / _surf_ten ,0.6) *
          std::pow(input._friction / (2.0 * input._dia ), 0.4)), 1.0 / 1.2);
}

void
MoskitoDFHK::Det_bubbly_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 1;
  input._C0 = _C0b;
  input._vd = LVar.vd_b;
}

void
MoskitoDFHK::Det_db_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 2;
  input._C0 = _C0db;
  input._vd = LVar.vd_b;
}

void
MoskitoDFHK::Det_slug_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 3;
  if (input._dir == 1)
    input._C0 = _C0s_u;
  else
    input._C0 = interpol(_C0b, _C0s_d, LVar.v_gb, LVar.v_sg);

  input._vd = LVar.vd_mix;
}

void
MoskitoDFHK::Det_churn_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 4;
  if (input._dir == 1)
    input._C0 = interpol(_C0s_u, _C0c_u, LVar.v_ms, input._v_m);
  else
    input._C0 = _C0c_d;

  input._vd = LVar.vd_mix;
}

void
MoskitoDFHK::Det_annular_flow(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const
{
  input._FlowPat = 5;
  if (input._dir == 1)
    input._C0 = interpol(_C0c_u, _C0a, LVar.v_gc, LVar.v_sg);
  else
    input._C0 = interpol(_C0c_d, _C0a, LVar.v_gc, LVar.v_sg);
  input._vd = 0.0;
}

Real
MoskitoDFHK::interpol(const Real & C0_1, const Real & C0_2, const Real & v_denom, const Real & v_num) const
{
  return C0_1 * (1 - std::exp(-0.1 * v_denom / (v_num - v_denom))) +  C0_2 * std::exp(-0.1 * v_denom / (v_num - v_denom));
}

// Internal calculation of vfrac; not used
// void
// MoskitoDFHK::HKvfrac(MoskitoDFGVar & input, MoskitoHKLVar & LVar) const // Add for HK_evaluation ???
// {
//   input._vfrac = LVar.v_sg / (input._C0 * fabs(input._v_m) + input._vd * input._dir);
//
//   if (input._vfrac < 0.0 || isnan(input._vfrac))
//     input._vfrac = 0.0;
// }
