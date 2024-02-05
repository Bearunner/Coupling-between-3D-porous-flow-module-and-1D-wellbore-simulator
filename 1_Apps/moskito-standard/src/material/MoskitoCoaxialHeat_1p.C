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

#include "MoskitoCoaxialHeat_1p.h"

registerMooseObject("MoskitoApp", MoskitoCoaxialHeat_1p);

// template <>
// InputParameters
// validParams<MoskitoCoaxialHeat_1p>()
// {
//     InputParameters params = validParams<Material>();

InputParameters
MoskitoCoaxialHeat_1p::validParams()
{
  InputParameters params = Material::validParams();

    params.addClassDescription("Materials for the Lateral heat transfer between "
          "inner and outer pipes");
    params.addRequiredCoupledVar("temperature_inner_pipe", "Temperature of fluid inside inner pipe (K)");
    params.addRequiredCoupledVar("flowrate_inner_pipe", "Mixture flow rate inside inner pipe (m^3/s)");
    params.addRequiredCoupledVar("pressure_inner_pipe", "Pressure inside inner pipe (Pa)");
    params.addRequiredCoupledVar("temperature_outer_pipe", "Temperature of fluid inside outer pipe (K)");
    params.addRequiredCoupledVar("flowrate_outer_pipe", "Mixture flow rate inside outer pipe (m^3/s)");
    params.addRequiredCoupledVar("pressure_outer_pipe", "Pressure inside outer pipe (Pa)");
    params.addRequiredParam<Real>("inner_pipe_outer_radius",
          "outer radius of the inner pipe (m)");
    params.addRequiredParam<Real>("inner_pipe_wall_thickness",
          "wall thickness of the inner pipe (m)");
    params.addRequiredParam<Real>("outer_pipe_inner_radius",
          "inner radius of the outer pipe (m)");
    params.addRequiredParam<Real>("conductivity_inner_pipe",
          "Thermal conductivity of the inner pipe (W/(m*K))");
    params.addRequiredParam<UserObjectName>("eos_uo",
          "The name of the userobject for EOS");
    params.addRequiredParam<UserObjectName>("viscosity_uo",
          "The name of the userobject for viscosity Eq");

    return params;
}

MoskitoCoaxialHeat_1p::MoskitoCoaxialHeat_1p(const InputParameters & parameters)
  : Material(parameters),
    eos_uo(getUserObject<MoskitoEOS1P>("eos_uo")),
    viscosity_uo(getUserObject<MoskitoViscosity1P>("viscosity_uo")),
    _ohc(declareProperty<Real>("overall_heat_transfer_coeff")),
    _rio(getParam<Real>("inner_pipe_outer_radius")),
    _wt(getParam<Real>("inner_pipe_wall_thickness")),
    _roi(getParam<Real>("outer_pipe_inner_radius")),
    _ki(getParam<Real>("conductivity_inner_pipe")),
    _T_i(coupledValue("temperature_inner_pipe")),
    _flow_i(coupledValue("flowrate_inner_pipe")),
    _p_i(coupledValue("pressure_inner_pipe")),
    _T_o(coupledValue("temperature_outer_pipe")),
    _flow_o(coupledValue("flowrate_outer_pipe")),
    _p_o(coupledValue("pressure_outer_pipe")),
    _nusselt_i(declareProperty<Real>("Nusselt_number_inner_pipe")),
    _nusselt_o(declareProperty<Real>("Nusselt_number_outer_pipe"))
{
}

void
MoskitoCoaxialHeat_1p::computeQpProperties()
{
  Real j = 0.0 ;
  j += 1.0 / _rio / Conv_coeff_outer();
  j += log(_rio /(_rio - _wt))/_ki;
  j += 1.0 / (_rio - _wt) / Conv_coeff_inner();
  _ohc[_qp] = 1.0 / j ;
}

Real
MoskitoCoaxialHeat_1p::Conv_coeff_inner()
{
  Real pr_i, gama_i, area_i, u_i, rho_i, vis_i, cp_i, lambda_i, Re_i;

  area_i = PI * (_rio - _wt) * (_rio - _wt);
  u_i = _flow_i[_qp] / area_i;
  rho_i = eos_uo.rho_from_p_T(_p_i[_qp], _T_i[_qp]);
  vis_i = viscosity_uo.mu(_p_i[_qp], _T_i[_qp]);
  Re_i = rho_i * 2.0 * (_rio - _wt) * fabs(u_i) / vis_i;

  if (Re_i>0.0)
  {
     cp_i = eos_uo.cp(_p_i[_qp], _T_i[_qp]);
     lambda_i = eos_uo.lambda(_p_i[_qp], _T_i[_qp]);
     pr_i = vis_i * cp_i / lambda_i;
     if (Re_i<2300.0)
     {
         _nusselt_i[_qp]  = 4.364 ;
     }
     if (Re_i>2300.0 && Re_i<10000.0)
     {
         gama_i = (Re_i - 2300.0)/(10000.0 - 2300.0);
         _nusselt_i[_qp] = (1.0 - gama_i) * 4.364 + gama_i * 0.023 * pow(Re_i, 0.8) * pow(pr_i, 0.3);
         std::cout<<Re_i<< "   " << _nusselt_i[_qp] <<std::endl;
     }
     if (Re_i>10000.0)
     {
         _nusselt_i[_qp] = 0.023 * pow(Re_i, 0.8) * pow(pr_i, 0.3);
     }
  }
  return _nusselt_i[_qp] * lambda_i / 2.0 / (_rio - _wt);
}

Real
MoskitoCoaxialHeat_1p::Conv_coeff_outer()
{
  Real pr_o, gama_o, area_o, u_o, rho_o, vis_o, cp_o, lambda_o, Re_o ,hd;

  area_o = PI * (_roi * _roi - _rio * _rio );
  u_o = _flow_o[_qp] / area_o;
  rho_o = eos_uo.rho_from_p_T(_p_o[_qp], _T_o[_qp]);
  vis_o = viscosity_uo.mu(_p_o[_qp], _T_o[_qp]);
  cp_o = eos_uo.cp(_p_o[_qp], _T_o[_qp]);
  lambda_o = eos_uo.lambda(_p_o[_qp], _T_o[_qp]);
  hd = 4.0 * area_o / (2.0 * PI * (_roi + _rio));
  Re_o = rho_o * hd * fabs(u_o) / vis_o;

  if (Re_o>0.0)
  {
     pr_o = vis_o * cp_o / lambda_o;
     if (Re_o<2300.0)
     {
         _nusselt_o[_qp]  = 4.364 ;
     }
     if (2300.0<Re_o<10000.0)
     {
         gama_o = (Re_o - 2300.0)/(10000.0 - 2300.0);
         _nusselt_o[_qp] = (1.0 - gama_o) * 4.364 + gama_o * 0.023 * pow(Re_o, 0.8) * pow(pr_o, 0.3);
     }
     if (10000.0<Re_o)
     {
         _nusselt_o[_qp] = 0.023 * pow(Re_o, 0.8) * pow(pr_o, 0.3);
     }
  }

  return _nusselt_o[_qp] * lambda_o / hd;
}
