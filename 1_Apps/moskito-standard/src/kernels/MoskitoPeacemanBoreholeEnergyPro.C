//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MoskitoPeacemanBoreholeEnergyPro.h"
#include "Function.h"

registerMooseObject("MoskitoApp", MoskitoPeacemanBoreholeEnergyPro);

InputParameters
MoskitoPeacemanBoreholeEnergyPro::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  // params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable");
  params.addRequiredCoupledVar("p_res", "p_res");
  // params.addRequiredCoupledVar("T_res", "T_res");
  params.addRequiredCoupledVar("rho_res", "rho_res");
  params.addRequiredCoupledVar("mu_res", "mu_res");
  params.addRequiredCoupledVar("h_res", "h_res");
  // MooseEnum WT("production=-1 injection=1");
  // params.addRequiredParam<MooseEnum>(
  //     "well_type", WT, "production or injection");
  params.addParam<Real>("WI", "constant");
  // params.addParam<Real>("lambda_oh", "thermal resistivity of openhole section");
  // params.addParam<Real>("radius_ot", "outside radius of tubing");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", "1", "A function that describes the body force");
  params.addParam<PostprocessorName>("postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  // params.addParam<Real>("value_T", 1.0, "Coefficient to multiply by the body force term");
  // params.addParam<FunctionName>("function_T", "1", "A function that describes the body force");
  // params.addParam<PostprocessorName>("postprocessor_T", 1, "A postprocessor whose value is multiplied by the body force");
  params.declareControllable("value");
  // params.declareControllable("value_T");
  params.addClassDescription(
      "Approximates a borehole in the mesh using the Peaceman approach, ie "
      "using a number of point sinks with given radii whose positions are "
      "read from a file.  NOTE: if you are using MoskitoPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed");
  return params;
}

MoskitoPeacemanBoreholeEnergyPro::MoskitoPeacemanBoreholeEnergyPro(const InputParameters & parameters)
  : Kernel(parameters),
    _p(coupledValue("pressure")),
    _p_res(coupledValue("p_res")),
    // _T_res(coupledValue("T_res")),
    _rho_res(coupledValue("rho_res")),
    _mu_res(coupledValue("mu_res")),
    _h_res(coupledValue("h_res")),
    _q_var_number(coupled("flowrate")),
    _p_var_number(coupled("pressure")),
    // _rho(getMaterialProperty<Real>("density")),
    // _drho_dp(getMaterialProperty<Real>("drho_dp")),
    // _drho_dT(getMaterialProperty<Real>("drho_dT")),
    // _mu(getMaterialProperty<Real>("viscosity")),
    // _dmu_dp(getMaterialProperty<Real>("dmu_dp")),
    // _dmu_dT(getMaterialProperty<Real>("dmu_dT")),
    // _well_type(getParam<MooseEnum>("well_type")),
    // _area(getMaterialProperty<Real>("well_area")),
    _WI(this->template getParam<Real>("WI")),
    _dia(getMaterialProperty<Real>("well_diameter")),
    // _lambda_oh(this->template getParam<Real>("lambda_oh")),
    // _radius_ot(this->template getParam<Real>("radius_ot")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor"))
    // _scale_T(this->template getParam<Real>("value_T")),
    // _function_T(getFunction("function_T")),
    // _postprocessor_T(getPostprocessorValue("postprocessor_T"))
{
}

Real
MoskitoPeacemanBoreholeEnergyPro::computeQpResidual()
{
  Real r = 0.0;

  r += _WI * _h_res[_qp] * _rho_res[_qp] / _mu_res[_qp] * (_p[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));

  //kinetic terms
   // r += std::pow(0.5 * _dia[_qp], 2)/8 * std::pow(_WI, 3) * _rho_res[_qp] / std::pow(_mu_res[_qp], 3) * std::pow((_p[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])), 3);

  // conduction terms
   // r -=  2.0 * PI * _radius_ot * _lambda_oh * (_T_res[_qp] * _scale_T * _postprocessor_T * _function_T.value(_t, _q_point[_qp]) - _u[_qp]) / _area[_qp];

  return r * _test[_i][_qp];

}


Real
MoskitoPeacemanBoreholeEnergyPro::computeQpJacobian()
{
  Real j = 0.0;

  // conduction terms

  // j +=  2.0 * PI * _radius_ot * _lambda_oh  * _phi[_j][_qp] / _area[_qp];

  return j * _test[_i][_qp];
}

Real
MoskitoPeacemanBoreholeEnergyPro::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _q_var_number)
  {
    j += 0;
  }

  if (jvar == _p_var_number)
  {

    j += _WI * _h_res[_qp] * _rho_res[_qp] / _mu_res[_qp] * _phi[_j][_qp];

    //kinetic terms
    // j += 3 * std::pow(0.5 * _dia[_qp], 2)/8 * std::pow(_WI, 3) * _rho_res[_qp] / std::pow(_mu_res[_qp], 3) * std::pow((_p[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])), 2) * _phi[_j][_qp];

  }

  return j * _test[_i][_qp];
}
