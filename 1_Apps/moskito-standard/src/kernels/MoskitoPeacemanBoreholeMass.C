//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MoskitoPeacemanBoreholeMass.h"
#include "Function.h"

registerMooseObject("MoskitoApp", MoskitoPeacemanBoreholeMass);

InputParameters
MoskitoPeacemanBoreholeMass::validParams()
{
  InputParameters params = Kernel::validParams();

  // params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("p_res", "p_res");
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable");

  // MooseEnum WT("production=-1 injection=1");
  // params.addRequiredParam<MooseEnum>(
  //     "well_type", WT, "production or injection");

  params.addParam<Real>("WI", "constant");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", "1", "A function that describes the body force");
  params.addParam<PostprocessorName>("postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  params.declareControllable("value");
  params.addClassDescription(
      "Approximates a borehole in the mesh using the Peaceman approach, ie "
      "using a number of point sinks with given radii whose positions are "
      "read from a file.  NOTE: if you are using MoskitoPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed");
  return params;
}

MoskitoPeacemanBoreholeMass::MoskitoPeacemanBoreholeMass(const InputParameters & parameters)
  : Kernel(parameters),
    // _p(coupledValue("pressure")),
    _p_res(coupledValue("p_res")),
    _q_var_number(coupled("flowrate")),
    _T_var_number(coupled("temperature")),
    // _p_var_number(coupled("pressure")),
    _rho(getMaterialProperty<Real>("density")),
    _drho_dp(getMaterialProperty<Real>("drho_dp")),
    _drho_dT(getMaterialProperty<Real>("drho_dT")),
    _mu(getMaterialProperty<Real>("viscosity")),
    _dmu_dp(getMaterialProperty<Real>("dmu_dp")),
    _dmu_dT(getMaterialProperty<Real>("dmu_dT")),
    // _well_type(getParam<MooseEnum>("well_type")),
    _WI(this->template getParam<Real>("WI")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor"))
{
}

Real
MoskitoPeacemanBoreholeMass::computeQpResidual()
{

  return _WI * _rho[_qp] / _mu[_qp] * (_u[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])) * _test[_i][_qp];

}


Real
MoskitoPeacemanBoreholeMass::computeQpJacobian()
{
  Real j = 0.0;

  j += -_dmu_dp[_qp] * _phi[_j][_qp] * _rho[_qp] / (_mu[_qp] * _mu[_qp]) + _drho_dp[_qp] * _phi[_j][_qp] / _mu[_qp];
  j *= (_u[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));
  j += _rho[_qp] / _mu[_qp] * _phi[_j][_qp];

  return _WI * j * _test[_i][_qp];
}

Real
MoskitoPeacemanBoreholeMass::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _q_var_number)
  {

    j += 0;

  }

  if (jvar == _T_var_number)
  {

    j += -_dmu_dT[_qp] * _phi[_j][_qp] * _rho[_qp] / (_mu[_qp] * _mu[_qp]) + _drho_dT[_qp] * _phi[_j][_qp] / _mu[_qp];
    j *= (_u[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));

  }

  return _WI * j * _test[_i][_qp];
}
