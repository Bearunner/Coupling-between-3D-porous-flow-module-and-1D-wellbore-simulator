//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MoskitoPeacemanBoreholeMassPro.h"
#include "Function.h"

registerMooseObject("MoskitoApp", MoskitoPeacemanBoreholeMassPro);

InputParameters
MoskitoPeacemanBoreholeMassPro::validParams()
{
  InputParameters params = Kernel::validParams();
  // params.addRequiredCoupledVar("pressure", "Pressure nonlinear variable");
  params.addRequiredCoupledVar("flowrate", "Volumetric flowrate nonlinear variable");
  params.addRequiredCoupledVar("temperature", "Temperature nonlinear variable");
  params.addRequiredCoupledVar("p_res", "p_res");
  params.addRequiredCoupledVar("rho_res", "rho_res");
  params.addRequiredCoupledVar("mu_res", "mu_res");
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

MoskitoPeacemanBoreholeMassPro::MoskitoPeacemanBoreholeMassPro(const InputParameters & parameters)
  : Kernel(parameters),
    _p_res(coupledValue("p_res")),
    _rho_res(coupledValue("rho_res")),
    _mu_res(coupledValue("mu_res")),
    _q_var_number(coupled("flowrate")),
    _T_var_number(coupled("temperature")),
    _WI(this->template getParam<Real>("WI")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor"))
{
}

Real
MoskitoPeacemanBoreholeMassPro::computeQpResidual()
{

  return _WI * _rho_res[_qp] / _mu_res[_qp] * (_u[_qp] - _p_res[_qp] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])) * _test[_i][_qp];

}


Real
MoskitoPeacemanBoreholeMassPro::computeQpJacobian()
{
  Real j = 0.0;

  j += _rho_res[_qp] / _mu_res[_qp] * _phi[_j][_qp];

  return _WI * j * _test[_i][_qp];
}

Real
MoskitoPeacemanBoreholeMassPro::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _q_var_number)
  {
    j += 0;
  }

  if (jvar == _T_var_number)
  {

    j += 0;
  }

  return j * _test[_i][_qp];
}
