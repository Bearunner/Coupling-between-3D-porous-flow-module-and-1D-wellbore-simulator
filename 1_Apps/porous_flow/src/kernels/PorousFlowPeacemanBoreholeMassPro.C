//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPeacemanBoreholeMassPro.h"
#include "Function.h"

registerMooseObject("PorousFlowApp", PorousFlowPeacemanBoreholeMassPro);

InputParameters
PorousFlowPeacemanBoreholeMassPro::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<UserObjectName>(
       "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("porepressure", "porepressure");
  params.addRequiredCoupledVar("temperature", "temperature");
  // MooseEnum WT("production=-1 injection=1");
  // params.addRequiredParam<MooseEnum>(
  //     "well_type", WT, "production or injection");
  params.addParam<unsigned int>(
      "fluid_phase",
      0,
      "The fluid phase whose pressure (and potentially mobility, enthalpy, etc) "
      "controls the flux to the line sink.  For p_or_t=temperature, and without "
      "any use_*, this parameter is irrelevant");
  params.addRequiredCoupledVar("p_well", "p_well");
  // params.addRequiredCoupledVar("rho_well", "rho_well");
  // params.addRequiredCoupledVar("mu_well", "mu_well");
  params.addParam<Real>("WI", "constant");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("function", "1", "A function that describes the body force");
  params.addParam<PostprocessorName>("postprocessor", 1, "A postprocessor whose value is multiplied by the body force");
  params.declareControllable("value");
  params.addClassDescription(
      "Approximates a borehole in the mesh using the Peaceman approach, ie "
      "using a number of point sinks with given radii whose positions are "
      "read from a file.  NOTE: if you are using PorousFlowPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed");
  return params;
}

PorousFlowPeacemanBoreholeMassPro::PorousFlowPeacemanBoreholeMassPro(const InputParameters & parameters)
  : Kernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _p_var_number(coupled("porepressure")),
    _T_var_number(coupled("temperature")),
    _ph(getParam<unsigned int>("fluid_phase")),
    // _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
    _fluid_density_node(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")),
    _dfluid_density_node_dvar(getMaterialProperty<std::vector<std::vector<Real>>>(
        "dPorousFlow_fluid_phase_density_nodal_dvar")),
    _fluid_viscosity(getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_nodal")),
    _dfluid_viscosity_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_nodal_dvar")),
    _p_well(coupledValue("p_well")),
    // _rho_well(coupledValue("rho_well")),
    // _mu_well(coupledValue("mu_well")),
    // _well_type(getParam<MooseEnum>("well_type")),
    _WI(this->template getParam<Real>("WI")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor"))
{
}

// Real
// PorousFlowPeacemanBoreholeMassPro::ptqp() const
// {
//   return (_pp)[_qp][_ph];
// }

// Real
// PorousFlowPeacemanBoreholeMassPro::rhoqp() const
// {
//   return (_fluid_density)[_qp][_ph];
// }

Real
PorousFlowPeacemanBoreholeMassPro::computeQpResidual()
{
  // const Real pp = ptqp();
  // const Real rho = rhoqp();

    Real r = 0.0;

    for (unsigned int ph = 0; ph < _num_phases; ++ph)
    {

     r += _WI * _fluid_density_node[_i][_ph] / _fluid_viscosity[_i][_ph] * (_u[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));

    }

  return r * _test[_i][_qp];
}


Real
PorousFlowPeacemanBoreholeMassPro::computeQpJacobian()
{
  return computeQpJac(_var.number());
}


Real
PorousFlowPeacemanBoreholeMassPro::computeQpJac(unsigned int jvar) const
{
  // If the variable is not a valid PorousFlow variable, set the Jacobian to 0
  // if (_dictator.notPorousFlowVariable(jvar))
  //   return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);


  Real j = 0.0;

  if (jvar == _p_var_number)

  {
    for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {

  j += -_dfluid_viscosity_dvar[_i][ph][pvar] * _phi[_j][_qp] * _fluid_density_node[_i][ph] / std::pow(_fluid_viscosity[_i][ph], 2) + _dfluid_density_node_dvar[_i][ph][pvar] * _phi[_j][_qp] / _fluid_viscosity[_i][ph];
  j *= (_u[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));
  j += _fluid_density_node[_i][ph] / _fluid_viscosity[_i][ph] * _phi[_j][_qp];

  }

  }

  return _WI * j * _test[_i][_qp];
}

Real
PorousFlowPeacemanBoreholeMassPro::computeQpOffDiagJacobian(unsigned int jvar)
{

  // if (_dictator.notPorousFlowVariable(jvar))
  //   return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);

  Real j = 0.0;


  if (jvar == _T_var_number)

  {

  for (unsigned int ph = 0; ph < _num_phases; ++ph)

  {

  j += -_dfluid_viscosity_dvar[_i][ph][pvar] * _phi[_j][_qp] * _fluid_density_node[_i][ph] / std::pow(_fluid_viscosity[_i][ph], 2) + _dfluid_density_node_dvar[_i][ph][pvar] * _phi[_j][_qp] / _fluid_viscosity[_i][ph];
  j *= (_u[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));

  }

  }

  return _WI * j * _test[_i][_qp];
}
