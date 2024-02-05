//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPeacemanBoreholeEnergyPro.h"
#include "Function.h"

registerMooseObject("PorousFlowApp", PorousFlowPeacemanBoreholeEnergyPro);

InputParameters
PorousFlowPeacemanBoreholeEnergyPro::validParams()
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
  params.addRequiredCoupledVar("pp", "pp");
  params.addRequiredCoupledVar("p_well", "p_well");
  // params.addRequiredCoupledVar("T_well", "T_well");
  // params.addRequiredCoupledVar("rho_well", "rho_well");
  // params.addRequiredCoupledVar("mu_well", "mu_well");
  params.addParam<Real>("WI", "constant");
  // params.addParam<Real>("dx", "dx");
  // params.addParam<Real>("dy", "dy");
  // params.addParam<Real>("lambda_oh", "thermal resistivity of openhole section");
  // params.addParam<Real>("radius_ot", "outside radius of tubing");
  // params.addParam<Real>("well_radius", "constant");
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
      "read from a file.  NOTE: if you are using PorousFlowPorosity that depends on volumetric "
      "strain, you should set strain_at_nearest_qp=true in your GlobalParams, to ensure the "
      "Porosity Material uses the volumetric strain at the Dirac quadpoints, and can therefore be "
      "computed");
  return params;
}

PorousFlowPeacemanBoreholeEnergyPro::PorousFlowPeacemanBoreholeEnergyPro(const InputParameters & parameters)
  : Kernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _p_var_number(coupled("porepressure")),
    _T_var_number(coupled("temperature")),
    _ph(getParam<unsigned int>("fluid_phase")),
    _pp(coupledValue("pp")), // or material???
    // _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_nodal")),
    _fluid_density_node(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")),
    _dfluid_density_node_dvar(getMaterialProperty<std::vector<std::vector<Real>>>(
        "dPorousFlow_fluid_phase_density_nodal_dvar")),
    _fluid_viscosity(getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_nodal")),
    _dfluid_viscosity_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_nodal_dvar")),
    _enthalpy(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_enthalpy_nodal")),
    _denthalpy_dvar(getMaterialProperty<std::vector<std::vector<Real>>>(
        "dPorousFlow_fluid_phase_enthalpy_nodal_dvar")),
    _p_well(coupledValue("p_well")),
    // _T_well(coupledValue("T_well")),
    // _rho_well(coupledValue("rho_well")),
    // _mu_well(coupledValue("mu_well")),
    // _well_type(getParam<MooseEnum>("well_type")),
    _WI(this->template getParam<Real>("WI")),
    // _dx(this->template getParam<Real>("dx")),
    // _dy(this->template getParam<Real>("dy")),
    // _lambda_oh(this->template getParam<Real>("lambda_oh")),
    // _radius_ot(this->template getParam<Real>("radius_ot")),
    // _well_radius(this->template getParam<Real>("well_radius")),
    _scale(this->template getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor"))
    // _scale_T(this->template getParam<Real>("value")),
    // _function_T(getFunction("function")),
    // _postprocessor_T(getPostprocessorValue("postprocessor"))
{
}

// Real
// PorousFlowPeacemanBoreholeEnergyPro::ptqp() const
// {
//   return (_pp)[_qp][_ph];
// }

// Real
// PorousFlowPeacemanBoreholeEnergyPro::rhoqp() const
// {
//   return (_fluid_density)[_qp][_ph];
// }

Real
PorousFlowPeacemanBoreholeEnergyPro::computeQpResidual()
{
  // const Real pp = ptqp();
  // const unsigned int num_nodes = _test.size();

    Real r = 0.0;

    for (unsigned int ph = 0; ph < _num_phases; ++ph)
    {
      // for (unsigned n = 0; n < num_nodes; ++n)
  {
     r += _WI * _enthalpy[_i][_ph] *  _fluid_density_node[_i][_ph] / _fluid_viscosity[_i][_ph] * (_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));

      //kinetic term

     // r += std::pow(_well_radius, 2)/8 * std::pow(_WI, 3) * _fluid_density_node[_i][_ph] / std::pow(_fluid_viscosity[_i][_ph], 3) * std::pow((_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])), 3);

     // conduction term

     // r -=  2.0 * PI * _radius_ot * _lambda_oh * (_T_well[_qp] * _scale_T * _postprocessor_T * _function_T.value(_t, _q_point[_qp]) - _u[_qp]) / (_dx * _dy);

}
    }

  return  r * _test[_i][_qp];
}


Real
PorousFlowPeacemanBoreholeEnergyPro::computeQpJacobian()
{
  return computeQpJac(_var.number());
}

Real
PorousFlowPeacemanBoreholeEnergyPro::computeQpJac(unsigned int jvar) const
{
  // const Real pp = ptqp();

  // If the variable is not a valid PorousFlow variable, set the Jacobian to 0

  // if (_dictator.notPorousFlowVariable(jvar))
  //   return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  const unsigned int num_nodes = _test.size();


  Real j = 0.0;


  if (jvar == _T_var_number)

  {

    for (unsigned int ph = 0; ph < _num_phases; ++ph)

  {
    for (unsigned n = 0; n < num_nodes; ++n)
  {

  j += _denthalpy_dvar[_i][_ph][pvar] * _fluid_density_node[_i][_ph] / _fluid_viscosity[_i][_ph] * _phi[_j][_qp];
  j += _enthalpy[_i][_ph] * (_dfluid_density_node_dvar[_i][_ph][pvar] * _phi[_j][_qp] / _fluid_viscosity[_i][_ph] - _dfluid_viscosity_dvar[_i][_ph][pvar] * _phi[_j][_qp] * _fluid_density_node[_i][_ph] / std::pow(_fluid_viscosity[_i][_ph], 2));
  j *= _WI * (_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));

  //kinetic term

  // j += std::pow(_well_radius, 2)/8 * std::pow(_WI, 3) * (_dfluid_density_node_dvar[_i][_ph][pvar] * _phi[_j][_qp] / std::pow(_fluid_viscosity[_i][_ph], 3) - 3 * _fluid_density_node[_i][_ph] / std::pow(_fluid_viscosity[_i][_ph], 4) * _dfluid_viscosity_dvar[_i][_ph][pvar] * _phi[_j][_qp]) * std::pow((_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])), 3);

  // conduction term

  // j +=  2.0 * PI * _radius_ot * _lambda_oh  * _phi[_j][_qp] / (_dx * _dy);

  }
  }

  }

  return  j * _test[_i][_qp];
}

Real
PorousFlowPeacemanBoreholeEnergyPro::computeQpOffDiagJacobian(unsigned int jvar)
{

  // const Real pp = ptqp();

  // If the variable is not a valid PorousFlow variable, set the Jacobian to 0

  // if (_dictator.notPorousFlowVariable(jvar))
  //   return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  const unsigned int num_nodes = _test.size();


  Real j = 0.0;

  if (jvar == _p_var_number)

  {

    for (unsigned int ph = 0; ph < _num_phases; ++ph)

  {
    for (unsigned n = 0; n < num_nodes; ++n)
  {

  j += _denthalpy_dvar[_i][_ph][pvar] * _fluid_density_node[_i][_ph] / _fluid_viscosity[_i][_ph] * _phi[_j][_qp];
  j += _enthalpy[_i][_ph] * (_dfluid_density_node_dvar[_i][_ph][pvar] * _phi[_j][_qp] / _fluid_viscosity[_i][_ph] - _dfluid_viscosity_dvar[_i][_ph][pvar] * _phi[_j][_qp] * _fluid_density_node[_i][_ph] / std::pow(_fluid_viscosity[_i][_ph], 2));
  j *= (_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp]));
  j += _enthalpy[_i][_ph] * _fluid_density_node[_i][_ph] / _fluid_viscosity[_i][_ph] * _phi[_j][_qp];
  j *= _WI;

   //kinetic term

  // j += std::pow(_well_radius, 2)/8 * std::pow(_WI, 3) * ((_dfluid_density_node_dvar[_i][_ph][pvar] * _phi[_j][_qp] / std::pow(_fluid_viscosity[_i][_ph], 3)) - 3 * _fluid_density_node[_i][_ph] / std::pow(_fluid_viscosity[_i][_ph], 4) * _dfluid_viscosity_dvar[_i][_ph][pvar] * _phi[_j][_qp]) * std::pow((_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])), 3);
  // j += std::pow(_well_radius, 2)/8 * std::pow(_WI, 3) * 3 * _fluid_density_node[_i][_ph] / std::pow(_fluid_viscosity[_i][_ph], 3) * _phi[_j][_qp] * std::pow((_pp[_i] - _p_well[_i] * _scale * _postprocessor * _function.value(_t, _q_point[_qp])), 2);

  }
  }
  }

  return  j * _test[_i][_qp];

}


// Real
// PorousFlowPeacemanBoreholeEnergyPro::h(unsigned nodenum, unsigned phase) const
// {
//   return _enthalpy[nodenum][phase];
// }
//
// Real
// PorousFlowPeacemanBoreholeEnergyPro::dh(unsigned nodenum, unsigned phase, unsigned pvar) const
// {
//   Real dh = _denthalpy_dvar[nodenum][phase][pvar];
//   return dh;
// }
