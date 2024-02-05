//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowMassFractionS.h"

registerMooseObject("PorousFlowApp", PorousFlowMassFractionS);

InputParameters
PorousFlowMassFractionS::validParams()
{
  InputParameters params = PorousFlowMaterialVectorBase::validParams();
  params.addCoupledVar("temperature", "Fluid temperature variable.  Note, the default is suitable if your "
                       "simulation is using Kelvin units, but probably not for Celsius");
  params.addCoupledVar("C1", 0, "C1");
  params.addCoupledVar("C2", 0, "C2");
  params.addCoupledVar("C3", 0, "C3");
  params.addCoupledVar("C4", 0, "C4");
  params.addCoupledVar("C5", 0, "C5");
  params.addCoupledVar("C6", 0, "C6");
  params.addCoupledVar("C7", 0, "C7");
  params.addCoupledVar("C8", 0, "C8");
  params.addCoupledVar("C9", 0, "C9");
  params.addCoupledVar("C10", 0, "C10");
  params.addParam<Real>("C1_temperature_threshold", 1e99, "temperature_threshold for C1 temperature-response (K)");
  params.addParam<Real>("C3_temperature_threshold", 1e99, "temperature_threshold for C3 temperature-response (K)");
  params.addParam<Real>("C5_temperature_threshold", 1e99, "temperature_threshold for C5 temperature-response (K)");
  params.addParam<Real>("C7_temperature_threshold", 1e99, "temperature_threshold for C7 temperature-response (K)");
  params.addParam<Real>("C9_temperature_threshold", 1e99, "temperature_threshold for C9 temperature-response (K)");
  params.addCoupledVar(
      "mass_fraction_vars",
      "List of variables that represent the mass fractions.  Format is 'f_ph0^c0 "
      "f_ph0^c1 f_ph0^c2 ... f_ph0^c(N-2) f_ph1^c0 f_ph1^c1 fph1^c2 ... "
      "fph1^c(N-2) ... fphP^c0 f_phP^c1 fphP^c2 ... fphP^c(N-2)' where "
      "N=num_components and P=num_phases, and it is assumed that "
      "f_ph^c(N-1)=1-sum(f_ph^c,{c,0,N-2}) so that f_ph^c(N-1) need not be given.  If no "
      "variables are provided then num_phases=1=num_components.");
  params.addPrivateParam<std::string>("pf_material_type", "mass_fraction");
  // params.addParam<FunctionName>("initial_temperature", 0.0, "initial temperature)");
  // params.addParam<FunctionName>("initial_C1", 0.0, "initial C1)");
  params.addClassDescription("This Material forms a std::vector<std::vector ...> of mass-fractions "
                             "out of the individual mass fractions");
  return params;
}

PorousFlowMassFractionS::PorousFlowMassFractionS(const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),

    _mass_frac(_nodal_material
                   ? declareProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_nodal")
                   : declareProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_qp")),
    _grad_mass_frac(_nodal_material ? nullptr
                                    : &declareProperty<std::vector<std::vector<RealGradient>>>(
                                          "PorousFlow_grad_mass_frac_qp")),
    _dmass_frac_dvar(_nodal_material ? declareProperty<std::vector<std::vector<std::vector<Real>>>>(
                                           "dPorousFlow_mass_frac_nodal_dvar")
                                     : declareProperty<std::vector<std::vector<std::vector<Real>>>>(
                                           "dPorousFlow_mass_frac_qp_dvar")),
    _num_passed_mf_vars(coupledComponents("mass_fraction_vars")),
    _T(_nodal_material ? coupledDofValues("temperature")
                       : coupledValue("temperature")),
    _C1(_nodal_material ? coupledDofValues("C1")
                        : coupledValue("C1")),
    _C2(_nodal_material ? coupledDofValues("C2")
                        : coupledValue("C2")),
    _C3(_nodal_material ? coupledDofValues("C3")
                        : coupledValue("C3")),
    _C4(_nodal_material ? coupledDofValues("C4")
                        : coupledValue("C4")),
    _C5(_nodal_material ? coupledDofValues("C5")
                        : coupledValue("C5")),
    _C6(_nodal_material ? coupledDofValues("C6")
                        : coupledValue("C6")),
    _C7(_nodal_material ? coupledDofValues("C7")
                        : coupledValue("C7")),
    _C8(_nodal_material ? coupledDofValues("C8")
                        : coupledValue("C8")),
    _C9(_nodal_material ? coupledDofValues("C9")
                        : coupledValue("C9")),
    _C10(_nodal_material ? coupledDofValues("C10")
                        : coupledValue("C10")),
    _C1_T_threshold(getParam<Real>("C1_temperature_threshold")),
    _C3_T_threshold(getParam<Real>("C3_temperature_threshold")),
    _C5_T_threshold(getParam<Real>("C5_temperature_threshold")),
    _C7_T_threshold(getParam<Real>("C7_temperature_threshold")),
    _C9_T_threshold(getParam<Real>("C9_temperature_threshold")),
    _new_C1(_nodal_material
                      ? declareProperty<Real>("new_C1_nodal")
                      : declareProperty<Real>("new_C1_qp")),
    _new_C2(_nodal_material
                      ? declareProperty<Real>("new_C2_nodal")
                      : declareProperty<Real>("new_C2_qp")),
    _new_C3(_nodal_material
                      ? declareProperty<Real>("new_C3_nodal")
                      : declareProperty<Real>("new_C3_qp")),
    _new_C4(_nodal_material
                      ? declareProperty<Real>("new_C4_nodal")
                      : declareProperty<Real>("new_C4_qp")),
    _new_C5(_nodal_material
                      ? declareProperty<Real>("new_C5_nodal")
                      : declareProperty<Real>("new_C5_qp")),
    _new_C6(_nodal_material
                      ? declareProperty<Real>("new_C6_nodal")
                      : declareProperty<Real>("new_C6_qp")),
    _new_C7(_nodal_material
                      ? declareProperty<Real>("new_C7_nodal")
                      : declareProperty<Real>("new_C7_qp")),
    _new_C8(_nodal_material
                      ? declareProperty<Real>("new_C8_nodal")
                      : declareProperty<Real>("new_C8_qp")),
    _new_C9(_nodal_material
                      ? declareProperty<Real>("new_C9_nodal")
                      : declareProperty<Real>("new_C9_qp")),
    _new_C10(_nodal_material
                      ? declareProperty<Real>("new_C10_nodal")
                      : declareProperty<Real>("new_C10_qp")),
    _C1_indicator(_nodal_material
                      ? declareProperty<Real>("C1_indicator_nodal")
                      : declareProperty<Real>("C1_indicator_qp")),
    _C3_indicator(_nodal_material
                      ? declareProperty<Real>("C3_indicator_nodal")
                      : declareProperty<Real>("C3_indicator_qp")),
    _C5_indicator(_nodal_material
                      ? declareProperty<Real>("C5_indicator_nodal")
                      : declareProperty<Real>("C5_indicator_qp")),
    _C7_indicator(_nodal_material
                      ? declareProperty<Real>("C7_indicator_nodal")
                      : declareProperty<Real>("C7_indicator_qp")),
    _C9_indicator(_nodal_material
                      ? declareProperty<Real>("C9_indicator_nodal")
                      : declareProperty<Real>("C9_indicator_qp"))
                      // _initial_temperature(getFunction("initial_temperature")),
                      // _T_ini(_nodal_material
                      //                   ? declareProperty<Real>("initial_temperature_nodal")
                      //                   : declareProperty<Real>("initial_temperature_qp")),
                      // _initial_C1(getFunction("initial_C1")),
                      // _C1_ini(_nodal_material
                      //                   ? declareProperty<Real>("initial_C1_nodal")
                      //                   : declareProperty<Real>("initial_C1_qp")),
                      // _C1_indicator_ini(_nodal_material
                      //                   ? declareProperty<Real>("C1_indicator_ini_nodal")
                      //                   : declareProperty<Real>("C1_indicator_ini_qp")),
                      // _new_C1(declareProperty<Real>("new_C1")),
                      // _new_C2(declareProperty<Real>("new_C2")),
                      // _C1_indicator(declareProperty<Real>("C1_indicator")),
                      // _K(declareProperty<std::vector<std::vector<Real>>>("color_indicator"))
{
  if (_num_phases < 1 || _num_components < 1)
    mooseError("PorousFlowMassFractionS: The Dictator proclaims that the number of phases is ",
               _num_phases,
               " and the number of components is ",
               _num_components,
               ", and stipulates that you should not use PorousFlowMassFractionS in this case");

  if (_num_passed_mf_vars != _num_phases * (_num_components - 1))
    paramError("mass_fraction_vars",
               "This value must be equal to the Dictator's num_phases (",
               _num_phases,
               ") multiplied by num_components-1 (",
               _num_components - 1,
               ")");

  _mf_vars_num.resize(_num_passed_mf_vars);
  _mf_vars.resize(_num_passed_mf_vars);
  _grad_mf_vars.resize(_num_passed_mf_vars);
  for (unsigned i = 0; i < _num_passed_mf_vars; ++i)
  {
    _mf_vars_num[i] = coupled("mass_fraction_vars", i);
    _mf_vars[i] = (_nodal_material ? &coupledDofValues("mass_fraction_vars", i)
                                   : &coupledValue("mass_fraction_vars", i));
    _grad_mf_vars[i] = &coupledGradient("mass_fraction_vars", i);
    // _mf_vars[1] = (_nodal_material ? &coupledDofValues("mass_fraction_vars", 1)
    //                                : &coupledValue("mass_fraction_vars", 1));
    // _mf_vars[2] = (_nodal_material ? &coupledDofValues("mass_fraction_vars", 2)
    //                               : &coupledValue("mass_fraction_vars", 2));
  }
}

void
PorousFlowMassFractionS::initQpStatefulProperties()
{
    // _C1_ini[_qp] = _initial_C1.value(_t, _q_point[0]);
    // _T_ini[_qp] = _initial_temperature.value(_t, _q_point[0]);
    //
    // if (_T_ini[_qp]  < 368.15) {
    //   _C1_indicator_ini[_qp] = 1.0;
    // }
    // else if (_T_ini[_qp] >=368.15 && _T_ini[_qp] <=378.15){
    //   _C1_indicator_ini[_qp] = 1.0 /(1.0 + pow (2.71828182846, 10.0 * (_T_ini[_qp] - 373.15)));
    // }
    // else {
    //   _C1_indicator_ini[_qp] = 0.0;
    // }
    //
    // _new_C1[_qp] = _C1_indicator_ini[_qp] * _C1_ini[_qp];
    // _new_C2[_qp] = 0.0;

    // all we need to do is compute _mass_frac for _nodal_materials
    // but the following avoids code duplication
    computeQpProperties();
}

void
PorousFlowMassFractionS::computeQpProperties()
{
  if (_T[_qp]  < (_C1_T_threshold - 0.15)) {
    _C1_indicator[_qp] = 1.0;
  }
  else if (_T[_qp] >=(_C1_T_threshold - 0.15) && _T[_qp] <=(_C1_T_threshold + 0.15)){
    _C1_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C1_T_threshold)));
  }
  else {
    _C1_indicator[_qp] = 0.0;
  }

  if (_T[_qp]  < (_C3_T_threshold - 0.15)) {
    _C3_indicator[_qp] = 1.0;
  }
  else if (_T[_qp] >=(_C3_T_threshold - 0.15) && _T[_qp] <=(_C3_T_threshold + 0.15)){
    _C3_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C3_T_threshold)));
  }
  else {
    _C3_indicator[_qp] = 0.0;
  }

  if (_T[_qp]  < (_C5_T_threshold - 0.15)) {
    _C5_indicator[_qp] = 1.0;
  }
  else if (_T[_qp] >=(_C5_T_threshold - 0.15) && _T[_qp] <=(_C5_T_threshold + 0.15)){
    _C5_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C5_T_threshold)));
  }
  else {
    _C5_indicator[_qp] = 0.0;
  }

  if (_T[_qp]  < (_C7_T_threshold - 0.15)) {
    _C7_indicator[_qp] = 1.0;
  }
  else if (_T[_qp] >=(_C7_T_threshold - 0.15) && _T[_qp] <=(_C7_T_threshold + 0.15)){
    _C7_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C7_T_threshold)));
  }
  else {
    _C7_indicator[_qp] = 0.0;
  }

  if (_T[_qp]  < (_C9_T_threshold - 0.15)) {
    _C9_indicator[_qp] = 1.0;
  }
  else if (_T[_qp] >=(_C9_T_threshold - 0.15) && _T[_qp] <=(_C9_T_threshold + 0.15)){
    _C9_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C9_T_threshold)));
  }
  else {
    _C9_indicator[_qp] = 0.0;
  }

 _new_C1[_qp] = _C1_indicator[_qp] * _C1[_qp];
 _new_C2[_qp] = _C1[_qp] - _new_C1[_qp] + _C2[_qp];

 _new_C3[_qp] = _C3_indicator[_qp] * _C3[_qp];
 _new_C4[_qp] = _C3[_qp] - _new_C3[_qp] + _C4[_qp];

 _new_C5[_qp] = _C5_indicator[_qp] * _C5[_qp];
 _new_C6[_qp] = _C5[_qp] - _new_C5[_qp] + _C6[_qp];

 _new_C7[_qp] = _C7_indicator[_qp] * _C7[_qp];
 _new_C8[_qp] = _C7[_qp] - _new_C7[_qp] + _C8[_qp];

 _new_C9[_qp] = _C9_indicator[_qp] * _C9[_qp];
 _new_C10[_qp] = _C9[_qp] - _new_C9[_qp] + _C10[_qp];

  // size all properties correctly
  _mass_frac[_qp].resize(_num_phases);
  _dmass_frac_dvar[_qp].resize(_num_phases);
  if (!_nodal_material)
    (*_grad_mass_frac)[_qp].resize(_num_phases);
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _mass_frac[_qp][ph].resize(_num_components);
    _dmass_frac_dvar[_qp][ph].resize(_num_components);
    for (unsigned int comp = 0; comp < _num_components; ++comp)
      _dmass_frac_dvar[_qp][ph][comp].assign(_num_var, 0.0);
    if (!_nodal_material)
      (*_grad_mass_frac)[_qp][ph].resize(_num_components);
  }

  // compute the values and derivatives
  unsigned int i = 0;
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    Real total_mass_frac = 0;
    if (!_nodal_material)
      (*_grad_mass_frac)[_qp][ph][_num_components - 1] = 0.0;
    for (unsigned int comp = 0; comp < _num_components - 1; ++comp)
    {
      _mass_frac[_qp][ph][comp] = (*_mf_vars[i])[_qp];

      // _new_C1[_qp] = _C1_indicator[_qp] * (*_mf_vars[1])[_qp];
      // _new_C2[_qp] = (*_mf_vars[1])[_qp] - _new_C1[_qp] + (*_mf_vars[2])[_qp];
      // if (_T[_qp] < 368.15) {
      //   _C1_indicator[_qp][ph][comp] = 1.0;
      // }
      // else if (_T[_qp] >=368.15 && _T[_qp] <=378.15){
      //   _C1_indicator[_qp][ph][comp] = 1.0 /(1.0 + pow (2.71828182846, 10.0 * (_T[_qp] - 373.15)));
      // }
      // else {
      //   _C1_indicator[_qp][ph][comp] = 0.0;
      // }

      // _new_C1[_qp] = _C1_indicator[_qp][ph][comp] * (*_mf_vars[1])[_qp];
      // _new_C2[_qp] = (*_mf_vars[1])[_qp] - _new_C1[_qp] + (*_mf_vars[2])[_qp];

      total_mass_frac += _mass_frac[_qp][ph][comp];
      if (!_nodal_material)
      {
        (*_grad_mass_frac)[_qp][ph][comp] = (*_grad_mf_vars[i])[_qp];
        (*_grad_mass_frac)[_qp][ph][_num_components - 1] -= (*_grad_mf_vars[i])[_qp];
      }
      if (_dictator.isPorousFlowVariable(_mf_vars_num[i]))
      {
        // _mf_vars[i] is a PorousFlow variable
        const unsigned int pf_var_num = _dictator.porousFlowVariableNum(_mf_vars_num[i]);
        _dmass_frac_dvar[_qp][ph][comp][pf_var_num] = 1.0;
        _dmass_frac_dvar[_qp][ph][_num_components - 1][pf_var_num] = -1.0;
      }
      i++;
    }
    _mass_frac[_qp][ph][_num_components - 1] = 1.0 - total_mass_frac;
  }

}
