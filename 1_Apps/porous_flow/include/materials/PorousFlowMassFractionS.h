//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowMaterialVectorBase.h"
#include "Function.h"

/**
 * Material designed to form a std::vector<std::vector>
 * of mass fractions from the individual mass fraction variables
 */
class PorousFlowMassFractionS : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  PorousFlowMassFractionS(const InputParameters & parameters);

protected:
  /// Mass fraction matrix at quadpoint or nodes
  MaterialProperty<std::vector<std::vector<Real>>> & _mass_frac;

  /// Gradient of the mass fraction matrix at the quad points
  MaterialProperty<std::vector<std::vector<RealGradient>>> * const _grad_mass_frac;

  /// Derivative of the mass fraction matrix with respect to the porous flow variables
  MaterialProperty<std::vector<std::vector<std::vector<Real>>>> & _dmass_frac_dvar;

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /**
   * Builds the mass-fraction variable matrix at the quad point
   * @param qp the quad point
   */
  void build_mass_frac(unsigned int qp);

  /**
   * Number of mass-fraction variables provided by the user
   * This needs to be _num_phases*(_num_components - 1), since the
   * mass fraction of the final component in each phase is
   * determined as 1 - sum_{components}(mass fraction of all other components in the phase)
   */
  const unsigned int _num_passed_mf_vars;

  const VariableValue & _T;
  const VariableValue & _C1;
  const VariableValue & _C2;
  const VariableValue & _C3;
  const VariableValue & _C4;
  const VariableValue & _C5;
  const VariableValue & _C6;
  const VariableValue & _C7;
  const VariableValue & _C8;
  const VariableValue & _C9;
  const VariableValue & _C10;
  Real _C1_T_threshold;
  Real _C3_T_threshold;
  Real _C5_T_threshold;
  Real _C7_T_threshold;
  Real _C9_T_threshold;
  MaterialProperty<Real> & _new_C1;
  MaterialProperty<Real> & _new_C2;
  MaterialProperty<Real> & _new_C3;
  MaterialProperty<Real> & _new_C4;
  MaterialProperty<Real> & _new_C5;
  MaterialProperty<Real> & _new_C6;
  MaterialProperty<Real> & _new_C7;
  MaterialProperty<Real> & _new_C8;
  MaterialProperty<Real> & _new_C9;
  MaterialProperty<Real> & _new_C10;
  MaterialProperty<Real> & _C1_indicator;
  MaterialProperty<Real> & _C3_indicator;
  MaterialProperty<Real> & _C5_indicator;
  MaterialProperty<Real> & _C7_indicator;
  MaterialProperty<Real> & _C9_indicator;
  // const Function & _initial_temperature;
  // MaterialProperty<Real> & _T_ini;
  // const Function & _initial_C1;
  // MaterialProperty<Real> & _C1_ini;
  // MaterialProperty<Real> & _C1_indicator_ini;
  // MaterialProperty<std::vector<std::vector<Real>>> & _C1_indicator;
  // MaterialProperty<std::vector<std::vector<Real>>> & _K;

  /// The variable number of the mass-fraction variables
  std::vector<unsigned int> _mf_vars_num;

  /// The mass-fraction variables
  std::vector<const VariableValue *> _mf_vars;

  /// The gradient of the mass-fraction variables
  std::vector<const VariableGradient *> _grad_mf_vars;
  // std::vector<const VariableValue *> _mf_vars[1];
  // std::vector<const VariableValue *> _mf_vars[2];

};
