//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "PorousFlowDictator.h"
// #include "PorousFlowLineSink.h"

/**
 * Approximates a borehole by a sequence of Dirac Points
 */
class PorousFlowPeacemanBoreholeEnergyPro : public Kernel
{
public:
  /**
   * Creates a new PorousFlowPeacemanBoreholeEnergyPro
   * This reads the file containing the lines of the form
   * radius x y z
   * that defines the borehole geometry.
   * It also calculates segment-lengths and rotation matrices
   * needed for computing the borehole well constant
   */
  static InputParameters validParams();

  PorousFlowPeacemanBoreholeEnergyPro(const InputParameters & parameters);

protected:
  // void virtual initialSetup() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // virtual Real h(unsigned nodenum, unsigned phase) const;
  // virtual Real dh(unsigned nodenum, unsigned phase, unsigned pvar) const;

  Real computeQpJac(unsigned int jvar) const;
  // Real ptqp() const;

  const PorousFlowDictator & _dictator;
  const unsigned int _num_phases;
  unsigned _p_var_number;
  unsigned _T_var_number;

  const unsigned int _ph;
  // const MaterialProperty<std::vector<Real>> & _pp;
  const VariableValue & _pp;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _fluid_density_node;

  /// Derivative of the fluid density for each phase wrt PorousFlow variables (at the qp)
  const MaterialProperty<std::vector<std::vector<Real>>> & _dfluid_density_node_dvar;
  /// PorousFlowDictator UserObject
  const MaterialProperty<std::vector<Real>> & _fluid_viscosity;

  /// Derivative of viscosity wrt PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dfluid_viscosity_dvar;

  const MaterialProperty<std::vector<Real>> & _enthalpy;

  /// Derivative of the enthalpy wrt PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _denthalpy_dvar;
// Real rhoqp() const;

// const MaterialProperty<std::vector<Real>> & _fluid_density;

const VariableValue & _p_well;
// const VariableValue & _T_well;
// const VariableValue & _rho_well;
// const VariableValue & _mu_well;
// MooseEnum _well_type;
const Real & _WI;
// const Real & _dx;
// const Real & _dy;
// const Real & _lambda_oh;
// const Real & _radius_ot;
// const Real & _well_radius;
  /// Scale factor
const Real & _scale;
/// Optional function value
const Function & _function;
/// Optional Postprocessor value
const PostprocessorValue & _postprocessor;
// const Real & _scale_T;
/// Optional function value
// const Function & _function_T;
/// Optional Postprocessor value
// const PostprocessorValue & _postprocessor_T;
// const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;

  /// Whether there is a quadpoint permeability material (for error checking)
  // const bool _has_permeability;

  /// Whether there is a quadpoint thermal conductivity material (for error checking)
  // const bool _has_thermal_conductivity;

  /// Permeability or conductivity of porous material

  // Real computeQpBaseOutflow() const override;

};
