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
class PorousFlowPeacemanBoreholeEnergy : public Kernel
{
public:
  /**
   * Creates a new PorousFlowPeacemanBoreholeEnergy
   * This reads the file containing the lines of the form
   * radius x y z
   * that defines the borehole geometry.
   * It also calculates segment-lengths and rotation matrices
   * needed for computing the borehole well constant
   */
  static InputParameters validParams();

  PorousFlowPeacemanBoreholeEnergy(const InputParameters & parameters);

protected:
  // void virtual initialSetup() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // Real computeQpJac(unsigned int jvar) const;

  const PorousFlowDictator & _dictator;
  // const unsigned int _num_phases;
// Real ptqp() const;
// Real rhoqp() const;
// const unsigned int _ph;

// const MaterialProperty<std::vector<Real>> & _pp;
// const MaterialProperty<std::vector<Real>> & _fluid_density;
unsigned _p_var_number;
unsigned _T_var_number;
const VariableValue & _pp;
const VariableValue & _p_well;
// const VariableValue & _T_well;
const VariableValue & _rho_well;
const VariableValue & _mu_well;
const VariableValue & _h_well;
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

};
