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
// #include "MoskitoLineSink.h"

/**
 * Approximates a borehole by a sequence of Dirac Points
 */
class MoskitoPeacemanBoreholeEnergyPro : public Kernel
{
public:
  /**
   * Creates a new MoskitoPeacemanBoreholeEnergyPro
   * This reads the file containing the lines of the form
   * radius x y z
   * that defines the borehole geometry.
   * It also calculates segment-lengths and rotation matrices
   * needed for computing the borehole well constant
   */
  static InputParameters validParams();

  MoskitoPeacemanBoreholeEnergyPro(const InputParameters & parameters);

  // void virtual initialSetup() override;

protected:

virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

const VariableValue & _p;
const VariableValue & _p_res;
// const VariableValue & _T_res;
const VariableValue & _rho_res;
const VariableValue & _mu_res;
const VariableValue & _h_res;
// Variable numberings
unsigned _q_var_number;
unsigned _p_var_number;

// const MaterialProperty<Real> & _rho;
//
// const MaterialProperty<Real> & _drho_dp;
// // The first derivative of density wrt temperature
// const MaterialProperty<Real> & _drho_dT;
//
// const MaterialProperty<Real> & _mu;
//
// const MaterialProperty<Real> & _dmu_dp;
// // The first derivative of density wrt temperature
// const MaterialProperty<Real> & _dmu_dT;
// const MaterialProperty<Real> & _area;
// MooseEnum _well_type;
const Real & _WI;
const MaterialProperty<Real> & _dia;

// const Real & _lambda_oh;
// const Real & _radius_ot;
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
