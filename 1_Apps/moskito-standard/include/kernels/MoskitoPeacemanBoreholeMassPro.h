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
class MoskitoPeacemanBoreholeMassPro : public Kernel
{
public:
  /**
   * Creates a new MoskitoPeacemanBoreholeMassPro
   * This reads the file containing the lines of the form
   * radius x y z
   * that defines the borehole geometry.
   * It also calculates segment-lengths and rotation matrices
   * needed for computing the borehole well constant
   */
  static InputParameters validParams();

  MoskitoPeacemanBoreholeMassPro(const InputParameters & parameters);

  // void virtual initialSetup() override;


protected:

virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

// const VariableValue & _p;
const VariableValue & _p_res;
const VariableValue & _rho_res;
const VariableValue & _mu_res;

// Variable numberings
unsigned _q_var_number;
unsigned _T_var_number;

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

// MooseEnum _well_type;
const Real & _WI;
  /// Scale factor
const Real & _scale;

/// Optional function value
const Function & _function;

/// Optional Postprocessor value
const PostprocessorValue & _postprocessor;

  /// Whether there is a quadpoint permeability material (for error checking)
  // const bool _has_permeability;

  /// Whether there is a quadpoint thermal conductivity material (for error checking)
  // const bool _has_thermal_conductivity;

  /// Permeability or conductivity of porous material

  // Real computeQpBaseOutflow() const override;

};
