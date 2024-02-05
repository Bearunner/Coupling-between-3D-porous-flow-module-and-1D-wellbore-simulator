/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#pragma once

#include "Kernel.h"

class MoskitoEnergy_1p1c : public Kernel
{
public:
  static InputParameters validParams();

  MoskitoEnergy_1p1c(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // The coupled flow_rate
  const VariableValue & _q;
  // The gradient of the coupled flow_rate
  // const VariableGradient & _grad_q;
  const VariableValue & _p;
  // The gradient of the coupled pressure
  // const VariableGradient & _grad_p;

  // Variable numberings
  unsigned _q_var_number;
  unsigned _p_var_number;

  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
  // The sign of well flow direction
  const MaterialProperty<Real> & _well_sign;
  // The thermal conductivity of casing and fluid
  // const MaterialProperty<Real> & _lambda;
  // The specific heat at constant pressure
  const MaterialProperty<Real> & _cp;
  // The density
  const MaterialProperty<Real> & _rho;
  // The first derivative of density wrt pressure
  const MaterialProperty<Real> & _drho_dp;
  // The first derivative of density wrt enthalpy
  const MaterialProperty<Real> & _drho_dT;
  // Enthalpy from P and T
  const MaterialProperty<Real> & _h;
  // The gravity acceleration as a vector
  const bool & _add_g;
  Real _gfac;
  const MaterialProperty<RealVectorValue> & _gravity;
};
