
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

// class MoskitoMomentum_2p1c;
//
// template <>
// InputParameters validParams<MoskitoMomentum_2p1c>();

class MoskitoMomentum_2p1c : public Kernel
{
public:
  static InputParameters validParams();
  MoskitoMomentum_2p1c(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // The required values for massrate coupling
  const VariableValue & _m;
  const VariableGradient & _grad_m;
  unsigned _m_var_number;

  // The required values for enthalpy coupling
  const VariableGradient & _grad_h;
  unsigned _h_var_number;

  const VariableGradient & _grad_c;
  unsigned _c_var_number;

  // Density derivatives wrt pressure and enthalpy
  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _drho_dp;
  const MaterialProperty<Real> & _drho_dp2;
  const MaterialProperty<Real> & _drho_dh;
  const MaterialProperty<Real> & _drho_dh2;
  const MaterialProperty<Real> & _drho_dph;

  const MaterialProperty<Real> & _drho_dc;
  const MaterialProperty<RealVectorValue> & _drho_m_dcvect;
  const MaterialProperty<RealVectorValue> & _dgamma_dcvect;

  // The pipe Moody friction factor
  const MaterialProperty<Real> & _f;
  // The gravity acceleration as a vector
  const MaterialProperty<RealVectorValue> & _gravity;
  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The wetted perimeter of pipe
  const MaterialProperty<Real> & _perimeter;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
  // The flow direction
  const MaterialProperty<Real> & _well_sign;
  //
  // // The gamma derivatives
  const MaterialProperty<Real> & _dgamma_dp;
  const MaterialProperty<Real> & _dgamma_dh;
  const MaterialProperty<Real> & _dgamma_dm;
  const MaterialProperty<Real> & _dgamma_dc;
};
