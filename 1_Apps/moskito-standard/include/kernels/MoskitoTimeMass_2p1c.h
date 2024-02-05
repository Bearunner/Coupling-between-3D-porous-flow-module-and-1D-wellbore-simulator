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

#include "TimeKernel.h"

// class MoskitoTimeMass_2p1c;
//
// template <>
// InputParameters validParams<MoskitoTimeMass_2p1c>();

class MoskitoTimeMass_2p1c : public TimeKernel
{
public:
  static InputParameters validParams();
  MoskitoTimeMass_2p1c(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // The required values for pressure coupling
  const VariableValue & _p_dot;
  const VariableValue & _dp_dot;
  const unsigned int _p_var_number;

  // The required values for enthalpy coupling
  const VariableValue & _h_dot;
  const VariableValue & _dh_dot;
  const unsigned int _h_var_number;

  const VariableValue & _c_dot;
  const VariableValue & _dc_dot;
  const unsigned int _c_var_number;

  // Density derivatives wrt pressure and enthalpy
  const MaterialProperty<Real> & _drho_dp;
  const MaterialProperty<Real> & _drho_dp2;
  const MaterialProperty<Real> & _drho_dh;
  const MaterialProperty<Real> & _drho_dh2;
  const MaterialProperty<Real> & _drho_dph;
  const MaterialProperty<Real> & _drho_dc;
  const MaterialProperty<Real> & _drho_cdot;
};
