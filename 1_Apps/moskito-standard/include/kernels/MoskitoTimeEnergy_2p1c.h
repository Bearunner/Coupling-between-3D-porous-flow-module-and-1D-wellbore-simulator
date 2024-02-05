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

// class MoskitoTimeEnergy_2p1c;
//
// template <>
// InputParameters validParams<MoskitoTimeEnergy_2p1c>();

class MoskitoTimeEnergy_2p1c : public TimeKernel
{
public:
  static InputParameters validParams();
  MoskitoTimeEnergy_2p1c(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // required values for pressure and flowrate coupling
  const VariableValue & _m;
  const VariableValue & _p_dot;
  const VariableValue & _m_dot;
  const VariableValue & _c_dot;
  const VariableValue & _dp_dot;
  const VariableValue & _dm_dot;
  const VariableValue & _dc_dot;
  const unsigned int _p_var_number;
  const unsigned int _m_var_number;
  const unsigned int _c_var_number;

  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The density
  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _dgamma_dcdot;
  const MaterialProperty<Real> & _drho_cdot;
  // The first derivative of density wrt pressure
  const MaterialProperty<Real> & _drho_dp;
  // The first derivative of density wrt enthalpy
  const MaterialProperty<Real> & _drho_dh;
  const MaterialProperty<Real> & _drho_dc;
  // The second derivative of density wrt pressure
  const MaterialProperty<Real> & _drho_dp2;
  // The second derivative of density wrt enthalpy
  const MaterialProperty<Real> & _drho_dh2;
  // The second derivative of density wrt enthalpy and pressure
  const MaterialProperty<Real> & _drho_dph;

  // The gamma first derivatives
  const MaterialProperty<Real> & _dgamma_dp;
  // The gamma first derivatives
  const MaterialProperty<Real> & _dgamma_dh;
  // The gamma first derivatives
  const MaterialProperty<Real> & _dgamma_dm;
  const MaterialProperty<Real> & _dgamma_dc;
};
