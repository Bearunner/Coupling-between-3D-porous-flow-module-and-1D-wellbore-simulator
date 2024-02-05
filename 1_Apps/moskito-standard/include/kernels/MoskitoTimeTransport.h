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

// class MoskitoTimeTransport;
//
// template <>
// InputParameters validParams<MoskitoTimeTransport>();

class MoskitoTimeTransport : public TimeKernel
{
public:
  static InputParameters validParams();
  MoskitoTimeTransport(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // required values for temperature coupling
  const VariableValue & _h_dot;
  const VariableValue & _dh_dot;
  const unsigned int _h_var_number;
  const VariableValue & _p_dot;
  const VariableValue & _dp_dot;
  const unsigned int _p_var_number;

  // The first derivative of density wrt pressure
  //const MaterialProperty<Real> & _drho_dp;
  // The first derivative of density wrt temperature
  //const MaterialProperty<Real> & _drho_dT;

  const MaterialProperty<Real> & _aux1;

  const MaterialProperty<Real> & _aux2;

  const MaterialProperty<Real> & _aux1_p;

  const MaterialProperty<Real> & _aux2_p;

  const MaterialProperty<Real> & _aux1_h;

  const MaterialProperty<Real> & _aux2_h;

  const MaterialProperty<Real> & _aux1_c;

  const MaterialProperty<Real> & _aux2_c;


};
