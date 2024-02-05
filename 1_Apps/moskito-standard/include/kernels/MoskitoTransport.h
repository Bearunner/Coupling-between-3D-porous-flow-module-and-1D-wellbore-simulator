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

// class MoskitoTransport;
//
// template <>
// InputParameters validParams<MoskitoTransport>();

class MoskitoTransport : public Kernel
{
public:
  static InputParameters validParams();
  MoskitoTransport(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // The gradient of the coupled flow_rate
  const VariableGradient & _grad_m;
  unsigned _m_var_number;

  // The gradient of the coupled temperature
  const VariableGradient & _grad_p;
  unsigned _p_var_number;

  const VariableGradient & _grad_h;

  // Variable numberings
  unsigned _h_var_number;

  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
  // The sign of well flow direction
  const MaterialProperty<Real> & _well_sign;

  const MaterialProperty<Real> & _aux1;

  const MaterialProperty<Real> & _aux2;

  const MaterialProperty<Real> & _aux1_p;

  const MaterialProperty<Real> & _aux2_p;

  const MaterialProperty<Real> & _u_l_p;

  const MaterialProperty<Real> & _u_g_p;

  const MaterialProperty<Real> & _aux1_h;

  const MaterialProperty<Real> & _aux2_h;

  const MaterialProperty<Real> & _u_l_h;

  const MaterialProperty<Real> & _u_g_h;

  const MaterialProperty<Real> & _aux1_c;

  const MaterialProperty<Real> & _aux2_c;

  const MaterialProperty<Real> & _u_l_c;

  const MaterialProperty<Real> & _u_g_c;

  const MaterialProperty<Real> & _u_l_m;

  const MaterialProperty<Real> & _u_g_m;

  const MaterialProperty<Real> & _u_g;

  const MaterialProperty<Real> & _u_l;

/*
  // The area of pipe
  const MaterialProperty<Real> & _area;
  // The unit vector of well direction
  const MaterialProperty<RealVectorValue> & _well_dir;
  // The sign of well flow direction
  const MaterialProperty<Real> & _well_sign;
  // The density
  const MaterialProperty<Real> & _rho;
  // The first derivative of density wrt pressure
  const MaterialProperty<Real> & _drho_dp;
  // The first derivative of density wrt temperature
  const MaterialProperty<Real> & _drho_dT;
  */
};
