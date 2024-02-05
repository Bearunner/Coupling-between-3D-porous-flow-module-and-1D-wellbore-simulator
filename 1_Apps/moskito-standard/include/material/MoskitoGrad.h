/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*  Co-developed by Sebastian Held                                        */
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

#include "Material.h"
#include "MoskitoEOS2P.h"
#include "MoskitoDriftFlux.h"
#include "MoskitoViscosity2P.h"


// class MoskitoGrad;
//
// template <>
// InputParameters validParams<MoskitoGrad>();

class MoskitoGrad : public Material
{
public:
  static InputParameters validParams();
  MoskitoGrad(const InputParameters & parameters);
  virtual void computeQpProperties() override;


protected:

  MaterialProperty<Real> & _drho_m_dc;
  MaterialProperty<RealVectorValue> & _drho_m_dcvect;
  MaterialProperty<Real> & _dgamma_dc;
  MaterialProperty<RealVectorValue> & _dgamma_dcvect;
  MaterialProperty<Real> & _dkappa_dc;
  MaterialProperty<RealVectorValue> & _dkappa_dcvect;
  MaterialProperty<Real> & _domega_dc;
  MaterialProperty<RealVectorValue> & _domega_dcvect;
  MaterialProperty<Real> & _T_c;
  MaterialProperty<RealVectorValue> & _T_cvect;
  MaterialProperty<Real> & _drho_cdot;
  MaterialProperty<Real> & _dgamma_dcdot;

  const MaterialProperty<Real> & _drho_dc1;
  const MaterialProperty<Real> & _drho_dc2;
  const MaterialProperty<Real> & _drho_dc3;
  const MaterialProperty<Real> & _drho_dc4;
  const MaterialProperty<Real> & _drho_dc5;
  const MaterialProperty<Real> & _drho_dc6;
  const MaterialProperty<Real> & _drho_dc7;
  const MaterialProperty<Real> & _drho_dc8;
  const MaterialProperty<Real> & _dgamma_dc1;
  const MaterialProperty<Real> & _dgamma_dc2;
  const MaterialProperty<Real> & _dgamma_dc3;
  const MaterialProperty<Real> & _dgamma_dc4;
  const MaterialProperty<Real> & _dgamma_dc5;
  const MaterialProperty<Real> & _dgamma_dc6;
  const MaterialProperty<Real> & _dgamma_dc7;
  const MaterialProperty<Real> & _dgamma_dc8;
  const MaterialProperty<Real> & _dkappa_dc1;
  const MaterialProperty<Real> & _dkappa_dc2;
  const MaterialProperty<Real> & _dkappa_dc3;
  const MaterialProperty<Real> & _dkappa_dc4;
  const MaterialProperty<Real> & _dkappa_dc5;
  const MaterialProperty<Real> & _dkappa_dc6;
  const MaterialProperty<Real> & _dkappa_dc7;
  const MaterialProperty<Real> & _dkappa_dc8;
  const MaterialProperty<Real> & _domega_dc1;
  const MaterialProperty<Real> & _domega_dc2;
  const MaterialProperty<Real> & _domega_dc3;
  const MaterialProperty<Real> & _domega_dc4;
  const MaterialProperty<Real> & _domega_dc5;
  const MaterialProperty<Real> & _domega_dc6;
  const MaterialProperty<Real> & _domega_dc7;
  const MaterialProperty<Real> & _domega_dc8;
  const MaterialProperty<Real> & _T_c1;
  const MaterialProperty<Real> & _T_c2;
  const MaterialProperty<Real> & _T_c3;
  const MaterialProperty<Real> & _T_c4;
  const MaterialProperty<Real> & _T_c5;
  const MaterialProperty<Real> & _T_c6;
  const MaterialProperty<Real> & _T_c7;
  const MaterialProperty<Real> & _T_c8;

  const VariableValue & _c;
  const VariableValue & _c2;
  const VariableValue & _c3;
  const VariableValue & _c4;
  const VariableValue & _c5;
  const VariableValue & _c6;
  const VariableValue & _c7;
  const VariableValue & _c8;

  const VariableValue & _c_dot;
  const VariableValue & _c2_dot;
  const VariableValue & _c3_dot;
  const VariableValue & _c4_dot;
  const VariableValue & _c5_dot;
  const VariableValue & _c6_dot;
  const VariableValue & _c7_dot;
  const VariableValue & _c8_dot;


  const VariableGradient & _grad_c;
  const VariableGradient & _grad_c2;
  const VariableGradient & _grad_c3;
  const VariableGradient & _grad_c4;
  const VariableGradient & _grad_c5;
  const VariableGradient & _grad_c6;
  const VariableGradient & _grad_c7;
  const VariableGradient & _grad_c8;
};
