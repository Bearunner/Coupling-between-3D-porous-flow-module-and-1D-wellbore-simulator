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

#include "MoskitoFluidWellGeneral.h"
#include "MoskitoEOS1P.h"
#include "MoskitoViscosity1P.h"

// class MoskitoFluidWell_1p1c;
//
// template <>
// InputParameters validParams<MoskitoFluidWell_1p1c>();

class MoskitoFluidWell_1p1c : public MoskitoFluidWellGeneral
{
public:
  static InputParameters validParams();
  MoskitoFluidWell_1p1c(const InputParameters & parameters);
  virtual void computeQpProperties() override;

protected:
  // Userobject to equation of state
  const MoskitoEOS1P & eos_uo;
  // Userobject to Viscosity Eq
  const MoskitoViscosity1P & viscosity_uo;
  // The convective heat transfer factor of fluid in coaxial configuration
  MaterialProperty<Real> & _hf;
  // The vescosity
  MaterialProperty<Real> & _vis;
  MaterialProperty<Real> & _dmu_dp;
  MaterialProperty<Real> & _dmu_dT;
  // The constant thermal conductivity of fluid
  MaterialProperty<Real> & _lambda;
  // The specific heat at constant pressure
  MaterialProperty<Real> & _cp;
  // The density
  MaterialProperty<Real> & _rho;
  // The first derivative of density wrt pressure
  MaterialProperty<Real> & _drho_dp;
  // The first derivative of density wrt temperature
  MaterialProperty<Real> & _drho_dT;
  // Enthalpy from P and T
  MaterialProperty<Real> & _h;

  // The coupled temperature
  const VariableValue & _T;
  // The coupled flow rate
  const VariableValue & _flow;

  MaterialProperty<RealVectorValue> & _wv;
  MaterialProperty<RealVectorValue> & _dwv_dT;
  // Derivative of Dracy velocity wrt pressure to multiply by _phi
  MaterialProperty<RealVectorValue> & _dwv_dp_phi;
  // Derivative of Dracy velocity wrt pressure to multiply by _grad_phi
  MaterialProperty<RankTwoTensor> & _dwv_dp_gradphi;

  // function for calculating convective heat transfer coeff
  Real Conv_coeff();
};
