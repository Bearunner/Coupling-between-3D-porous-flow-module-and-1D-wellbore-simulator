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

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "MoskitoEOS2P.h"


/*class Moskito_2pmc;*/

/*template <>*/
/*InputParameters validParams<Moskito_2pmc>();*/



class Moskito_2pmc : public DerivativeMaterialInterface<Material>//, public NewtonIteration
{
public:
  static InputParameters validParams();
  Moskito_2pmc(const InputParameters & parameters);
  virtual void computeQpProperties() override;


protected:
  // Userobject to equation of state
  const MoskitoEOS2P & eos_2p;

  // Density of gas
  MaterialProperty<Real> & _rho_gas;
  // Density of liquid
  MaterialProperty<Real> & _rho_liq;
  // Density of mixture
  MaterialProperty<Real> & _rho_mix;
  // Profile-adjusted density of mixture


  // The first derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_mix_dp;
  // The second derivative of mixture density wrt pressure
  MaterialProperty<Real> & _drho_mix_dp2;


  // The coupled enthalpy
  const VariableValue & _h;
  const VariableValue & _m;
  const VariableValue & _p;
  VariableName _p_name;

  //const VariableValue & _c;
  //VariableName _c_name;

  // The gradient of the coupled variables
  const VariableGradient & _grad_m;
  const VariableGradient & _grad_h;
  const VariableGradient & _grad_p;

  //MaterialProperty<Real> & F;
  MaterialProperty<Real> & _d_rho_mixd_p;
  MaterialProperty<Real> & _new_gradient;
};
