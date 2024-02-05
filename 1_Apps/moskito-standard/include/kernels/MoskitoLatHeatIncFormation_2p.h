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

#include "Kernel.h"

// class MoskitoLatHeatIncFormation_2p;
//
// template <>
// InputParameters validParams<MoskitoLatHeatIncFormation_2p>();

class MoskitoLatHeatIncFormation_2p : public Kernel
{
public:
  static InputParameters validParams();
  MoskitoLatHeatIncFormation_2p(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  unsigned _p_var_number;
  unsigned _c_var_number;

  // Thermal wellbore resistivity
  const MaterialProperty<Real> & _lambda;
  // formation Temperature
  const MaterialProperty<Real> & _Tform;
  const MaterialProperty<Real> & _area;
  const MaterialProperty<Real> & _T;
  const MaterialProperty<Real> & _T_h;
  const MaterialProperty<Real> & _T_p;
  const MaterialProperty<Real> & _T_c;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
};
