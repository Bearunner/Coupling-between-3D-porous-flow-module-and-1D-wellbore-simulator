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

// class MoskitoLatHeatExcFormation_1p;
//
// template <>
// InputParameters validParams<MoskitoLatHeatExcFormation_1p>();

class MoskitoLatHeatExcFormation_1p : public Kernel
{
public:
  static InputParameters validParams();
  MoskitoLatHeatExcFormation_1p(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // Radius tubing outer
  const MaterialProperty<Real> & _rto;
  // Thermal wellbore resistivity
  const MaterialProperty<Real> & _Uto;
  // Formation temperature at the formation-casing interface
  const VariableValue & _Twf;
  unsigned _Twf_var_number;
  // Diameter filled with liquid = _rti
  const MaterialProperty<Real> & _area;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
};
