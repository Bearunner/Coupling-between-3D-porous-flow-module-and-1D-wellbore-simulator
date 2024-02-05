/**************************************************************************/
/*  Moskito - THMC sImulator for GEoscience Research                        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of Moskito App                                        */
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

#ifndef MOSKITOSOLUTETIMEKERNELS2_H
#define MOSKITOSOLUTETIMEKERNELS2_H

#include "TimeDerivative.h"



class MoskitoSoluteTimeKernelS2 : public TimeDerivative
{
public:
  static InputParameters validParams();
  MoskitoSoluteTimeKernelS2(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const VariableValue & _T_old;
  Real _T_threshold;
  const VariableValue & _u_old;
  const VariableValue & _original_tracer;
  // imported props from materials
  const MaterialProperty<Real> & _scale_factor;
  const MaterialProperty<Real> & _TimeKernelS;
  // const MaterialProperty<Real> & _new_C6_old;
  const MaterialProperty<RealVectorValue> & _SUPG_p;
  const MaterialProperty<bool> & _SUPG_ind;
};

#endif // MOSKITOSOLUTETIMEKERNELS2_H