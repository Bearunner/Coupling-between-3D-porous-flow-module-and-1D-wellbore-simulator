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

#ifndef MOSKITOSOLUTEMASSCHANGERATE_H
#define MOSKITOSOLUTEMASSCHANGERATE_H

#include "AuxKernel.h"



class MoskitoSoluteMassChangeRate : public AuxKernel
{
public:
  static InputParameters validParams();
  MoskitoSoluteMassChangeRate(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:

  // const VariableValue & _T_old;
  const VariableValue & _C_dot;
  // Real _T_threshold;
  // // imported props from materials
  // const VariableValue & _u_old;
  // const MaterialProperty<Real> & _scale_factor;
  // const MaterialProperty<Real> & _TimeKernelS;
  // const MaterialProperty<Real> & _new_C5_old;
  // const MaterialProperty<RealVectorValue> & _SUPG_p;
  // const MaterialProperty<bool> & _SUPG_ind;
};

#endif // MOSKITOSOLUTEMASSCHANGERATE_H
