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

#ifndef MOSKITOSOLUTESOURCEKERNELS_H
#define MOSKITOSOLUTESOURCEKERNELS_H

#include "Kernel.h"


class Function;



class MoskitoSoluteSourceKernelS : public Kernel
{
public:
  static InputParameters validParams();
  MoskitoSoluteSourceKernelS(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  const Real & _scale;
  const MaterialProperty<Real> & _n;
  const Function & _function;
  const MaterialProperty<Real> & _K;

  // imported props from materials
  const MaterialProperty<Real> & _scale_factor;
  const MaterialProperty<RealVectorValue> & _SUPG_p;
  const MaterialProperty<bool> & _SUPG_ind;
};

#endif  //MOSKITOSOLUTESOURCEKERNELS_H
