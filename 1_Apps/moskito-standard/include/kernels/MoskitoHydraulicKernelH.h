/**************************************************************************/
/*  MOSKITO - THMC sImulator for GEoscience Research                        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                        */
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

#ifndef MOSKITOHYDRAULICKERNELH_H
#define MOSKITOHYDRAULICKERNELH_H

#include "Kernel.h"
#include "RankTwoTensor.h"



class MoskitoHydraulicKernelH : public Kernel
{
public:
  static InputParameters validParams();
  MoskitoHydraulicKernelH(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const MaterialProperty<Real> & _scale_factor;
  const MaterialProperty<RankTwoTensor> & _k_vis;
  const MaterialProperty<Real> & _rho_f;
  const MaterialProperty<Real> & _drho_dp_f;
  const MaterialProperty<Real> & _drho_dT_f;
  const MaterialProperty<Real> & _mu_f;
  const MaterialProperty<Real> & _dmu_dp_f;
  const MaterialProperty<Real> & _dmu_dT_f;
  const MaterialProperty<RealVectorValue> & _g;
  unsigned int _temperature_var;
};

#endif // MOSKITOHYDRAULICKERNELH_H
