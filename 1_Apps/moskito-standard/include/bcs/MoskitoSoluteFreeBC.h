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

#ifndef MOSKITOSOLUTEFREEBC_H
#define MOSKITOSOLUTEFREEBC_H

#include "IntegratedBC.h"
#include "RankTwoTensor.h"

/* The residual is simply -test*k*grad_u*normal which is the term
 * you get from integration by parts.
 *
 * See also: Griffiths, David F. "The 'no boundary condition' outflow
 * boundary condition.", International Journal for Numerical Methods
 * in Fluids, vol. 24, no. 4, 1997, pp. 393-411.
 */
class MoskitoSoluteFreeBC : public IntegratedBC
{
public:
 static InputParameters validParams();
 MoskitoSoluteFreeBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


  // const unsigned int _dim;
  // const unsigned int _coupled_components;
  // const std::vector<const VariableValue *> _velocity;

  const MaterialProperty<RealVectorValue> & _SUPG_p;
  const MaterialProperty<bool> & _SUPG_ind;
  const MaterialProperty<bool> & _av_ind;
  const MaterialProperty<RealVectorValue> & _av;
  const MaterialProperty<RealVectorValue> * _dav_dp_phi;
  const MaterialProperty<RankTwoTensor> * _dav_dp_gradphi;
  unsigned int _pressure_var;
};

#endif // MOSKITOSOLUTEFREEBC_H
