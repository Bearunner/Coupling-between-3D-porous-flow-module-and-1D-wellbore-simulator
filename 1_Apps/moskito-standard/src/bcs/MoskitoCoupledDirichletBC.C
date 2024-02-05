/**************************************************************************/
/*  TIGER - THMC sImulator for GEoscience Research                        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of TIGER App                                        */
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

#include "MoskitoCoupledDirichletBC.h"

registerMooseObject("MoskitoApp", MoskitoCoupledDirichletBC);

// template <>
// InputParameters
// validParams<MoskitoCoupledDirichletBC>()
// {
//   InputParameters params = validParams<NodalBC>();

InputParameters
MoskitoCoupledDirichletBC::validParams()
{
  InputParameters params = NodalBC::validParams();

  params.addRequiredCoupledVar("coupled_var", "Value on the Boundary");
  return params;
}

MoskitoCoupledDirichletBC::MoskitoCoupledDirichletBC(const InputParameters & parameters)
  : NodalBC(parameters),
  _coupled_var(coupledValue("coupled_var")),
  _coupled_var_number(coupled("coupled_var"))
{
}

Real
MoskitoCoupledDirichletBC::computeQpResidual()
{
  return _u[_qp] - _coupled_var[_qp];
}

Real
MoskitoCoupledDirichletBC::computeQpJacobian()
{
  return 1.0;
}

Real
MoskitoCoupledDirichletBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real j = 0.0;

  if (jvar == _coupled_var_number)
  {
    j = -1.0;
  }

  return j;
}
