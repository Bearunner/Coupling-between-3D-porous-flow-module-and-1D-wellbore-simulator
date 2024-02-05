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

#include "MoskitoAdvectionBC.h"

registerMooseObject("MoskitoApp", MoskitoAdvectionBC);


InputParameters
MoskitoAdvectionBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Boundary conditions for outflow/outflow of advected quantities:"
                           "\n phi * velocity * normal, where phi is the advected quantitiy");
  params.addRequiredCoupledVar("velocity_vector",
                             "The components of the velocity vector up to problem dimension");
  return params;
}

MoskitoAdvectionBC::MoskitoAdvectionBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
   _dim(_mesh.dimension()),
   _coupled_components(coupledComponents("velocity_vector")),
   _velocity(coupledValues("velocity_vector"))
{
if (_dim > _coupled_components)
  paramError(
      "velocity_vector",
      "Number of components of velocity_vector must be at least equal to the mesh dimension");
if (_coupled_components > 3)
  paramError("velocity_vector",
             "You cannot supply more than 3 components for the velocity vector");
}

Real
MoskitoAdvectionBC::computeQpResidual()
{
  RealVectorValue vel;
  for (unsigned int j = 0; j < _coupled_components; ++j)
    vel(j) = (*_velocity[j])[_qp];
  if (vel * _normals[_qp] > 0)
    return _test[_i][_qp] * _u[_qp] * vel * _normals[_qp];
  return 0;
}

Real
MoskitoAdvectionBC::computeQpJacobian()
{
  RealVectorValue vel;
  for (unsigned int j = 0; j < _coupled_components; ++j)
    vel(j) = (*_velocity[j])[_qp];
  if (vel * _normals[_qp] > 0)
    return _test[_i][_qp] * _phi[_j][_qp] * vel * _normals[_qp];
  return 0;
}
