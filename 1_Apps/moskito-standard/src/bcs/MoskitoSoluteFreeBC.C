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

#include "MoskitoSoluteFreeBC.h"

registerMooseObject("MoskitoApp", MoskitoSoluteFreeBC);


InputParameters
MoskitoSoluteFreeBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Boundary conditions for outflow/outflow of advected quantities:"
                           "\n phi * velocity * normal, where phi is the advected quantitiy");
  params.addCoupledVar("pressure", "Pore pressure nonlinear variable");

  // params.addRequiredCoupledVar("velocity_vector",
  //                            "The components of the velocity vector up to problem dimension");
  return params;
}

MoskitoSoluteFreeBC::MoskitoSoluteFreeBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
   // _dim(_mesh.dimension()),
   // _coupled_components(coupledComponents("velocity_vector")),
   // _velocity(coupledValues("velocity_vector"))
   _SUPG_p(getMaterialProperty<RealVectorValue>("solute_petrov_supg_p_function")),
   _SUPG_ind(getMaterialProperty<bool>("solute_supg_indicator")),
   _av_ind(getMaterialProperty<bool>("solute_av_dv_indicator")),
   _av(getMaterialProperty<RealVectorValue>("solute_advection_velocity")),
   _pressure_var(coupled("pressure"))
{
// if (_dim > _coupled_components)
//   paramError(
//       "velocity_vector",
//       "Number of components of velocity_vector must be at least equal to the mesh dimension");
// if (_coupled_components > 3)
//   paramError("velocity_vector",
//              "You cannot supply more than 3 components for the velocity vector");
if (parameters.isParamSetByUser("pressure"))
{
  // _dav_dT = &getMaterialProperty<RealVectorValue>("d_well_velocity_dT");
  _dav_dp_phi = &getMaterialProperty<RealVectorValue>("d_well_velocity_dp_phi");
  _dav_dp_gradphi = &getMaterialProperty<RankTwoTensor>("d_well_velocity_dp_gradphi");
}
else
{
  // _dav_dT = NULL;
  _dav_dp_phi = NULL;
  _dav_dp_gradphi = NULL;
}
}

Real
MoskitoSoluteFreeBC::computeQpResidual()
{
  // RealVectorValue vel;
  // for (unsigned int j = 0; j < _coupled_components; ++j)
  //   vel(j) = (*_velocity[j])[_qp];
  // if (vel * _normals[_qp] > 0)
  //   return _test[_i][_qp] * _u[_qp] * vel * _normals[_qp];
  // return 0;
  Real res = 0.0;
  res = _test[_i][_qp] * _av[_qp] * _u[_qp] * _normals[_qp];

  if (_SUPG_ind[_qp])
    // res +=  -_SUPG_p[_qp] * _grad_test[_i][_qp] * _av[_qp] * _grad_u[_qp] * _normals[_qp]; //not 100% sure about this part!!!
    res +=  -_SUPG_p[_qp] * _test[_i][_qp] * _av[_qp] * _grad_u[_qp] * _normals[_qp]; //not 100% sure about this part!!!

  else
    res += 0.0;

  return res;
}

Real
MoskitoSoluteFreeBC::computeQpJacobian()
{
  // RealVectorValue vel;
  // for (unsigned int j = 0; j < _coupled_components; ++j)
  //   vel(j) = (*_velocity[j])[_qp];
  // if (vel * _normals[_qp] > 0)
  //   return _test[_i][_qp] * _phi[_j][_qp] * vel * _normals[_qp];
  // return 0;
  Real j = 0.0;

  j = _test[_i][_qp] * _av[_qp] * _phi[_j][_qp] * _normals[_qp];

  if (_SUPG_ind[_qp])
    // j += -_SUPG_p[_qp] * _grad_test[_i][_qp] * _av[_qp] * _grad_phi[_j][_qp] * _normals[_qp];//not 100% sure about this part!!!
    j += -_SUPG_p[_qp] * _test[_i][_qp] * _av[_qp] * _grad_phi[_j][_qp] * _normals[_qp];//not 100% sure about this part!!!

  else
    j += 0.0;

  return j;

}


Real
MoskitoSoluteFreeBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Real test = 0.0 , j = 0.0;
  //
  // if (jvar == _pressure_var)
  // {
  //   if (_SUPG_ind[_qp])
  //     test = _test[_i][_qp] + _SUPG_p[_qp] * _grad_test[_i][_qp];
  //   else
  //     test = _test[_i][_qp];
  //
  //   j  = (*_dav_dp_phi)[_qp] * _phi[_j][_qp] * _grad_u[_qp];
  //   j += (*_dav_dp_gradphi)[_qp] * _grad_phi[_j][_qp] * _grad_u[_qp];
  //   j *= _scale_factor[_qp] * test;
    // j  = (*_dav_dp_phi)[_qp] * _phi[_j][_qp] * _u[_qp] * _scale_factor[_qp] * _grad_test[_i][_qp];
    // j += (*_dav_dp_gradphi)[_qp] * _grad_phi[_j][_qp] * _u[_qp] * _scale_factor[_qp] * _grad_test[_i][_qp];
    Real j = 0.0;

    if (jvar == _pressure_var)
    {
      j = _test[_i][_qp] * ((*_dav_dp_phi)[_qp] * _phi[_j][_qp] + (*_dav_dp_gradphi)[_qp] * _grad_phi[_j][_qp]) * _u[_qp] * _normals[_qp];

      if (_SUPG_ind[_qp])
        // j += -_SUPG_p[_qp] * _grad_test[_i][_qp] * ((*_dav_dp_phi)[_qp] * _phi[_j][_qp] + (*_dav_dp_gradphi)[_qp] * _grad_phi[_j][_qp]) * _grad_u[_qp] * _normals[_qp];
        j += -_SUPG_p[_qp] * _test[_i][_qp] * ((*_dav_dp_phi)[_qp] * _phi[_j][_qp] + (*_dav_dp_gradphi)[_qp] * _grad_phi[_j][_qp]) * _grad_u[_qp] * _normals[_qp];

      else
        j += 0.0;
  }

  return j;
}
