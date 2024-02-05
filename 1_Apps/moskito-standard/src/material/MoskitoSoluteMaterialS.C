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

#include "MoskitoSoluteMaterialS.h"
#include "MooseMesh.h"
#include "libmesh/quadrature.h"
#include "Material.h"
#include "RankTwoTensor.h"
#include <cfloat>
#include "Function.h"

registerMooseObject("MoskitoApp", MoskitoSoluteMaterialS);

InputParameters
MoskitoSoluteMaterialS::validParams()
{
  InputParameters params = Material::validParams();

  MooseEnum Advection
        ("pure_diffusion well_velocity user_velocity well_user_velocities",
        "well_velocity");
  params.addParam<MooseEnum>("advection_type", Advection,
        "Type of the velocity to simulate advection [pure_diffusion "
        "well_velocity user_velocity well_user_velocities]");
    params.addParam<bool>("output_Pe_Cr_numbers", false ,
        "calcuate Peclet and Courant numbers");
  params.addParam<bool>("has_supg", false ,
        "Is Streameline Upwinding / Petrov Galerkin (SU/PG) activated?");
  params.addParam<Real>("supg_coeficient_scale", 1.0 ,
        "The user defined factor to scale SU/PG coefficent (tau)");
  params.addParam<FunctionName>("user_velocity", 0.0,
        "a vector function to define the velocity field");
  params.addParam<UserObjectName>("supg_uo", "",
        "The name of the userobject for SU/PG");
  params.addRequiredParam<Real>("diffusion", "Molecular diffusion of component in water (m^2/s) - something like 2e-9");
  params.addParam<Real>("dispersion_longitudinal", 0, "Longitudinal dispersivity (m)");
  params.addParam<Real>("dispersion_transverse", 0, "Transverse dispersivity (m)");
  params.addParam<Real>("formation_factor", 1, "The formation factor, 0.1 for clays, 0.7 for sand, depending on the tortuosity of the porous medium");
  params.addCoupledVar("temperature", 273.15, "temperature nonlinear variable (K)");
  params.addCoupledVar("C1", 0, "C1");
  params.addCoupledVar("C2", 0, "C2");
  params.addCoupledVar("C3", 0, "C3");
  params.addCoupledVar("C4", 0, "C4");
  params.addCoupledVar("C5", 0, "C5");
  params.addCoupledVar("C6", 0, "C6");
  params.addParam<Real>("C1_temperature_threshold", 0, "temperature_threshold for C1 temperature-response (K)");
  params.addParam<Real>("C3_temperature_threshold", 0, "temperature_threshold for C3 temperature-response (K)");
  params.addParam<Real>("C5_temperature_threshold", 0, "temperature_threshold for C5 temperature-response (K)");
  // params.addParam<FunctionName>("initial_temperature", 273.15, "initial temperature)");
  // params.addParam<FunctionName>("initial_C1", "initial C1)");
  // params.addParam<FunctionName>("initial_C3", "initial C3)");
  // params.addParam<FunctionName>("initial_C5", "initial C5)");
  params.addParam<Real>("porosity", 1.0, "porosity (temporal and spatial function)");
  params.addParam<FunctionName>("scale_factor", 1.0, "The scale factor for non-3D "
      "elements ( particularlly lower dimensional elements): if mesh is 3D"
      ", aperture for 2D elements (fractures) and diameter for 1D elements"
      " (wells) should be used; if mesh is 2D, height for 2D elements (2D"
      " matrix) and aperture times height for 1D elements (fractures) "
      "should be used; and if mesh is 1D, area for 1D elements (pipes or "
      "wells) should be used); ParsedFunctions can be used as well.");

  params.addClassDescription("Solute material for solute kernels");
  return params;
}

MoskitoSoluteMaterialS::MoskitoSoluteMaterialS(const InputParameters & parameters)
  : Material(parameters),
    _at(getParam<MooseEnum>("advection_type")),
    _T(coupledValue("temperature")),
    _C1(coupledValue("C1")),
    _C2(coupledValue("C2")),
    _C3(coupledValue("C3")),
    _C4(coupledValue("C4")),
    _C5(coupledValue("C5")),
    _C6(coupledValue("C6")),
    _C1_T_threshold(getParam<Real>("C1_temperature_threshold")),
    _C3_T_threshold(getParam<Real>("C3_temperature_threshold")),
    _C5_T_threshold(getParam<Real>("C5_temperature_threshold")),
    _new_C1(declareProperty<Real>("new_C1")),
    _new_C2(declareProperty<Real>("new_C2")),
    _new_C3(declareProperty<Real>("new_C3")),
    _new_C4(declareProperty<Real>("new_C4")),
    _new_C5(declareProperty<Real>("new_C5")),
    _new_C6(declareProperty<Real>("new_C6")),
    _C1_indicator(declareProperty<Real>("C1_indicator")),
    _C3_indicator(declareProperty<Real>("C3_indicator")),
    _C5_indicator(declareProperty<Real>("C5_indicator")),
    _has_PeCr(getParam<bool>("output_Pe_Cr_numbers")),
    _has_supg(getParam<bool>("has_supg")),
    _supg_scale(getParam<Real>("supg_coeficient_scale")),
    _TimeKernelS(declareProperty<Real>("TimeKernel_S")),
    _SUPG_ind(declareProperty<bool>("solute_supg_indicator")),
    _av_ind(declareProperty<bool>("solute_av_dv_indicator")),
    _av(declareProperty<RealVectorValue>("solute_advection_velocity")),
    _SUPG_p(declareProperty<RealVectorValue>("solute_petrov_supg_p_function")),
    _rot_mat(declareProperty<RankTwoTensor>("lowerD_rotation_matrix")),
    _diffusion_molecular(getParam<Real>("diffusion")),
    _disp_l(getParam<Real>("dispersion_longitudinal")),
    _disp_t(getParam<Real>("dispersion_transverse")),
    _formation_factor(getParam<Real>("formation_factor")),
    _n0(getParam<Real>("porosity")),
    _n(declareProperty<Real>("porosity")),
    _dispersion_tensor(declareProperty<RankTwoTensor>("dispersion_tensor")),
    _diffusion_factor(declareProperty<Real>("diffusion_factor")),
    _diffdisp(declareProperty<RankTwoTensor>("diffusion_dispersion")),
    _Fo(declareProperty<Real>("neumann_number")),
    _PeDisp(declareProperty<Real>("peclet_number_dispersive")),
    _scale_factor(declareProperty<Real>("scale_factor")),
    _scale_factor0(getFunction("scale_factor"))
    // _K(declareProperty<Real>("color_indicator")),
    // _initial_temperature(getFunction("initial_temperature")),
    // _T_ini(declareProperty<Real>("initial_temperature")),
    // _initial_C1(getFunction("initial_C1")),
    // _C1_ini(declareProperty<Real>("initial_C1")),
    // _C1_indicator_ini(declareProperty<Real>("C1_indicator_ini")),
    // _initial_C3(getFunction("initial_C3")),
    // _C3_ini(declareProperty<Real>("initial_C3")),
    // _C3_indicator_ini(declareProperty<Real>("C3_indicator_ini")),
    // _initial_C5(getFunction("initial_C5")),
    // _C5_ini(declareProperty<Real>("initial_C5")),
    // _C5_indicator_ini(declareProperty<Real>("C5_indicator_ini")),
{
  _Pe = (_has_PeCr || _has_supg) ?
              &declareProperty<Real>("solute_peclet_number") : NULL;
  _Cr = (_has_PeCr || _has_supg) ?
              &declareProperty<Real>("solute_courant_number") : NULL;
  _vel_func = (_at == AT::user_velocity || _at == AT::well_user_velocities) ?
              &getFunction("user_velocity") : NULL;
  _supg_uo = (parameters.isParamSetByUser("supg_uo")) ?
              &getUserObject<MoskitoSUPG>("supg_uo") : NULL;
  _wv = (_at == AT::well_velocity || _at == AT::well_user_velocities) ?
              &getMaterialProperty<RealVectorValue>("well_velocity_vector") : NULL;
}

void
MoskitoSoluteMaterialS::initQpStatefulProperties()
{
  // _C1_ini[_qp] = _initial_C1.value(_t, _q_point[_qp]);
  // _T_ini[_qp] = _initial_temperature.value(_t, _q_point[_qp]);
  //
  // if (_T_ini[_qp]  < 358.0) {
  //   _C1_indicator_ini[_qp] = 1.0;
  // }
  // else if (_T_ini[_qp] >=358.0 && _T_ini[_qp] <=358.3){
  //   _C1_indicator_ini[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T_ini[_qp] - 358.15)));
  // }
  // else {
  //   _C1_indicator_ini[_qp] = 0.0;
  // }
  //
  // _new_C1[_qp] = _C1_indicator_ini[_qp] * _C1_ini[_qp];
  // _new_C2[_qp] = 0.0;
  //
  //
  // _C3_ini[_qp] = _initial_C3.value(_t, _q_point[_qp]);
  // // _T_ini[_qp] = _initial_temperature.value(_t, _q_point[_qp]);
  //
  // if (_T_ini[_qp]  < 363.0) {
  //   _C3_indicator_ini[_qp] = 1.0;
  // }
  // else if (_T_ini[_qp] >=363.0 && _T_ini[_qp] <=363.3){
  //   _C3_indicator_ini[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T_ini[_qp] - 363.15)));
  // }
  // else {
  //   _C3_indicator_ini[_qp] = 0.0;
  // }
  //
  // _new_C3[_qp] = _C3_indicator_ini[_qp] * _C3_ini[_qp];
  // _new_C4[_qp] = 0.0;
  computeQpProperties();
}

void
MoskitoSoluteMaterialS::computeQpProperties()
{
  _n[_qp] = _n0;
    // Chemical kernel for calculating the time derivative, n0 is porosity
 _TimeKernelS[_qp] = _n[_qp];

 if (_T[_qp]  < (_C1_T_threshold - 0.15)) {
   _C1_indicator[_qp] = 1.0;
 }
 else if (_T[_qp] >=(_C1_T_threshold - 0.15) && _T[_qp] <=(_C1_T_threshold + 0.15)){
   _C1_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C1_T_threshold)));
 }
 else {
   _C1_indicator[_qp] = 0.0;
 }

 if (_T[_qp]  < (_C3_T_threshold - 0.15)) {
   _C3_indicator[_qp] = 1.0;
 }
 else if (_T[_qp] >=(_C3_T_threshold - 0.15) && _T[_qp] <=(_C3_T_threshold + 0.15)){
   _C3_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C3_T_threshold)));
 }
 else {
   _C3_indicator[_qp] = 0.0;
 }

 if (_T[_qp]  < (_C5_T_threshold - 0.15)) {
   _C5_indicator[_qp] = 1.0;
 }
 else if (_T[_qp] >=(_C5_T_threshold - 0.15) && _T[_qp] <=(_C5_T_threshold + 0.15)){
   _C5_indicator[_qp] = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T[_qp] - _C5_T_threshold)));
 }
 else {
   _C5_indicator[_qp] = 0.0;
 }

_new_C1[_qp] = _C1_indicator[_qp] * _C1[_qp];
_new_C2[_qp] = _C1[_qp] - _new_C1[_qp] + _C2[_qp];

_new_C3[_qp] = _C3_indicator[_qp] * _C3[_qp];
_new_C4[_qp] = _C3[_qp] - _new_C3[_qp] + _C4[_qp];

_new_C5[_qp] = _C5_indicator[_qp] * _C5[_qp];
_new_C6[_qp] = _C5[_qp] - _new_C5[_qp] + _C6[_qp];

	Real h_n = _current_elem->hmin();
	_Fo[_qp] = _diffusion_molecular * _dt / (h_n*h_n);

  _diffusion_factor[_qp] = _diffusion_molecular * _n[_qp] * _formation_factor;

  switch (_at)
  {
    case AT::pure_diffusion:
      _av[_qp].zero();
      _av_ind[_qp] = false;
    break;
    case AT::well_velocity:
      _av[_qp] = (*_wv)[_qp];
      _av_ind[_qp] = true;
      break;
    case AT::user_velocity:
      _av[_qp] = _vel_func->vectorValue(_t, _q_point[_qp]);
      _av_ind[_qp] = false;
      break;
    case AT::well_user_velocities:
      _av[_qp] = (*_wv)[_qp] + _vel_func->vectorValue(_t, _q_point[_qp]);
      _av_ind[_qp] = true;
      break;
  }

  RealVectorValue wellLocal = _rot_mat[_qp].transpose() * _av[_qp];

 _diffdisp[_qp] = DispersionTensorCalculator(wellLocal, _disp_l, _disp_t, _current_elem->dim(), _mesh.dimension(), _diffusion_factor[_qp]);


  if (_current_elem->dim() < _mesh.dimension())
    _diffdisp[_qp].rotate(_rot_mat[_qp]);

  Real lambda = _diffdisp[_qp].trace() / (_current_elem->dim() * _TimeKernelS[_qp]);

  if (_has_PeCr && !_has_supg)
    _supg_uo->PeCrNrsCalculator(lambda, _dt, _current_elem, _av[_qp], (*_Pe)[_qp], (*_Cr)[_qp]);

  if (_has_supg)
  {
    // should be multiplied by the gradient of the test function to build the Petrov Galerkin P function
    _supg_uo->SUPGCalculator(lambda, _dt, _current_elem, _av[_qp], _SUPG_p[_qp], (*_Pe)[_qp], (*_Cr)[_qp]);
    _SUPG_p[_qp] *= _supg_scale;

    if (_SUPG_p[_qp].norm() != 0.0)
      _SUPG_ind[_qp] = true;
    else
      _SUPG_ind[_qp] = false;
  }
  else
    _SUPG_ind[_qp] = false;

    _scale_factor[_qp] = Scaling();

    if (_current_elem->dim() < _mesh.dimension())
    _rot_mat[_qp] = lowerDRotationMatrix(_current_elem->dim());
    else
    _rot_mat[_qp] = RankTwoTensor::Identity();
}

RankTwoTensor
MoskitoSoluteMaterialS::DispersionTensorCalculator(const RealVectorValue & w_v, Real const & disp_l, Real const & disp_t, int dim, int dimMesh, Real diffusion_factor)
{

  if (w_v.norm() != 0)
  {

    Real d00 = (1/w_v.norm()) * ((disp_t * (w_v(1) * w_v(1) + w_v(2) * w_v(2)) + disp_l * w_v(0) * w_v(0)));
    d00 += diffusion_factor;

    Real d01 = (1/w_v.norm()) * (((disp_l - disp_t) * w_v(0) * w_v(1)));
    Real d02 = (1/w_v.norm()) * (((disp_l - disp_t) * w_v(0) * w_v(2)));

    Real d11 = 0;
    if (dimMesh >= dim && dim > 1)
    {
      d11 += (1/w_v.norm()) * ((disp_t * (w_v(0) * w_v(0) + w_v(2) * w_v(2)) + disp_l * w_v(1) * w_v(1)));
      d11 += diffusion_factor;
    }

    Real d12 = (1/w_v.norm()) * (((disp_l - disp_t) * w_v(1) * w_v(2)));

    Real d22 = 0;
    if (dim == dimMesh && dim > 2)
    {
      d22 += (1/w_v.norm()) * ((disp_t * (w_v(0) * w_v(0) + w_v(1) * w_v(1)) + disp_l * w_v(2) * w_v(2)));
      d22 += diffusion_factor;
    }

    _dispersion_ten = (RankTwoTensor(d00, d11, d22, d12, d02, d01));
  }
  else
  _dispersion_ten = (RankTwoTensor(diffusion_factor, diffusion_factor, diffusion_factor, 0., 0., 0.));

  return _dispersion_ten;
}

Real
MoskitoSoluteMaterialS::Scaling()
{
  Real scale_factor = 1.0;

  switch (_mesh.dimension())
  {
    case 1 ... 2:
      scale_factor = _scale_factor0.value(_t, _q_point[_qp]);
      break;
    case 3:
      if (_current_elem->dim() == 2)
        // fracture aperture
        scale_factor = _scale_factor0.value(_t, _q_point[_qp]);
      else if (_current_elem->dim() == 1)
       // radius of well
        scale_factor = PI * _scale_factor0.value(_t, _q_point[_qp]) * _scale_factor0.value(_t, _q_point[_qp]) / 4.0;
      break;
  }

  return scale_factor;
}

RankTwoTensor
MoskitoSoluteMaterialS::lowerDRotationMatrix(int dim)
{
  RealVectorValue xp, yp, zp;
  xp = _current_elem->point(1) - _current_elem->point(0);

  switch (dim)
  {
    case 1:
      for (unsigned int i = 0; i < 3; ++i)
        yp(i) = 0.0;
      if (std::fabs(xp(0)) > 0.0 && std::fabs(xp(1)) + std::fabs(xp(2)) < DBL_MIN)
        yp(2) = 1.0;
      else if (std::fabs(xp(1)) > 0.0 && std::fabs(xp(0)) + std::fabs(xp(2)) < DBL_MIN)
        yp(0) = 1.0;
      else if (std::fabs(xp(2)) > 0.0 && std::fabs(xp(0)) + std::fabs(xp(1)) < DBL_MIN)
        yp(1) = 1.0;
      else
      {
        for (unsigned int i = 0; i < 3; ++i)
          if (std::fabs(xp(i)) > 0.0)
          {
            yp(i) = -xp(i);
            break;
          }
      }
      zp = xp.cross(yp);
      yp = zp.cross(xp);
      break;

    case 2:
      yp = _current_elem->point(2) - _current_elem->point(1);
      zp = xp.cross(yp);
      if (!((std::fabs(zp(0)) + std::fabs(zp(1)))/zp.norm() < DBL_MIN)) //horizontal fracture check
        xp = RealVectorValue(0.,0.,1.).cross(zp);
      else
        xp = RealVectorValue(1.,0.,0.);
      yp = zp.cross(xp);
      break;
  }

  RankTwoTensor _rm;
  for (unsigned int i = 0; i < 3; ++i)
  {
    (_rm)(i, 0) = xp(i) / xp.norm();
    (_rm)(i, 1) = yp(i) / yp.norm();
    (_rm)(i, 2) = zp(i) / zp.norm();
  }

  return _rm;
}
