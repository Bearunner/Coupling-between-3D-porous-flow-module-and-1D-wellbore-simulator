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

#include "MoskitoSoluteTimeKernelS.h"

registerMooseObject("MoskitoApp", MoskitoSoluteTimeKernelS);

InputParameters
MoskitoSoluteTimeKernelS::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addCoupledVar("temperature", "temperature nonlinear variable (K)");
  params.addRequiredParam<Real>("temperature_threshold", "temperature_threshold for tracer temperature-response (K)");
  return params;
}

MoskitoSoluteTimeKernelS::MoskitoSoluteTimeKernelS(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _T_old(coupledValueOld("temperature")),
    _T_threshold(getParam<Real>("temperature_threshold")),
    _u_old(valueOld()),
    _scale_factor(getMaterialProperty<Real>("scale_factor")),
    _TimeKernelS(getMaterialProperty<Real>("TimeKernel_S")),
    // _new_C5_old(getMaterialPropertyOld<Real>("new_C5")),
    _SUPG_p(getMaterialProperty<RealVectorValue>("solute_petrov_supg_p_function")),
    _SUPG_ind(getMaterialProperty<bool>("solute_supg_indicator"))
{
}

Real
MoskitoSoluteTimeKernelS::computeQpResidual()
{
  Real tracer_indicator = 0.0;

  if (_T_old[_i] < (_T_threshold - 0.15)) {
    tracer_indicator = 1.0;
  }
  else if (_T_old[_i] >=(_T_threshold - 0.15) && _T_old[_i] <= (_T_threshold + 0.15)){
    tracer_indicator = 1.0 /(1.0 + pow (2.71828182846, 50.0 * (_T_old[_i] - _T_threshold)));
  }
  else {
    tracer_indicator = 0.0;
  }

  Real test = 0.0;

  if (_SUPG_ind[_qp])
    test = _test[_i][_qp] + _SUPG_p[_qp] * _grad_test[_i][_qp];
  else
    test = _test[_i][_qp];

  return _scale_factor[_qp] * _TimeKernelS[_qp] * test * (_u[_qp] - tracer_indicator * _u_old[_qp]) / _dt;
  // return _scale_factor[_qp] * _TimeKernelS[_qp] * test * _u_dot[_qp]; ##initial
  // return _scale_factor[_qp] * _TimeKernelS[_qp] * test * (_u[_qp] - _u_old[_qp]) / _dt;
}

Real
MoskitoSoluteTimeKernelS::computeQpJacobian()
{
  Real test = 0.0, j = 0.0;

  if (_SUPG_ind[_qp])
    test = _test[_i][_qp] + _SUPG_p[_qp] * _grad_test[_i][_qp];
  else
    test = _test[_i][_qp];

  // j  = _TimeKernelS[_qp]  * _phi[_j][_qp] * _du_dot_du[_qp]; ##initial
  j  = _TimeKernelS[_qp]  * _phi[_j][_qp] / _dt;

  j *= _scale_factor[_qp] * test;

  return j;
}
