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

#include "MoskitoSoluteSourceKernelS2.h"
#include "Function.h"

registerMooseObject("MoskitoApp", MoskitoSoluteSourceKernelS2);

InputParameters
MoskitoSoluteSourceKernelS2::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("C1", "Concentration of original color");
  params.addParam<Real>("value", 1.0, "Constant solute source (sink) (mol/s/m^3) "
        "(positive is a source, and negative is a sink) or a multiplier "
        "for the the provided function");
  params.addParam<FunctionName>("function", "1.0", "solute source (sink) as "
        "a function (mol/s/m^3) (positive is a source, and negative is a sink)");
  return params;
}

MoskitoSoluteSourceKernelS2::MoskitoSoluteSourceKernelS2(const InputParameters & parameters)
  : Kernel(parameters),
    _C1(coupledValue("C1")),
    _new_C1(getMaterialProperty<Real>("new_C1")),
    _n(getMaterialProperty<Real>("porosity")),
    _scale(getParam<Real>("value")),
    _function(getFunction("function")),
    _K(getMaterialProperty<Real>("color_indicator")),
    _scale_factor(getMaterialProperty<Real>("scale_factor")),
    _SUPG_p(getMaterialProperty<RealVectorValue>("solute_petrov_supg_p_function")),
    _SUPG_ind(getMaterialProperty<bool>("solute_supg_indicator"))
{
}

Real
MoskitoSoluteSourceKernelS2::computeQpResidual()
{
  Real factor = -_scale * _function.value(_t, _q_point[_qp])  * _K[_qp] * _n[_qp] * 0.0 / _dt;
  // Real factor = -_scale * _function.value(_t, _q_point[_qp]) * (_C1[_qp] - _new_C1[_qp] ) * _K[_qp] * _n[_qp] / _dt;

  Real test = 0.0;

  if (_SUPG_ind[_qp])
    test = _test[_i][_qp] + _SUPG_p[_qp] * _grad_test[_i][_qp];
  else
    test = _test[_i][_qp];

  return _scale_factor[_qp] * test * factor;
}
