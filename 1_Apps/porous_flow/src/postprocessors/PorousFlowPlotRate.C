//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPlotRate.h"
#include "PorousFlowSumQuantity.h"

registerMooseObject("PorousFlowApp", PorousFlowPlotRate);

InputParameters
PorousFlowPlotRate::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>(
      "uo", "PorousFlowSumQuantity user object name that holds the required information");
  params.addClassDescription("Extracts the value from the PorousFlowSumQuantity UserObject");
  return params;
}

PorousFlowPlotRate::PorousFlowPlotRate(const InputParameters & parameters)
  : GeneralPostprocessor(parameters), _total_mass(getUserObject<PorousFlowSumQuantity>("uo"))
{
}

PorousFlowPlotRate::~PorousFlowPlotRate() {}

void
PorousFlowPlotRate::initialize()
{
}

void
PorousFlowPlotRate::execute()
{
}

PostprocessorValue
PorousFlowPlotRate::getValue()
{
  return _total_mass.getValue() / _dt;
}
