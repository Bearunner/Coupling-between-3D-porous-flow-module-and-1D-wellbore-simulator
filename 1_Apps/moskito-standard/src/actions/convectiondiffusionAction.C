//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "convectiondiffusionAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

registerMooseAction("MoskitoApp", convectiondiffusionAction, "add_kernel");

registerMooseAction("MoskitoApp", convectiondiffusionAction, "add_material");

// template <>
// InputParameters
// validParams<convectiondiffusionAction>()
// {
//   InputParameters params = validParams<Action>();

InputParameters
convectiondiffusionAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "variables", "The names of the convection and diffusion variables in the simulation");
  params.addParam<int>("index_CO2", 100, "CO2 index in variable vector");
  params.addParam<int>("index_CH4", 100, "CH4 index in variable vector");
  params.addParam<int>("index_N2", 100, "N2 index in variable vector");
  params.addParam<int>("index_H2S", 100, "H2S index in variable vector");
  params.addParam<int>("index_NaCl", 100, "NaCl index in variable vector");
  params.addParam<int>("index_KCl", 100, "KCl index in variable vector");
  params.addParam<int>("index_CaCl2", 100, "CaCl2 index in variable vector");
  params.addParam<int>("index_MgCl2", 100, "MgCl2 index in variable vector");
  params.addParam<int>("index_presure", 100, "pressure index in variable vector");
  params.addParam<int>("index_enthalpy", 100, "enthalpy index in variable vector");
  params.addParam<int>("index_massrate", 100, "mass rate index in variable vector");
  //params.addRequiredParam<MaterialPropertyName>("aux1_name", "the name of aux1");

  return params;
}

convectiondiffusionAction::convectiondiffusionAction(const InputParameters & params) : Action(params) {}

void
convectiondiffusionAction::act()
{
  std::vector<NonlinearVariableName> variables =
      getParam<std::vector<NonlinearVariableName>>("variables");

 int index_CO2 = getParam<int>("index_CO2");
 int index_CH4 = getParam<int>("index_CH4");
 int index_N2 = getParam<int>("index_N2");
 int index_H2S = getParam<int>("index_H2S");
 int index_NaCl = getParam<int>("index_NaCl");
 int index_KCl = getParam<int>("index_KCl");
 int index_CaCl2 = getParam<int>("index_CaCl2");
 int index_MgCl2 = getParam<int>("index_MgCl2");
 int index_presure = getParam<int>("index_presure");
 int index_enthalpy = getParam<int>("index_enthalpy");
 int index_massrate = getParam<int>("index_massrate");

  // Do some error checking
  mooseAssert(variables.size() == 5, "Expected 1 variables, received " << variables.size());

  // Setup our Diffusion Kernel on the "u" variable
  if (_current_task == "add_kernel")
  {
    if(index_CO2 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoTransport");
        params.set<NonlinearVariableName>("variable") = variables[index_CO2-1];
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        params.set<MaterialPropertyName>("aux1_name") = "aux1";
        params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p";
        params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h";
        params.set<MaterialPropertyName>("aux1_c_name") = "m1";
        params.set<MaterialPropertyName>("aux2_name") = "aux2";
        params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p";
        params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h";
        params.set<MaterialPropertyName>("aux2_c_name") = "m2";
        params.set<MaterialPropertyName>("u_g_c_name") = "m3";
        params.set<MaterialPropertyName>("u_l_c_name") = "m4";
        _problem->addKernel("MoskitoTransport", "dif_u", params);
      }


      {
        InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
        params.set<NonlinearVariableName>("variable") = variables[index_CO2-1];
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<MaterialPropertyName>("aux1_name") = "aux1";
        params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p";
        params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h";
        params.set<MaterialPropertyName>("aux1_c_name") = "m1";
        params.set<MaterialPropertyName>("aux2_name") = "aux2";
        params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p";
        params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h";
        params.set<MaterialPropertyName>("aux2_c_name") = "m2";
        _problem->addKernel("MoskitoTimeTransport", "dif_u_t", params);
      }
    }
    if(index_CH4 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoTransport");
        params.set<NonlinearVariableName>("variable") = variables[index_CH4-1];
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        params.set<MaterialPropertyName>("aux1_name") = "aux1_CH4";
        params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_CH4";
        params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_CH4";
        params.set<MaterialPropertyName>("aux1_c_name") = "s1";
        params.set<MaterialPropertyName>("aux2_name") = "aux2_CH4";
        params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_CH4";
        params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_CH4";
        params.set<MaterialPropertyName>("aux2_c_name") = "s2";
        params.set<MaterialPropertyName>("u_g_c_name") = "s3";
        params.set<MaterialPropertyName>("u_l_c_name") = "s4";
        _problem->addKernel("MoskitoTransport", "dif_u2", params);
      }


      {
        InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
        params.set<NonlinearVariableName>("variable") = variables[index_CH4-1];
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<MaterialPropertyName>("aux1_name") = "aux1_CH4";
        params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_CH4";
        params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_CH4";
        params.set<MaterialPropertyName>("aux1_c_name") = "s1";
        params.set<MaterialPropertyName>("aux2_name") = "aux2_CH4";
        params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_CH4";
        params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_CH4";
        params.set<MaterialPropertyName>("aux2_c_name") = "s2";
        _problem->addKernel("MoskitoTimeTransport", "dif_u_t2", params);
      }
    }
      if(index_N2 < 100)
      {
        {
          InputParameters params = _factory.getValidParams("MoskitoTransport");
          params.set<NonlinearVariableName>("variable") = variables[index_N2-1];
          params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
          params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
          params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
          params.set<MaterialPropertyName>("aux1_name") = "aux1_N2";
          params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_N2";
          params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_N2";
          params.set<MaterialPropertyName>("aux1_c_name") = "a1";
          params.set<MaterialPropertyName>("aux2_name") = "aux2_N2";
          params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_N2";
          params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_N2";
          params.set<MaterialPropertyName>("aux2_c_name") = "a2";
          params.set<MaterialPropertyName>("u_g_c_name") = "a3";
          params.set<MaterialPropertyName>("u_l_c_name") = "a4";
          _problem->addKernel("MoskitoTransport", "dif_u3", params);
        }


        {
          InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
          params.set<NonlinearVariableName>("variable") = variables[index_N2-1];
          params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
          params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
          params.set<MaterialPropertyName>("aux1_name") = "aux1_N2";
          params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_N2";
          params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_N2";
          params.set<MaterialPropertyName>("aux1_c_name") = "a1";
          params.set<MaterialPropertyName>("aux2_name") = "aux2_N2";
          params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_N2";
          params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_N2";
          params.set<MaterialPropertyName>("aux2_c_name") = "a2";
          _problem->addKernel("MoskitoTimeTransport", "dif_u_t3", params);
        }
      }
      if(index_H2S < 100)
      {
          {
            InputParameters params = _factory.getValidParams("MoskitoTransport");
            params.set<NonlinearVariableName>("variable") = variables[index_H2S-1];
            params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
            params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
            params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
            params.set<MaterialPropertyName>("aux1_name") = "aux1_H2S";
            params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_H2S";
            params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_H2S";
            params.set<MaterialPropertyName>("aux1_c_name") = "b1";
            params.set<MaterialPropertyName>("aux2_name") = "aux2_H2S";
            params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_H2S";
            params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_H2S";
            params.set<MaterialPropertyName>("aux2_c_name") = "b2";
            params.set<MaterialPropertyName>("u_g_c_name") = "b3";
            params.set<MaterialPropertyName>("u_l_c_name") = "b4";
            _problem->addKernel("MoskitoTransport", "dif_u4", params);
          }


          {
            InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
            params.set<NonlinearVariableName>("variable") = variables[index_H2S-1];
            params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
            params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
            params.set<MaterialPropertyName>("aux1_name") = "aux1_H2S";
            params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_H2S";
            params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_H2S";
            params.set<MaterialPropertyName>("aux1_c_name") = "b1";
            params.set<MaterialPropertyName>("aux2_name") = "aux2_H2S";
            params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_H2S";
            params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_H2S";
            params.set<MaterialPropertyName>("aux2_c_name") = "b2";
            _problem->addKernel("MoskitoTimeTransport", "dif_u_t4", params);
          }
        }
        if(index_NaCl < 100)
        {
            {
              InputParameters params = _factory.getValidParams("MoskitoTransport");
              params.set<NonlinearVariableName>("variable") = variables[index_NaCl-1];
              params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
              params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
              params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
              params.set<MaterialPropertyName>("aux1_name") = "aux1_NaCl";
              params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_NaCl";
              params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_NaCl";
              params.set<MaterialPropertyName>("aux1_c_name") = "c1";
              params.set<MaterialPropertyName>("aux2_name") = "aux2_NaCl";
              params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_NaCl";
              params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_NaCl";
              params.set<MaterialPropertyName>("aux2_c_name") = "c2";
              params.set<MaterialPropertyName>("u_g_c_name") = "c3";
              params.set<MaterialPropertyName>("u_l_c_name") = "c4";
              _problem->addKernel("MoskitoTransport", "dif_u5", params);
            }


            {
              InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
              params.set<NonlinearVariableName>("variable") = variables[index_NaCl-1];
              params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
              params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
              params.set<MaterialPropertyName>("aux1_name") = "aux1_NaCl";
              params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_NaCl";
              params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_NaCl";
              params.set<MaterialPropertyName>("aux1_c_name") = "c1";
              params.set<MaterialPropertyName>("aux2_name") = "aux2_NaCl";
              params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_NaCl";
              params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_NaCl";
              params.set<MaterialPropertyName>("aux2_c_name") = "c2";
              _problem->addKernel("MoskitoTimeTransport", "dif_u_t5", params);
            }
          }
          if(index_KCl < 100)
          {
              {
                InputParameters params = _factory.getValidParams("MoskitoTransport");
                params.set<NonlinearVariableName>("variable") = variables[index_KCl-1];
                params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
                params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
                params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
                params.set<MaterialPropertyName>("aux1_name") = "aux1_KCl";
                params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_KCl";
                params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_KCl";
                params.set<MaterialPropertyName>("aux1_c_name") = "d1";
                params.set<MaterialPropertyName>("aux2_name") = "aux2_KCl";
                params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_KCl";
                params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_KCl";
                params.set<MaterialPropertyName>("aux2_c_name") = "d2";
                params.set<MaterialPropertyName>("u_g_c_name") = "d3";
                params.set<MaterialPropertyName>("u_l_c_name") = "d4";
                _problem->addKernel("MoskitoTransport", "dif_u6", params);
              }


              {
                InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
                params.set<NonlinearVariableName>("variable") = variables[index_KCl-1];
                params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
                params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
                params.set<MaterialPropertyName>("aux1_name") = "aux1_KCl";
                params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_KCl";
                params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_KCl";
                params.set<MaterialPropertyName>("aux1_c_name") = "d1";
                params.set<MaterialPropertyName>("aux2_name") = "aux2_KCl";
                params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_KCl";
                params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_KCl";
                params.set<MaterialPropertyName>("aux2_c_name") = "d2";
                _problem->addKernel("MoskitoTimeTransport", "dif_u_t6", params);
              }
            }
            if(index_CaCl2 < 100)
            {
                {
                  InputParameters params = _factory.getValidParams("MoskitoTransport");
                  params.set<NonlinearVariableName>("variable") = variables[index_CaCl2-1];
                  params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
                  params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
                  params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
                  params.set<MaterialPropertyName>("aux1_name") = "aux1_CaCl2";
                  params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_CaCl2";
                  params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_CaCl2";
                  params.set<MaterialPropertyName>("aux1_c_name") = "e1";
                  params.set<MaterialPropertyName>("aux2_name") = "aux2_CaCl2";
                  params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_CaCl2";
                  params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_CaCl2";
                  params.set<MaterialPropertyName>("aux2_c_name") = "e2";
                  params.set<MaterialPropertyName>("u_g_c_name") = "e3";
                  params.set<MaterialPropertyName>("u_l_c_name") = "e4";
                  _problem->addKernel("MoskitoTransport", "dif_u7", params);
                }


                {
                  InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
                  params.set<NonlinearVariableName>("variable") = variables[index_CaCl2-1];
                  params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
                  params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
                  params.set<MaterialPropertyName>("aux1_name") = "aux1_CaCl2";
                  params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_CaCl2";
                  params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_CaCl2";
                  params.set<MaterialPropertyName>("aux1_c_name") = "e1";
                  params.set<MaterialPropertyName>("aux2_name") = "aux2_CaCl2";
                  params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_CaCl2";
                  params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_CaCl2";
                  params.set<MaterialPropertyName>("aux2_c_name") = "e2";
                  _problem->addKernel("MoskitoTimeTransport", "dif_u_t7", params);
                }
              }
              if(index_MgCl2 < 100)
              {
                  {
                    InputParameters params = _factory.getValidParams("MoskitoTransport");
                    params.set<NonlinearVariableName>("variable") = variables[index_MgCl2-1];
                    params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
                    params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
                    params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
                    params.set<MaterialPropertyName>("aux1_name") = "aux1_MgCl2";
                    params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_MgCl2";
                    params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_MgCl2";
                    params.set<MaterialPropertyName>("aux1_c_name") = "f1";
                    params.set<MaterialPropertyName>("aux2_name") = "aux2_MgCl2";
                    params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_MgCl2";
                    params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_MgCl2";
                    params.set<MaterialPropertyName>("aux2_c_name") = "f2";
                    params.set<MaterialPropertyName>("u_g_c_name") = "f3";
                    params.set<MaterialPropertyName>("u_l_c_name") = "f4";
                    _problem->addKernel("MoskitoTransport", "dif_u8", params);
                  }


                  {
                    InputParameters params = _factory.getValidParams("MoskitoTimeTransport");
                    params.set<NonlinearVariableName>("variable") = variables[index_MgCl2-1];
                    params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
                    params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
                    params.set<MaterialPropertyName>("aux1_name") = "aux1_MgCl2";
                    params.set<MaterialPropertyName>("aux1_p_name") = "aux1_p_MgCl2";
                    params.set<MaterialPropertyName>("aux1_h_name") = "aux1_h_MgCl2";
                    params.set<MaterialPropertyName>("aux1_c_name") = "f1";
                    params.set<MaterialPropertyName>("aux2_name") = "aux2_MgCl2";
                    params.set<MaterialPropertyName>("aux2_p_name") = "aux2_p_MgCl2";
                    params.set<MaterialPropertyName>("aux2_h_name") = "aux2_h_MgCl2";
                    params.set<MaterialPropertyName>("aux2_c_name") = "f2";
                    _problem->addKernel("MoskitoTimeTransport", "dif_u_t8", params);
                  }
                }

    ////////////////////////////////////////////////////////////////////////////

  }
  else if (_current_task == "add_material")
  {
    if(index_CO2 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "m1";
        params.set<MaterialPropertyName>("aux2_c_name") = "m2";
        params.set<MaterialPropertyName>("u_g_c_name") = "m3";
        params.set<MaterialPropertyName>("u_l_c_name") = "m4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "m5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "m6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "m7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "m8";
        params.set<MaterialPropertyName>("T_cc_name") = "m9";
        params.set<int>("component") = 0;
        params.set<std::vector<std::string>>("output_properties") = {"m1","m2","m3","m4","m5","m6","m7","m8","m9"};
        params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat1", params);
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    if(index_CH4 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "s1";
        params.set<MaterialPropertyName>("aux2_c_name") = "s2";
        params.set<MaterialPropertyName>("u_g_c_name") = "s3";
        params.set<MaterialPropertyName>("u_l_c_name") = "s4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "s5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "s6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "s7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "s8";
        params.set<MaterialPropertyName>("T_cc_name") = "s9";
        params.set<int>("component") = 1;
        params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat2", params);
      }
    }
    if(index_N2 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "a1";
        params.set<MaterialPropertyName>("aux2_c_name") = "a2";
        params.set<MaterialPropertyName>("u_g_c_name") = "a3";
        params.set<MaterialPropertyName>("u_l_c_name") = "a4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "a5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "a6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "a7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "a8";
        params.set<MaterialPropertyName>("T_cc_name") = "a9";
        params.set<int>("component") = 2;
        //params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        //params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat3", params);
      }
    }
    if(index_H2S < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "b1";
        params.set<MaterialPropertyName>("aux2_c_name") = "b2";
        params.set<MaterialPropertyName>("u_g_c_name") = "b3";
        params.set<MaterialPropertyName>("u_l_c_name") = "b4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "b5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "b6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "b7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "b8";
        params.set<MaterialPropertyName>("T_cc_name") = "b9";
        params.set<int>("component") = 3;
        //params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        //params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat4", params);
      }
    }
    if(index_NaCl < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "c1";
        params.set<MaterialPropertyName>("aux2_c_name") = "c2";
        params.set<MaterialPropertyName>("u_g_c_name") = "c3";
        params.set<MaterialPropertyName>("u_l_c_name") = "c4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "c5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "c6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "c7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "c8";
        params.set<MaterialPropertyName>("T_cc_name") = "c9";
        params.set<int>("component") = 4;
        //params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        //params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat5", params);
      }
    }
    if(index_KCl < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "d1";
        params.set<MaterialPropertyName>("aux2_c_name") = "d2";
        params.set<MaterialPropertyName>("u_g_c_name") = "d3";
        params.set<MaterialPropertyName>("u_l_c_name") = "d4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "d5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "d6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "d7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "d8";
        params.set<MaterialPropertyName>("T_cc_name") = "d9";
        params.set<int>("component") = 5;
        //params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        //params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat6", params);
      }
    }
    if(index_CaCl2 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        if(index_MgCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        }
        params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "e1";
        params.set<MaterialPropertyName>("aux2_c_name") = "e2";
        params.set<MaterialPropertyName>("u_g_c_name") = "e3";
        params.set<MaterialPropertyName>("u_l_c_name") = "e4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "e5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "e6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "e7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "e8";
        params.set<MaterialPropertyName>("T_cc_name") = "e9";
        params.set<int>("component") = 6;
        //params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        //params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat7", params);
      }
    }
    if(index_MgCl2 < 100)
    {
      {
        InputParameters params = _factory.getValidParams("MoskitoComponent");
        params.set<std::vector<VariableName>>("pressure") = {variables[index_presure-1]};
        params.set<std::vector<VariableName>>("enthalpy") = {variables[index_enthalpy-1]};
        params.set<std::vector<VariableName>>("massrate") = {variables[index_massrate-1]};
        if(index_CO2 < 100)
        {
          params.set<std::vector<VariableName>>("concentration") = {variables[index_CO2-1]};
        }
        if(index_N2 < 100)
        {
          params.set<std::vector<VariableName>>("c3") = {variables[index_N2-1]};
        }
        if(index_CH4 < 100)
        {
          params.set<std::vector<VariableName>>("c2") = {variables[index_CH4-1]};
        }
        if(index_H2S < 100)
        {
          params.set<std::vector<VariableName>>("c4") = {variables[index_H2S-1]};
        }
        if(index_KCl < 100)
        {
          params.set<std::vector<VariableName>>("c6") = {variables[index_KCl-1]};
        }
        if(index_CaCl2 < 100)
        {
          params.set<std::vector<VariableName>>("c7") = {variables[index_CaCl2-1]};
        }
        if(index_NaCl < 100)
        {
          params.set<std::vector<VariableName>>("c5") = {variables[index_NaCl-1]};
        }
        params.set<std::vector<VariableName>>("c8") = {variables[index_MgCl2-1]};
        //params.set<std::vector<VariableName>>("c2") = {variables[4]};
        //params.set<std::vector<VariableName>>("c3") = {variables[5]};
        params.set<UserObjectName>("eos_uo") = "eos";
        params.set<UserObjectName>("viscosity_uo") = "viscosity_2p";
        params.set<UserObjectName>("drift_flux_uo") = "df";
        params.set<MaterialPropertyName>("aux1_c_name") = "f1";
        params.set<MaterialPropertyName>("aux2_c_name") = "f2";
        params.set<MaterialPropertyName>("u_g_c_name") = "f3";
        params.set<MaterialPropertyName>("u_l_c_name") = "f4";
        params.set<MaterialPropertyName>("drho_m_dcc_name") = "f5";
        params.set<MaterialPropertyName>("dgamma_dcc_name") = "f6";
        params.set<MaterialPropertyName>("dkappa_dcc_name") = "f7";
        params.set<MaterialPropertyName>("domega_dcc_name") = "f8";
        params.set<MaterialPropertyName>("T_cc_name") = "f9";
        params.set<int>("component") = 7;
        //params.set<std::vector<std::string>>("output_properties") = {"s1","s2","s3","s4","s5","s6","s7","s8","s9"};
        //params.set<std::vector<OutputName>>("outputs") = {"exodus"};
        _problem->addMaterial("MoskitoComponent", "mat8", params);
      }
    }
  }


}
