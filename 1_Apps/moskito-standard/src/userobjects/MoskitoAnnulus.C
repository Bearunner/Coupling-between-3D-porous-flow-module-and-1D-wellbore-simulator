/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
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

#include "MoskitoAnnulus.h"

registerMooseObject("MoskitoApp", MoskitoAnnulus);

// template <>
// InputParameters
// validParams<MoskitoAnnulus>()
// {
//   InputParameters params = validParams<GeneralUserObject>();

InputParameters
MoskitoAnnulus::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addParam<Real>("emissivity_outer", 0.7,
        "Emissivity of outer surface");
  params.addParam<Real>("emissivity_inner", 0.7,
        "Emissivity of inner surface ()");
  params.addRequiredParam<Real>("thermal_conductivity",
        "Thermal conductivity of annulus fluid (W/(m*K))");
  params.addRequiredParam<Real>("density",
        "Density of annulus fluid (kg/m^3))");
  params.addRequiredParam<Real>("viscosity",
        "Dynamic viscosity of annulus fluid (kg/(m*s))");
  params.addRequiredParam<Real>("heat_capacity",
        "Specific Heat capacity of annulus fluid (J/(K*kg))");
  params.addRequiredParam<Real>("thermal_expansion",
        "Thermal volumetric expansion coefficient annulus fluid (1/K)");
  params.addParam<Real>("gravity", 9.81,
        "The gravity acceleration for the convective term (m/s^2)");
  MooseEnum hc_method
        ("Dropkin_Sommerscales Raithby_Hollands Churchill", "Dropkin_Sommerscales");
  params.addParam<MooseEnum>("convective_heat_method", hc_method,
        "Select a method for calculating convective heat resistivity coefficient");

  return params;
}

MoskitoAnnulus::MoskitoAnnulus(const InputParameters & parameters)
  : GeneralUserObject(parameters),
  _eo(getParam<Real>("emissivity_outer")),
  _ei(getParam<Real>("emissivity_inner")),
  _mu(getParam<Real>("viscosity")),
  _lambda(getParam<Real>("thermal_conductivity")),
  _cp(getParam<Real>("heat_capacity")),
  _beta(getParam<Real>("thermal_expansion")),
  _rho(getParam<Real>("density")),
  _g(getParam<Real>("gravity")),
  _hc_method(getParam<MooseEnum>("convective_heat_method"))
{
  _Pr = _cp * _mu / _lambda;
  _alpha = _lambda / _rho / _cp;
}

void
MoskitoAnnulus::execute(){}

void
MoskitoAnnulus::initialize(){}

void
MoskitoAnnulus::finalize(){}

Real
MoskitoAnnulus::GrashofNo(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const
{
  Real Gr = 0.0;
  Gr += std::pow(ro-ri,3.0) * _g * _rho * _rho * _beta * std::abs(Ti - To);
  Gr /= _mu * _mu;

  return Gr;
}

Real
MoskitoAnnulus::RayleighNo(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const
{
  Real Lc, Ra = 0.0;
  Lc = 2.0 * std::pow(std::log(ro/ri), 4.0/3.0);
  Lc /= std::pow(std::pow(ri, -0.6) + std::pow(ro, -0.6), 5.0/3.0);
  Ra = _rho * _g * _beta * std::pow(Lc, 3.0) * std::abs(Ti - To) / (_mu * _alpha);

  return Ra;
}

Real
MoskitoAnnulus::NusseltNo(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const
{
  Real Nu = 0.0;
  Real f_Pr = std::pow(1.0 + std::pow(0.5/_Pr, 9.0/16.0), -16.0/9.0);
  Nu = 0.364 * std::pow(RayleighNo(ri, ro, Ti, To)*f_Pr, 0.25) * std::sqrt(ro / ri);

  return Nu;
}

Real
MoskitoAnnulus::RadiativeHTCoefficient(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const
{
  Real F = (1.0 / _eo - 1.0) * ri / ro + 1.0 / _ei;
  Real hr = _Boltz / F * (Ti * Ti + To * To) * (Ti + To);

  return hr;
}

Real
MoskitoAnnulus::ConvectiveHTCoefficient(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const
{
  Real khc, hc = 0.0;

  switch (_hc_method)
  {
    case HC_Cases::Dropkin_Sommerscales:
    {
      Real Gr = GrashofNo(ri, ro, Ti, To);
      // if(Gr*_Pr>=5e4 && Gr*_Pr<=7.2e8)
        khc = 0.049 * std::pow(Gr*_Pr, 0.333) * std::pow(_Pr,0.074) * _lambda;
      // else
      //   mooseError(name(), " Dropkin_Sommerscales method is not valid because ",
      //             "Pr_No = ", _Pr, " and Gr_No = ", Gr, "are not in the validity",
      //             " range => 5e4<=Pr_No*Gr_No<=7.2e8, Ra_No<1e7");

      hc = khc / (ri * std::log(ro / ri));
    }
    break;

    case HC_Cases::Raithby_Hollands:
    {
      Real Ra = RayleighNo(ri, ro, Ti, To);
      // if(_Pr>=0.7 && _Pr<=6000 && Ra<1e7)
        khc = 0.384 * std::pow(_Pr*Ra/(0.861 + _Pr),0.25) * _lambda;
      // else
      //   mooseError(name(), " Raithby&Hollands method is not valid because either ",
      //             "Pr_No = ", _Pr, "or Ra_No = ", Ra, "are not in the validity",
      //             " range => 0.7<=Pr_No<=6000, Ra_No<1e7");

      hc = khc / (ri * std::log(ro / ri));
    }
    break;

    case HC_Cases::Churchill:
      hc = NusseltNo(ri, ro, Ti, To) * _lambda / (2.0 * ro);
    break;
  }

  return hc;
}

Real
MoskitoAnnulus::SurfaceTemperature(const Real & T0, const Real & fac, const Real & deltaT) const
{
  return T0 + fac * deltaT;
}


void
MoskitoAnnulus::CheckValidity(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const
{
  switch (_hc_method)
  {
    case HC_Cases::Dropkin_Sommerscales:
    {
      Real Gr = GrashofNo(ri, ro, Ti, To);
      if(!(Gr*_Pr>=5e4 && Gr*_Pr<=7.2e8))
        mooseError(name(), ": Dropkin_Sommerscales method is not valid because ",
                  "Pr_No = ", _Pr, " and Gr_No = ", Gr, " are not in the validity",
                  " range => 5e4<=Pr_No*Gr_No<=7.2e8, Ra_No<1e7");
    }
    break;

    case HC_Cases::Raithby_Hollands:
    {
      Real Ra = RayleighNo(ri, ro, Ti, To);
      if (!(_Pr>=0.7 && _Pr<=6000 && Ra<1e7))
        mooseError(name(), ": Raithby&Hollands method is not valid because either ",
                  "Pr_No = ", _Pr, " or Ra_No = ", Ra, " are not in the validity",
                  " range => 0.7<=Pr_No<=6000, Ra_No<1e7");
    }
  }
}
