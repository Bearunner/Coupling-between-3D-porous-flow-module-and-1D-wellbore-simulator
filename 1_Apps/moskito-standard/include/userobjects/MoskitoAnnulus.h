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

#pragma once

#include "GeneralUserObject.h"

// class MoskitoAnnulus;
//
// template <>
// InputParameters validParams<MoskitoAnnulus>();

class MoskitoAnnulus : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MoskitoAnnulus(const InputParameters & parameters);

  virtual void execute();
  virtual void initialize();
  virtual void finalize();

  Real RadiativeHTCoefficient(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const;
  Real ConvectiveHTCoefficient(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const;
  Real SurfaceTemperature(const Real & T0, const Real & fac, const Real & deltaT) const;
  void CheckValidity(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const;

protected:
  Real GrashofNo(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const;
  Real RayleighNo(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const;
  Real NusseltNo(const Real & ri, const Real & ro, const Real & Ti, const Real & To) const;

  const Real & _eo;
  const Real & _ei;
  const Real & _mu;
  const Real & _lambda;
  const Real & _cp;
  const Real & _beta;
  const Real & _rho;
  const Real & _g;
  MooseEnum _hc_method;
  enum HC_Cases {Dropkin_Sommerscales, Raithby_Hollands, Churchill};


  Real _Pr;
  Real _alpha;;
  const Real _Boltz = 5.670374419e-8;

};
