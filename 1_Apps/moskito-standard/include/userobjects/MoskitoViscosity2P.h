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

#ifndef MOSKITOVISCOSITY2P_H
#define MOSKITOVISCOSITY2P_H

#include "GeneralUserObject.h"
#include "MoskitoViscosity1P.h"

// class MoskitoViscosity2P;
//
// template <>
// InputParameters validParams<MoskitoViscosity2P>();

  /*
  all these mixing approaches are based on M.M. Awad, Y.S. Muzychka, 2008
  "Effective property models for homogeneous two-phase flows"
  */

class MoskitoViscosity2P : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MoskitoViscosity2P(const InputParameters & parameters);
  virtual ~MoskitoViscosity2P();

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  // mixture viscosity from pressure, temperature and mass fraction
  Real mixture_mu(Real pressure, Real temperature, Real mass_fraction) const;

protected:
  // Maxwell Eucken model 1 method
  Real ME1_calc(Real mu_l, Real mu_g, Real x) const;
  // Maxwell Eucken model 2 method
  Real ME2_calc(Real mu_l, Real mu_g, Real x) const;
  // Effective Medium Theory method
  Real EMT_calc(Real mu_l, Real mu_g, Real x) const;

  // Userobject to viscosity equation for gas
  const MoskitoViscosity1P & gas;
  // Userobject to viscosity equation for liquid
  const MoskitoViscosity1P & liquid;
  // mixing approach
  MooseEnum _mt;
  enum MT {Series, Parallel, ME1, ME2, EMT, Mean_ME12};
};

#endif /* MOSKITOVISCOSITY2P_H */
