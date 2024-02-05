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

#ifndef MOSKITODRIFTFLUX_H
#define MOSKITODRIFTFLUX_H

#include "GeneralUserObject.h"
#include "MoskitoDFGVar.h"

// class MoskitoDriftFlux;
//
// template <>
// InputParameters validParams<MoskitoDriftFlux>();

class MoskitoDriftFlux : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MoskitoDriftFlux(const InputParameters & parameters);
  virtual ~MoskitoDriftFlux();

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  virtual void DFMCalculator(MoskitoDFGVar & input) const = 0;
protected:
  // Convert from Si units to American system
  const Real m_to_ft    = 3.2808398;
  const Real m3_to_ft3  = 35.3146667;
  const Real kg_to_lbm  = 2.2046225;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
};

#endif /* MOSKITODRIFTFLUX_H */
