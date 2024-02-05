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

#include "Material.h"

// class MoskitoFluidWellGeneral;
//
// template <>
// InputParameters validParams<MoskitoFluidWellGeneral>();

class MoskitoFluidWellGeneral : public Material
{
public:
  static InputParameters validParams();
  MoskitoFluidWellGeneral(const InputParameters & parameters);
  virtual void computeQpProperties() override;

protected:
  // Velocity in well
  MaterialProperty<Real> & _u;
  // Reynolds number in well
  MaterialProperty<Real> & _Re;
  // Moody friction coefficient
  MaterialProperty<Real> & _friction;
  // Well diameter
  MaterialProperty<Real> & _dia;
  // Well area
  MaterialProperty<Real> & _area;
  // Well wetted perimeter
  MaterialProperty<Real> & _perimeter;
  // unit vector along well towards bottomhole
  MaterialProperty<RealVectorValue> & _well_dir;
  // The gravity acceleration as a vector
  MaterialProperty<RealVectorValue> & _gravity;
  // thermal conductivity of casing and fluid
  MaterialProperty<Real> & _lambda;
  // Direction of flow, the positive sign is production and vice versa
  MaterialProperty<Real> & _well_sign;
  // Hudraulic diameter
  MaterialProperty<Real> & _H_dia;

  // The coupled pressure
  const VariableValue & _P;
  const VariableValue & _T;
  const VariableValue & _q;
  MaterialProperty<Real> & _well_P;
  MaterialProperty<Real> & _well_T;
  MaterialProperty<Real> & _well_q;

  // function to calculate friction factor using Moody chart
  void MoodyFrictionFactor(Real & friction, Real rel_roughness, Real ReNo, MooseEnum roughness_type);

  // function for calculating the unit vector of well orientation
  RealVectorValue WellUnitVector();

  // local variables
  RealVectorValue _g;
  Real _lambda0;
  Real _thickness;
  Real _d;
  Real _rel_roughness;
  Real _u_f;
  Real _u_area;
  Real _u_perimeter;
  bool _f_defined;
  bool _area_defined;
  bool _perimeter_defined;
  MooseEnum _roughness_type;
  MooseEnum _well_direction;
  MooseEnum _well_type;
  const Real PI = 3.141592653589793238462643383279502884197169399375105820974944592308;
};
