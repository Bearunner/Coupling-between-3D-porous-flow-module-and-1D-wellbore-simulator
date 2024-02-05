/**************************************************************************/
/*  MOSKITO - Multiphysics cOupled Simulator toolKIT for wellbOres        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of MOSKITO App                                      */
/*  Co-developed by Sebastian Held                                        */
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
#include "MoskitoAnnulus.h"
#include "NewtonIteration.h"

// class MoskitoLatHeat_Exc_Formation_1p;
//
// template <>
// InputParameters validParams<MoskitoLatHeat_Exc_Formation_1p>();

class MoskitoLatHeat_Exc_Formation_1p : public Material, public NewtonIteration
{
public:
  static InputParameters validParams();
  MoskitoLatHeat_Exc_Formation_1p(const InputParameters & parameters);
  virtual void computeQpProperties() override;

  Real ResistivityNoAnnulus(const int & begin, const int & end, const bool & hf);

  virtual Real computeReferenceResidual(const Real trial_value, const Real scalar) override;
  virtual Real computeResidual(const Real trial_value, const Real scalar) override;
  virtual Real computeDerivative(const Real trial_value, const Real scalar) override;
  virtual Real initialGuess(const Real trial_value) override;
  virtual Real minimumPermissibleValue(const Real trial_value) const override
  {
    return 1.0e-6;
  }
  virtual Real maximumPermissibleValue(const Real trial_value) const override
  {
    return 1.0e+6;
  }

protected:
  // Fluid temperature at the center of the pipe
  const VariableValue & _Tf;
  // Formation temperature at the formation-casing interface
  const VariableValue & _Twf;

  // tubing outer radius
  MaterialProperty<Real> & _Rto;
  // Well thermal resistivity
  MaterialProperty<Real> & _Uto;

  // Tolerance of finite difference derivation
  const Real _tol;

  // outer diamters of well complition
  std::vector<Real> _OD_vec;
  // conductivities vector
  std::vector<Real> _lambda_vec;

  // imported from other materials
  const MaterialProperty<Real> & _Dti;
  const MaterialProperty<RealVectorValue> & _gravity;
  const MaterialProperty<RealVectorValue> & _well_dir;
  const MaterialProperty<Real> & _Re;
  const MaterialProperty<Real> * _hf;
  const bool & _add_hf;

  // Annulus outer and inner radiuses
  Real _rai = 0.0, _rao = 0.0;
  unsigned int _annulus_loc = 0;
  bool _annulus_ind = false;
  const MoskitoAnnulus * _annulus_uo;
};
