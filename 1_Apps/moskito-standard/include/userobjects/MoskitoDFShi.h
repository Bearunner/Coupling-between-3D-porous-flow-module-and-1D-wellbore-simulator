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

#ifndef MOSKITODFSHI_H
#define MOSKITODFSHI_H

#include "MoskitoDriftFlux.h"

// class MoskitoDFShi;
//
// template <>
// InputParameters validParams<MoskitoDFShi>();

// forward declaration
class MoskitoShiLVar;

class MoskitoDFShi : public MoskitoDriftFlux
{
public:
  static InputParameters validParams();
  MoskitoDFShi(const InputParameters & parameters);

  virtual void DFMCalculator(MoskitoDFGVar & input) const override;

  // Preparation of Shi individual parameters
  void Shiinitialisation(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const;
  // final calculator
  void Shicalculator(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const;

protected:
  // Calculation of Kutateladze for K-function but also v_sgf
  void cal_Kutateladze(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const;
  // Calcuation of v_c and v_sgf
  void cal_velocities(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const;
  // Calulation of parameters for C0
  Real cal_auxvar_nu(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const;
  // Calcuation of transition function which continously interpolates between bubbly and flooding film conditions
  Real cal_transition_fct(MoskitoDFGVar & input, MoskitoShiLVar & LVar) const;
  // Void fraction from Material

private:
  // Surface tension of the liquid phase (kg/s²)
  Real _surf_ten;
  // wall friction factor; Pan suggests to have a const. value of 0.008 (dimensionless)
  Real _cw = 0.008;
  // _cKu is a dimensionless factor to compute the Kutateladze number; Richter et al. set it to 75, Pan et al. had better results with 142
  Real _cKu;
  // a1 and a2 are fitting parameters for the transition function K
  Real _tran_param_a1;
  Real _tran_param_a2;
  // Fv is fitting parameter of C0; purpose is to tune the sensitivity of C0 towards liquid-gas ratio
  Real _C0_Fv;
  // _C0_cMax is a C0 fitting parameter; should be 1.0 < _C0_cMax < 1.5; apprtoaching 1.0 for larger well diameters; in ECLIPSE it is set to 1.2
  Real _C0_cMax;
  // parameterization vector (m0,n1,n2) to fit the inclination function; 1.27<m0<1.85, 0.21<n1<0.24, 0.95<n2<1.08
  // std::vector<Real> _Shi_incl_triple
};

class MoskitoShiLVar
{
public:
  // Bond number is a parameter to compute the Kutateladze number, related to the well diameter, the densities of the phases and the surface tension
  Real Bond = 0.0;
  // Kutateladzenumber used for computing C0 and the transition function K
  Real Kutateladze = 0.0;
  // Superficial flodding gas velocity, is used to copmpute C0
  Real v_sgf = 0.0;
  // charcteristic velocity, used for computing drift velocity
  Real v_c = 0.0;
};

#endif /* MOSKITODFHK_H */
