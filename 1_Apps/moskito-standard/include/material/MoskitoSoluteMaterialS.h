/**************************************************************************/
/*  Moskito - Hydro-thermal sImulator GEothermal Reservoirs                 */
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

#ifndef MOSKITOSOLUTEMATERIALS_H
#define MOSKITOSOLUTEMATERIALS_H

#include "RankTwoTensor.h"
#include "Material.h"
#include "MoskitoSUPG.h"
#include "Function.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

class MoskitoSoluteMaterialS : public Material
{
public:
  static InputParameters validParams();
  MoskitoSoluteMaterialS(const InputParameters & parameters);

private:
  // enum to select type of advection velocity
  enum AT {pure_diffusion, well_velocity, user_velocity, well_user_velocities};
  MooseEnum _at;
  const VariableValue & _T;
  const VariableValue & _C1;
  const VariableValue & _C2;
  const VariableValue & _C3;
  const VariableValue & _C4;
  const VariableValue & _C5;
  const VariableValue & _C6;
  Real _C1_T_threshold;
  Real _C3_T_threshold;
  Real _C5_T_threshold;
  MaterialProperty<Real> & _new_C1;
  MaterialProperty<Real> & _new_C2;
  MaterialProperty<Real> & _new_C3;
  MaterialProperty<Real> & _new_C4;
  MaterialProperty<Real> & _new_C5;
  MaterialProperty<Real> & _new_C6;
  MaterialProperty<Real> & _C1_indicator;
  MaterialProperty<Real> & _C3_indicator;
  MaterialProperty<Real> & _C5_indicator;
  // boolean selecting mode for upwinding and critical numbers output
  bool _has_PeCr;
  bool _has_supg;
  // userdefined factor to manually modify upwinding coefficient
  Real _supg_scale;
  // userdefined velocity vector function for advection
  const Function * _vel_func;

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties() override;
  RankTwoTensor DispersionTensorCalculator(const RealVectorValue & well_v, Real const & dispersion_l, Real const & dispersion_tr, int dim, int dimMesh, Real diffusion_factor);
  // computes scaling factor for lower dimensional elements
  Real Scaling();
  RankTwoTensor lowerDRotationMatrix(int dim);

  // Peclet number upon request
  MaterialProperty<Real> * _Pe;
  // Courant number upon request
  MaterialProperty<Real> * _Cr;
  // coefficient for solute time kernel
  MaterialProperty<Real> & _TimeKernelS;
  // indicator to inform kernels for considering upwinding
  MaterialProperty<bool> & _SUPG_ind;
  // indicator to inform kernels for considering derivative of well velocity
  MaterialProperty<bool> & _av_ind;
  // advection velocity
  MaterialProperty<RealVectorValue> & _av;
  // upwinding coefficient
  MaterialProperty<RealVectorValue> & _SUPG_p;
  // Material for rotation matrix for local cordinates
  MaterialProperty<RankTwoTensor> & _rot_mat;
  // imported well velocity
  const MaterialProperty<RealVectorValue> * _wv;

  // userobject to calculate upwinding
  const MoskitoSUPG * _supg_uo;

  // molecular diffusion as input parameter
  Real _diffusion_molecular;
  // initial molecular dispersivity
  Real _disp_l;
  Real _disp_t;
  // Formation factor
  Real _formation_factor;
  // MaterialProperty<Real> & _K;
  // const Function & _initial_temperature;
  // MaterialProperty<Real> & _T_ini;
  // const Function & _initial_C1;
  // MaterialProperty<Real> & _C1_ini;
  // MaterialProperty<Real> & _C1_indicator_ini;
  // const Function & _initial_C3;
  // MaterialProperty<Real> & _C3_ini;
  // MaterialProperty<Real> & _C3_indicator_ini;
  // const Function & _initial_C5;
  // MaterialProperty<Real> & _C5_ini;
  // MaterialProperty<Real> & _C5_indicator_ini;
  Real _n0;
  MaterialProperty<Real> & _n;
  // Tensors for internal calculation, will be transferred
  RankTwoTensor _dispersion_ten = RankTwoTensor();
  RankTwoTensor _diffusion_ten = RankTwoTensor();
  // Tensor for Handover to the Kernels and output as AuxVariables
  MaterialProperty<RankTwoTensor> & _dispersion_tensor;
  // Relativ Diffusion depending on porosity
  MaterialProperty<Real> & _diffusion_factor;
  // Tensor for Handover of combined diffusion and dispersion to Kernels
  MaterialProperty<RankTwoTensor> & _diffdisp;

  // Neumann number Fo
  MaterialProperty<Real> & _Fo;
  // Peclet number
  MaterialProperty<Real> & _PeDisp;

  // scaling factor
  MaterialProperty<Real> & _scale_factor;
  // Initial scaling factor
  const Function & _scale_factor0;

};

#endif /* MOSKITOSOLUTEMATERIALS_H */
