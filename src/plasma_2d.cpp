// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author: nick $
// $Id: euler.cpp 17621 2011-02-14 16:53:09Z nick $
//
// euler.cpp: These are the MASA class member functions and constructors
//          For the Euler Equations
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h>

#ifdef HAVE_METAPHYSICL

#include <ad_masa.h>

const unsigned int NDIM = 2;

using namespace MASA;

/* ------------------------------------------------
 *
 *         2D PERIODIC AMBIPOLAR TERNARY MIXTURE
 *
 *
 *
 * -----------------------------------------------
 */

template <typename Scalar>
MASA::periodic_ternary_2d<Scalar>::periodic_ternary_2d()
{
  this->mmsname = "periodic_ternary_2d";
  this->dimension=2;
  this->register_var("u0", &u0);
  this->register_var("dux", &dux);
  this->register_var("duy", &duy);

  this->register_var("kux", &kux);
  this->register_var("kuy", &kuy);
  this->register_var("offset_ux", &offset_ux);
  this->register_var("offset_uy", &offset_uy);

  this->register_var("v0", &v0);
  this->register_var("dvx", &dvx);
  this->register_var("dvy", &dvy);

  this->register_var("kvx", &kvx);
  this->register_var("kvy", &kvy);
  this->register_var("offset_vx", &offset_vx);
  this->register_var("offset_vy", &offset_vy);

  this->register_var("n0", &n0);

  this->register_var("X0", &X0);
  this->register_var("dX0x", &dX0x);
  this->register_var("dX0y", &dX0y);
  this->register_var("kx0", &kx0);
  this->register_var("ky0", &ky0);
  this->register_var("offset_x0", &offset_x0);
  this->register_var("offset_y0", &offset_y0);

  this->register_var("X1", &X1);
  this->register_var("dX1x", &dX1x);
  this->register_var("dX1y", &dX1y);
  this->register_var("kx1", &kx1);
  this->register_var("ky1", &ky1);
  this->register_var("offset_x1", &offset_x1);
  this->register_var("offset_y1", &offset_y1);

  this->register_var("T0", &T0);
  this->register_var("dTx", &dTx);
  this->register_var("dTy", &dTy);
  this->register_var("kTx", &kTx);
  this->register_var("kTy", &kTy);
  this->register_var("offset_Tx", &offset_Tx);
  this->register_var("offset_Ty", &offset_Ty);

  this->register_var("mA", &mA);
  this->register_var("mI", &mI);
  this->register_var("mE", &mE);

  this->register_var("R", &R);

  this->register_var("CV_A", &CV_A);
  this->register_var("CV_I", &CV_I);
  this->register_var("CV_E", &CV_E);

  this->register_var("CP_A", &CP_A);
  this->register_var("CP_I", &CP_I);
  this->register_var("CP_E", &CP_E);

  this->register_var("formEnergy_I", &formEnergy_I);

  this->register_var("Lx", &Lx);
  this->register_var("Ly", &Ly);

  this->register_var("mu", &mu);
  this->register_var("muB", &muB);
  this->register_var("k_heat", &k_heat);
  this->register_var("D_A", &D_A);
  this->register_var("D_I", &D_I);
  this->register_var("D_E", &D_E);

  this->register_var("qe", &qe);
  this->register_var("kB", &kB);

  this->register_var("ZI", &ZI);
  this->register_var("ZE", &ZE);

  // init defaults
  this->init_var();

}

template <typename Scalar>
int MASA::periodic_ternary_2d<Scalar>::init_var()
{
  int err = 0;

  // randomly generated
  err += this->set_var("u0", 1.38);
  err += this->set_var("dux", 1.38);
  err += this->set_var("duy", 1.38);

  err += this->set_var("kux", 1.38);
  err += this->set_var("kuy", 1.38);
  err += this->set_var("offset_ux", 1.38);
  err += this->set_var("offset_uy", 1.38);

  err += this->set_var("v0", 1.38);
  err += this->set_var("dvx", 1.38);
  err += this->set_var("dvy", 1.38);

  err += this->set_var("kvx", 1.38);
  err += this->set_var("kvy", 1.38);
  err += this->set_var("offset_vx", 1.38);
  err += this->set_var("offset_vy", 1.38);

  err += this->set_var("n0", 1.38);

  err += this->set_var("X0", 1.38);
  err += this->set_var("dX0x", 1.38);
  err += this->set_var("dX0y", 1.38);
  err += this->set_var("kx0", 1.38);
  err += this->set_var("ky0", 1.38);
  err += this->set_var("offset_x0", 1.38);
  err += this->set_var("offset_y0", 1.38);

  err += this->set_var("X1", 1.38);
  err += this->set_var("dX1x", 1.38);
  err += this->set_var("dX1y", 1.38);
  err += this->set_var("kx1", 1.38);
  err += this->set_var("ky1", 1.38);
  err += this->set_var("offset_x1", 1.38);
  err += this->set_var("offset_y1", 1.38);

  err += this->set_var("T0", 1.38);
  err += this->set_var("dTx", 1.38);
  err += this->set_var("dTy", 1.38);
  err += this->set_var("kTx", 1.38);
  err += this->set_var("kTy", 1.38);
  err += this->set_var("offset_Tx", 1.38);
  err += this->set_var("offset_Ty", 1.38);

  err += this->set_var("mA", 1.38);
  err += this->set_var("mI", 1.38);
  err += this->set_var("mE", 1.38);

  err += this->set_var("R", 1.38);

  err += this->set_var("CV_A", 1.38);
  err += this->set_var("CV_I", 1.38);
  err += this->set_var("CV_E", 1.38);

  err += this->set_var("CP_A", 1.38);
  err += this->set_var("CP_I", 1.38);
  err += this->set_var("CP_E", 1.38);

  err += this->set_var("formEnergy_I", 1.38);

  err += this->set_var("Lx", 1.38);
  err += this->set_var("Ly", 1.38);

  err += this->set_var("mu", 1.38);
  err += this->set_var("muB", 1.38);
  err += this->set_var("k_heat", 1.38);
  err += this->set_var("D_A", 1.38);
  err += this->set_var("D_I", 1.38);
  err += this->set_var("D_E", 1.38);

  err += this->set_var("qe", 1.38);
  err += this->set_var("kB", 1.38);

  err += this->set_var("ZI", 1.38);
  err += this->set_var("ZE", 1.38);

  return err;
}


template <typename Scalar>
void MASA::periodic_ternary_2d<Scalar>::eval_q_state(Scalar x1,Scalar y1,std::vector<Scalar> &source)
{
  using std::cos;
  using std::sin;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());

  source.resize(6);

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;
  U[0] = eval_exact_u(x, y);
  U[1] = eval_exact_v(x, y);

  ADScalar nE = eval_exact_nE(x, y);
  ADScalar nI = eval_exact_nI(x, y);
  ADScalar nA = n0 - nI - nE;
  ADScalar T = eval_exact_T(x, y);
  ADScalar p = n0 * R * T;

  ADScalar heatCapacity = CV_I * nI + CV_E * nE + CV_A * nA;

  ADScalar rho = mI * nI + mE * nE + mA * nA;
  NumberVector<NDIM, ADScalar> rhoU;
  rhoU[0] = rho * U[0];
  rhoU[1] = rho * U[1];

  ADScalar rhoE = 0.5 * rho * (U.dot(U)) + heatCapacity * T + nI * formEnergy_I;

  source[0] = raw_value(divergence(rhoU));

  ADScalar YI = mI * nI / rho;
  ADScalar YE = mE * nE / rho;
  ADScalar YA = 1.0 - YI - YE;
  NumberVector<NDIM, ADScalar> gradYA = YA.derivatives();
  NumberVector<NDIM, ADScalar> gradYI = YI.derivatives();
  NumberVector<NDIM, ADScalar> gradYE = YE.derivatives();
  NumberVector<NDIM, ADScalar> V_A1 = - D_A / YA * gradYA;
  NumberVector<NDIM, ADScalar> V_I1 = - D_I / YI * gradYI;
  NumberVector<NDIM, ADScalar> V_E1 = - D_E / YE * gradYE;
  //
  // ADScalar mob_I = qe / kB * ZI / T * D_I;
  // ADScalar mob_E = qe / kB * ZE / T * D_E;
  // ADScalar mho = mob_I * ne * ZI + mob_E * ne * ZE;
  //
  // NumberVector<NDIM, ADScalar> ambE = - (V_I2 * ZI + V_E2 * ZE) * ne / mho;
  // NumberVector<NDIM, ADScalar> V_I1 = V_I2 + mob_I * ambE;
  // NumberVector<NDIM, ADScalar> V_E1 = V_E2 + mob_E * ambE;
  //
  NumberVector<NDIM, ADScalar> Vc = YI * V_I1 + YE * V_E1 + YA * V_A1;
  NumberVector<NDIM, ADScalar> V_I = V_I1 - Vc;
  NumberVector<NDIM, ADScalar> V_E = V_E1 - Vc;
  NumberVector<NDIM, ADScalar> V_A = V_A1 - Vc;

  source[4] = raw_value(divergence(mI * nI * (U + V_I)));
  source[5] = raw_value(divergence(mE * nE * (U + V_E)));

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar>> Identity = NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar>> Tau = mu * (GradU + transpose(GradU))
                                                         + (muB - 2./3. * mu)*divergence(U)*Identity;

  NumberVector<NDIM, Scalar> source_rhoU = raw_value(divergence(rho*U.outerproduct(U) - Tau) + p.derivatives());
  source[1] = source_rhoU[0];
  source[2] = source_rhoU[1];

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = - k_heat * T.derivatives()
                                   + nI * (CP_I * T + formEnergy_I) * V_I
                                   + nE * CP_E * T * V_E
                                   + nA * CP_A * T * V_A;

  source[3] = raw_value(divergence((rhoE + p) * U + q - Tau.dot(U)));

  return;
}

/* ------------------------------------------------
 *
 *
 *   Analytical terms
 *
 * -----------------------------------------------
 */

template <typename Scalar>
void MASA::periodic_ternary_2d<Scalar>::eval_exact_state(Scalar x,Scalar y,std::vector<Scalar> &state)
{
  state.resize(6);

  Scalar exact_u = eval_exact_u(x, y);
  Scalar exact_v = eval_exact_v(x, y);
  Scalar exact_nE = eval_exact_nE(x, y);
  Scalar exact_nI = eval_exact_nI(x, y);
  Scalar exact_nA = n0 - exact_nE - exact_nI;
  Scalar exact_T = eval_exact_T(x, y);
  Scalar exact_rho = mI * exact_nI + mE * exact_nE + mA * exact_nA;

  Scalar heatCapacity = CV_I * exact_nI + CV_E * exact_nE + CV_A * exact_nA;
  Scalar exact_rhoE = 0.5 * exact_rho * (exact_u * exact_u + exact_v * exact_v)
                      + heatCapacity * exact_T + exact_nI * formEnergy_I;

  state[0] = exact_rho;
  state[1] = exact_rho * exact_u;
  state[2] = exact_rho * exact_v;
  state[3] = exact_rhoE;
  state[4] = mI * exact_nI;
  state[5] = mE * exact_nE;

  return;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::periodic_ternary_2d<Scalar>::eval_exact_u(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_u = u0 + dux * sin(2.0 * pi * kux * (x / Lx - offset_ux))
                        + duy * sin(2.0 * pi * kuy * (y / Ly - offset_uy));

  return exact_u;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::periodic_ternary_2d<Scalar>::eval_exact_v(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_v = v0 + dvx * sin(2.0 * pi * kvx * (x / Lx - offset_vx))
                           + dvy * sin(2.0 * pi * kvy * (y / Ly - offset_vy));

  return exact_v;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::periodic_ternary_2d<Scalar>::eval_exact_nI(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_nI = n0 * (X0 + dX0x * cos(2.0 * pi * kx0 * (x / Lx - offset_x0))
                                  + dX0y * cos(2.0 * pi * ky0 * (y / Ly - offset_y0)));

  return exact_nI;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::periodic_ternary_2d<Scalar>::eval_exact_nE(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_nE = n0 * (X1 + dX1x * cos(2.0 * pi * kx1 * (x / Lx - offset_x1))
                                  + dX1y * cos(2.0 * pi * ky1 * (y / Ly - offset_y1)));

  return exact_nE;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::periodic_ternary_2d<Scalar>::eval_exact_T(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_T = T0 + dTx * cos(2.0 * pi * kTx * (x / Lx - offset_Tx))
                           + dTy * cos(2.0 * pi * kTy * (y / Ly - offset_Ty));

  return exact_T;
}

// /* ------------------------------------------------
//  *
//  *         2D PERIODIC AMBIPOLAR TERNARY MIXTURE
//  *
//  *
//  *
//  * -----------------------------------------------
//  */
//
// template <typename Scalar>
// MASA::periodic_ambipolar_ternary_2d<Scalar>::periodic_ambipolar_ternary_2d()
// {
//   this->mmsname = "periodic_ambipolar_ternary_2d";
//   this->dimension=2;
//   this->register_var("u0", &u0);
//   this->register_var("dux", &dux);
//   this->register_var("duy", &duy);
//
//   this->register_var("kux", &kux);
//   this->register_var("kuy", &kuy);
//   this->register_var("offset_ux", &offset_ux);
//   this->register_var("offset_uy", &offset_uy);
//
//   this->register_var("v0", &v0);
//   this->register_var("dvx", &dvx);
//   this->register_var("dvy", &dvy);
//
//   this->register_var("kvx", &kvx);
//   this->register_var("kvy", &kvy);
//   this->register_var("offset_vx", &offset_vx);
//   this->register_var("offset_vy", &offset_vy);
//
//   this->register_var("n0", &n0);
//   this->register_var("X0", &X0);
//   this->register_var("dX", &dX);
//   this->register_var("T0", &T0);
//   this->register_var("dT", &dT);
//
//   this->register_var("mA", &mA);
//   this->register_var("mI", &mI);
//   this->register_var("mE", &mE);
//
//   this->register_var("R", &R);
//
//   this->register_var("CV_A", &CV_A);
//   this->register_var("CV_I", &CV_I);
//   this->register_var("CV_E", &CV_E);
//
//   this->register_var("CP_A", &CP_A);
//   this->register_var("CP_I", &CP_I);
//   this->register_var("CP_E", &CP_E);
//
//   this->register_var("formEnergy_I", &formEnergy_I);
//
//   this->register_var("Lx", &Lx);
//   this->register_var("Ly", &Ly);
//   this->register_var("kx", &kx);
//   this->register_var("ky", &ky);
//   this->register_var("offset_x", &offset_x);
//   this->register_var("offset_y", &offset_y);
//
//   this->register_var("kTx", &kTx);
//   this->register_var("kTy", &kTy);
//   this->register_var("offset_Tx", &offset_Tx);
//   this->register_var("offset_Ty", &offset_Ty);
//
//   this->register_var("mu", &mu);
//   this->register_var("muB", &muB);
//   this->register_var("k_heat", &k_heat);
//   this->register_var("D_A", &D_A);
//   this->register_var("D_I", &D_I);
//   this->register_var("D_E", &D_E);
//
//   this->register_var("qe", &qe);
//   this->register_var("kB", &kB);
//
//   this->register_var("ZI", &ZI);
//   this->register_var("ZE", &ZE);
//
//   // init defaults
//   this->init_var();
//
// }
//
// template <typename Scalar>
// int MASA::periodic_ambipolar_ternary_2d<Scalar>::init_var()
// {
//   int err = 0;
//
//   // randomly generated
//   err += this->set_var("u0", 1.38);
//   err += this->set_var("dux", 1.38);
//   err += this->set_var("duy", 1.38);
//
//   err += this->set_var("kux", 1.38);
//   err += this->set_var("kuy", 1.38);
//   err += this->set_var("offset_ux", 1.38);
//   err += this->set_var("offset_uy", 1.38);
//
//   err += this->set_var("v0", 1.38);
//   err += this->set_var("dvx", 1.38);
//   err += this->set_var("dvy", 1.38);
//
//   err += this->set_var("kvx", 1.38);
//   err += this->set_var("kvy", 1.38);
//   err += this->set_var("offset_vx", 1.38);
//   err += this->set_var("offset_vy", 1.38);
//
//   err += this->set_var("n0", 1.38);
//   err += this->set_var("X0", 1.38);
//   err += this->set_var("dX", 1.38);
//   err += this->set_var("T0", 1.38);
//   err += this->set_var("dT", 1.38);
//
//   err += this->set_var("mA", 1.38);
//   err += this->set_var("mI", 1.38);
//   err += this->set_var("mE", 1.38);
//
//   err += this->set_var("R", 1.38);
//
//   err += this->set_var("CV_A", 1.38);
//   err += this->set_var("CV_I", 1.38);
//   err += this->set_var("CV_E", 1.38);
//
//   err += this->set_var("CP_A", 1.38);
//   err += this->set_var("CP_I", 1.38);
//   err += this->set_var("CP_E", 1.38);
//
//   err += this->set_var("formEnergy_I", 1.38);
//
//   err += this->set_var("Lx", 1.38);
//   err += this->set_var("Ly", 1.38);
//   err += this->set_var("kx", 1.38);
//   err += this->set_var("ky", 1.38);
//   err += this->set_var("offset_x", 1.38);
//   err += this->set_var("offset_y", 1.38);
//
//   err += this->set_var("kTx", 1.38);
//   err += this->set_var("kTy", 1.38);
//   err += this->set_var("offset_Tx", 1.38);
//   err += this->set_var("offset_Ty", 1.38);
//
//   err += this->set_var("mu", 1.38);
//   err += this->set_var("muB", 1.38);
//   err += this->set_var("k_heat", 1.38);
//   err += this->set_var("D_A", 1.38);
//   err += this->set_var("D_I", 1.38);
//   err += this->set_var("D_E", 1.38);
//
//   err += this->set_var("qe", 1.38);
//   err += this->set_var("kB", 1.38);
//
//   err += this->set_var("ZI", 1.38);
//   err += this->set_var("ZE", 1.38);
//
//   return err;
// }
//
//
// template <typename Scalar>
// Scalar MASA::periodic_ambipolar_ternary_2d<Scalar>::eval_q_state(Scalar x1,Scalar y1,int eq)
// {
//   using std::cos;
//   using std::sin;
//
//   typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
//   typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
//   typedef SecondDerivType ADScalar;
//
//   const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
//   const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
//
//   // Treat velocity as a vector
//   NumberVector<NDIM, ADScalar> U;
//   U[0] = u0 + dux * sin(2.0 * pi * kux * (x / Lx - offset_ux))
//             + duy * sin(2.0 * pi * kuy * (y / Ly - offset_uy));
//   U[1] = v0 + dvx * sin(2.0 * pi * kvx * (x / Lx - offset_vx))
//             + dvy * sin(2.0 * pi * kvy * (y / Ly - offset_vy));
//   ADScalar ne = n0 * (X0 + dX * sin(2.0 * pi * kx * (x / Lx - offset_x)) * sin(2.0 * pi * ky * (y / Ly - offset_y)));
//   ADScalar T = T0 + dT * sin(2.0 * pi * kTx * (x / Lx - offset_Tx)) * sin(2.0 * pi * kTy * (y / Ly - offset_Ty));
//   ADScalar p = n0 * R * T;
//
//   ADScalar heatCapacity = (CV_I + CV_E) * ne + CV_A * (n0 - 2.0 * ne);
//
//   ADScalar rho = (mI + mE) * ne + mA * (n0 - 2.0 * ne);
//   NumberVector<NDIM, ADScalar> rhoU;
//   rhoU[0] = rho * U[0];
//   rhoU[1] = rho * U[1];
//
//   ADScalar rhoE = 0.5 * rho * (U.dot(U)) + heatCapacity * T + ne * formEnergy_I;
//
//   Scalar source_rho = raw_value(divergence(rhoU));
//
//   ADScalar YA = 1.0 - (mI + mE) * ne / rho;
//   ADScalar YI = mI * ne / rho;
//   ADScalar YE = mE * ne / rho;
//   NumberVector<NDIM, ADScalar> gradYA = YA.derivatives();
//   NumberVector<NDIM, ADScalar> gradYI = YI.derivatives();
//   NumberVector<NDIM, ADScalar> gradYE = YE.derivatives();
//   NumberVector<NDIM, ADScalar> V_A2 = - D_A * gradYA / YA;
//   NumberVector<NDIM, ADScalar> V_I2 = - D_I * gradYI / YI;
//   NumberVector<NDIM, ADScalar> V_E2 = - D_E * gradYE / YE;
//
//   ADScalar mob_I = qe / kB * ZI / T * D_I;
//   ADScalar mob_E = qe / kB * ZE / T * D_E;
//   ADScalar mho = mob_I * ne * ZI + mob_E * ne * ZE;
//
//   NumberVector<NDIM, ADScalar> ambE = - (V_I2 * ZI + V_E2 * ZE) * ne / mho;
//   NumberVector<NDIM, ADScalar> V_I1 = V_I2 + mob_I * ambE;
//   NumberVector<NDIM, ADScalar> V_E1 = V_E2 + mob_E * ambE;
//
//   NumberVector<NDIM, ADScalar> Vc = mI * ne * V_I1 + mE * ne * V_E1 + mA * (n0 - 2.0 * ne) * V_A2;
//   NumberVector<NDIM, ADScalar> V_I = V_I1 - Vc;
//   NumberVector<NDIM, ADScalar> V_E = V_E1 - Vc;
//   NumberVector<NDIM, ADScalar> V_A = V_A2 - Vc;
//
//   Scalar source_rhoYI = raw_value(divergence(mI * ne * (U + V_I)));
//
//   // The shear strain tensor
//   NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);
//
//   // The identity tensor I
//   NumberVector<NDIM, NumberVector<NDIM, Scalar>> Identity = NumberVector<NDIM, Scalar>::identity();
//
//   // The shear stress tensor
//   NumberVector<NDIM, NumberVector<NDIM, ADScalar>> Tau = mu * (GradU + transpose(GradU))
//                                                          + (muB - 2./3. * mu)*divergence(U)*Identity;
//
//   NumberVector<NDIM, Scalar> source_rhoU = raw_value(divergence(rho*U.outerproduct(U) - Tau) + p.derivatives());
//
//   // Temperature flux
//   NumberVector<NDIM, ADScalar> q = - k_heat * T.derivatives()
//                                    + ne * (CP_I * T + formEnergy_I) * V_I
//                                    + ne * CP_E * T * V_E
//                                    + (n0 - 2.0 * ne) * CP_A * T * V_A;
//
//   Scalar source_rhoE = raw_value(divergence((rhoE + p) * U + q - Tau.dot(U)));
//
//   switch (eq) {
//     case 0: // rho
//       return source_rho;
//     break;
//     case 1: // rho * u
//       return source_rhoU[0];
//     break;
//     case 2: // rho * v
//       return source_rhoU[1];
//     break;
//     case 3: // rho * E
//       return source_rhoE;
//     break;
//     case 4: // rho * Y_I
//       return source_rhoYI;
//     break;
//     default:
//       // return 0.0;
//       return -1.0;
//     break;
//   }
// }
//
// /* ------------------------------------------------
//  *
//  *
//  *   Analytical terms
//  *
//  * -----------------------------------------------
//  */
//
// template <typename Scalar>
// Scalar MASA::periodic_ambipolar_ternary_2d<Scalar>::eval_exact_state(Scalar x,Scalar y,int eq)
// {
//   using std::cos;
//   using std::sin;
//
//   Scalar exact_u = u0 + dux * sin(2.0 * pi * kux * (x / Lx - offset_ux))
//                       + duy * sin(2.0 * pi * kuy * (y / Ly - offset_uy));
//   Scalar exact_v = v0 + dvx * sin(2.0 * pi * kvx * (x / Lx - offset_vx))
//                       + dvy * sin(2.0 * pi * kvy * (y / Ly - offset_vy));
//   Scalar exact_ne = n0 * (X0 + dX * sin(2.0 * pi * kx * (x / Lx - offset_x)) * sin(2.0 * pi * ky * (y / Ly - offset_y)));
//   Scalar exact_T = T0 + dT * sin(2.0 * pi * kTx * (x / Lx - offset_Tx)) * sin(2.0 * pi * kTy * (y / Ly - offset_Ty));
//   Scalar exact_rho = (mI + mE) * exact_ne + mA * (n0 - 2.0 * exact_ne);
//
//   Scalar heatCapacity = (CV_I + CV_E) * exact_ne + CV_A * (n0 - 2.0 * exact_ne);
//   Scalar exact_rhoE = 0.5 * exact_rho * (exact_u * exact_u + exact_v * exact_v)
//                       + heatCapacity * exact_T + exact_ne * formEnergy_I;
//
//   switch (eq) {
//     case 0: // rho
//       return exact_rho;
//     break;
//     case 1: // rho * u
//       return exact_rho * exact_u;
//     break;
//     case 2: // rho * v
//       return exact_rho * exact_v;
//     break;
//     case 3: // rho * E
//       return exact_rhoE;
//     break;
//     case 4: // rho * Y_I
//       return mI * exact_ne;
//     break;
//     default:
//       return -1.0;
//     break;
//   }
// }

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::periodic_ternary_2d);
// MASA_INSTANTIATE_ALL(MASA::periodic_ambipolar_ternary_2d);

#endif // HAVE_METAPHYSICL
