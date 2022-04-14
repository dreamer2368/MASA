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
MASA::ternary_2d_periodic<Scalar>::ternary_2d_periodic()
{
  this->mmsname = "ternary_2d_periodic";
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

  this->register_var("rho0", &rho0);
  this->register_var("drhox", &drhox);
  this->register_var("drhoy", &drhoy);
  this->register_var("krhox", &krhox);
  this->register_var("krhoy", &krhoy);
  this->register_var("offset_rhox", &offset_rhox);
  this->register_var("offset_rhoy", &offset_rhoy);

  this->register_var("Y0", &Y0);
  this->register_var("dY0x", &dY0x);
  this->register_var("dY0y", &dY0y);
  this->register_var("kx0", &kx0);
  this->register_var("ky0", &ky0);
  this->register_var("offset_x0", &offset_x0);
  this->register_var("offset_y0", &offset_y0);

  this->register_var("Y1", &Y1);
  this->register_var("dY1x", &dY1x);
  this->register_var("dY1y", &dY1y);
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
int MASA::ternary_2d_periodic<Scalar>::init_var()
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

  err += this->set_var("rho0", 1.38);
  err += this->set_var("drhox", 1.38);
  err += this->set_var("drhoy", 1.38);
  err += this->set_var("krhox", 1.38);
  err += this->set_var("krhoy", 1.38);
  err += this->set_var("offset_rhox", 1.38);
  err += this->set_var("offset_rhoy", 1.38);

  err += this->set_var("Y0", 1.38);
  err += this->set_var("dY0x", 1.38);
  err += this->set_var("dY0y", 1.38);
  err += this->set_var("kx0", 1.38);
  err += this->set_var("ky0", 1.38);
  err += this->set_var("offset_x0", 1.38);
  err += this->set_var("offset_y0", 1.38);

  err += this->set_var("Y1", 1.38);
  err += this->set_var("dY1x", 1.38);
  err += this->set_var("dY1y", 1.38);
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
void MASA::ternary_2d_periodic<Scalar>::eval_q_state(Scalar x1,Scalar y1,std::vector<Scalar> &source)
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

  ADScalar rho = eval_exact_rho(x, y);
  ADScalar YE = eval_exact_YE(x, y);
  ADScalar YI = eval_exact_YI(x, y);
  ADScalar T = eval_exact_T(x, y);

  ADScalar YA = 1.0 - YE - YI;
  ADScalar nA = rho * YA / mA;
  ADScalar nI = rho * YI / mI;
  ADScalar nE = rho * YE / mE;
  ADScalar nTotal = nA + nI + nE;
  ADScalar XI = nI / nTotal;
  ADScalar XE = nE / nTotal;
  ADScalar XA = nA / nTotal;

  ADScalar p = nTotal * R * T;
  ADScalar heatCapacity = CV_I * nI + CV_E * nE + CV_A * nA;
  ADScalar rhoE = 0.5 * rho * (U.dot(U)) + heatCapacity * T + nI * formEnergy_I;

  NumberVector<NDIM, ADScalar> rhoU;
  rhoU[0] = rho * U[0];
  rhoU[1] = rho * U[1];

  source[0] = raw_value(divergence(rhoU));

  NumberVector<NDIM, ADScalar> gradXA = XA.derivatives();
  NumberVector<NDIM, ADScalar> gradXI = XI.derivatives();
  NumberVector<NDIM, ADScalar> gradXE = XE.derivatives();
  NumberVector<NDIM, ADScalar> V_A1 = - D_A / XA * gradXA;
  NumberVector<NDIM, ADScalar> V_I1 = - D_I / XI * gradXI;
  NumberVector<NDIM, ADScalar> V_E1 = - D_E / XE * gradXE;
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
void MASA::ternary_2d_periodic<Scalar>::eval_exact_state(Scalar x,Scalar y,std::vector<Scalar> &state)
{
  state.resize(6);

  Scalar exact_u = eval_exact_u(x, y);
  Scalar exact_v = eval_exact_v(x, y);
  Scalar exact_rho = eval_exact_rho(x, y);
  Scalar exact_YE = eval_exact_YE(x, y);
  Scalar exact_YI = eval_exact_YI(x, y);
  Scalar exact_T = eval_exact_T(x, y);

  Scalar exact_YA = 1.0 - exact_YE - exact_YI;
  Scalar exact_nA = exact_rho * exact_YA / mA;
  Scalar exact_nI = exact_rho * exact_YI / mI;
  Scalar exact_nE = exact_rho * exact_YE / mE;

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
inputScalar MASA::ternary_2d_periodic<Scalar>::eval_exact_u(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_u = u0 + dux * sin(2.0 * pi * kux * (x / Lx - offset_ux))
                        + duy * sin(2.0 * pi * kuy * (y / Ly - offset_uy));

  return exact_u;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::ternary_2d_periodic<Scalar>::eval_exact_v(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_v = v0 + dvx * sin(2.0 * pi * kvx * (x / Lx - offset_vx))
                           + dvy * sin(2.0 * pi * kvy * (y / Ly - offset_vy));

  return exact_v;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::ternary_2d_periodic<Scalar>::eval_exact_rho(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_rho = rho0 + drhox * cos(2.0 * pi * krhox * (x / Lx - offset_rhox))
                               + drhoy * cos(2.0 * pi * krhoy * (y / Ly - offset_rhoy));

  return exact_rho;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::ternary_2d_periodic<Scalar>::eval_exact_YI(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_YI = Y0 + dY0x * cos(2.0 * pi * kx0 * (x / Lx - offset_x0))
                            + dY0y * cos(2.0 * pi * ky0 * (y / Ly - offset_y0));

  return exact_YI;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::ternary_2d_periodic<Scalar>::eval_exact_YE(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_YE = Y1 + dY1x * cos(2.0 * pi * kx1 * (x / Lx - offset_x1))
                            + dY1y * cos(2.0 * pi * ky1 * (y / Ly - offset_y1));

  return exact_YE;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::ternary_2d_periodic<Scalar>::eval_exact_T(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_T = T0 + dTx * cos(2.0 * pi * kTx * (x / Lx - offset_Tx))
                           + dTy * cos(2.0 * pi * kTy * (y / Ly - offset_Ty));

  return exact_T;
}

/* ------------------------------------------------
 *
 *         2D PERIODIC AMBIPOLAR TERNARY MIXTURE
 *
 *
 *
 * -----------------------------------------------
 */

template <typename Scalar>
MASA::ternary_2d_periodic_ambipolar<Scalar>::ternary_2d_periodic_ambipolar()
{
  this->mmsname = "ternary_2d_periodic_ambipolar";
  // init defaults
  this->init_var();
}

template <typename Scalar>
void MASA::ternary_2d_periodic_ambipolar<Scalar>::eval_q_state(Scalar x1,Scalar y1,std::vector<Scalar> &source)
{
  using std::cos;
  using std::sin;

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef SecondDerivType ADScalar;

  const ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  const ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());

  source.resize(5);

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;
  U[0] = this->eval_exact_u(x, y);
  U[1] = this->eval_exact_v(x, y);

  ADScalar rho = this->eval_exact_rho(x, y);
  // ADScalar YE = eval_exact_YE(x, y);
  ADScalar YI = this->eval_exact_YI(x, y);
  ADScalar T = this->eval_exact_T(x, y);

  ADScalar nI = rho * YI / this->mI;
  ADScalar nE = nI;
  ADScalar YE = this->mE * nI / rho;
  ADScalar nA = (rho - (this->mI + this->mE) * nI) / this->mA;
  ADScalar YA = 1.0 - (this->mI + this->mE) * nI / rho;
  ADScalar nTotal = nA + nI + nE;
  ADScalar XI = nI / nTotal;
  ADScalar XE = nE / nTotal;
  ADScalar XA = nA / nTotal;

  ADScalar p = nTotal * this->R * T;
  ADScalar heatCapacity = this->CV_I * nI + this->CV_E * nE + this->CV_A * nA;
  ADScalar rhoE = 0.5 * rho * (U.dot(U)) + heatCapacity * T + nI * this->formEnergy_I;

  NumberVector<NDIM, ADScalar> rhoU;
  rhoU[0] = rho * U[0];
  rhoU[1] = rho * U[1];

  source[0] = raw_value(divergence(rhoU));

  NumberVector<NDIM, ADScalar> gradXA = XA.derivatives();
  NumberVector<NDIM, ADScalar> gradXI = XI.derivatives();
  NumberVector<NDIM, ADScalar> gradXE = XE.derivatives();
  NumberVector<NDIM, ADScalar> V_A1 = - this->D_A / XA * gradXA;
  NumberVector<NDIM, ADScalar> V_I2 = - this->D_I / XI * gradXI;
  NumberVector<NDIM, ADScalar> V_E2 = - this->D_E / XE * gradXE;

  ADScalar mob_I = this->qe / this->kB * this->ZI / T * this->D_I;
  ADScalar mob_E = this->qe / this->kB * this->ZE / T * this->D_E;
  ADScalar mho = mob_I * nI * this->ZI + mob_E * nI * this->ZE;

  NumberVector<NDIM, ADScalar> ambE = - (V_I2 * this->ZI + V_E2 * this->ZE) * nI / mho;
  NumberVector<NDIM, ADScalar> V_I1 = V_I2 + mob_I * ambE;
  NumberVector<NDIM, ADScalar> V_E1 = V_E2 + mob_E * ambE;

  NumberVector<NDIM, ADScalar> Vc = YI * V_I1 + YE * V_E1 + YA * V_A1;
  NumberVector<NDIM, ADScalar> V_I = V_I1 - Vc;
  NumberVector<NDIM, ADScalar> V_E = V_E1 - Vc;
  NumberVector<NDIM, ADScalar> V_A = V_A1 - Vc;

  source[4] = raw_value(divergence(this->mI * nI * (U + V_I)));
  // source[5] = raw_value(divergence(mE * nE * (U + V_E)));

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar>> Identity = NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar>> Tau = this->mu * (GradU + transpose(GradU))
                                                         + (this->muB - 2./3. * this->mu)*divergence(U)*Identity;

  NumberVector<NDIM, Scalar> source_rhoU = raw_value(divergence(rho*U.outerproduct(U) - Tau) + p.derivatives());
  source[1] = source_rhoU[0];
  source[2] = source_rhoU[1];

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = - this->k_heat * T.derivatives()
                                   + nI * (this->CP_I * T + this->formEnergy_I) * V_I
                                   + nE * this->CP_E * T * V_E
                                   + nA * this->CP_A * T * V_A;

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
void MASA::ternary_2d_periodic_ambipolar<Scalar>::eval_exact_state(Scalar x,Scalar y,std::vector<Scalar> &state)
{
  state.resize(5);

  Scalar exact_u = this->eval_exact_u(x, y);
  Scalar exact_v = this->eval_exact_v(x, y);
  Scalar exact_rho = this->eval_exact_rho(x, y);
  // Scalar exact_YE = eval_exact_YE(x, y); // not used for ambipolar case.
  Scalar exact_YI = this->eval_exact_YI(x, y);
  Scalar exact_T = this->eval_exact_T(x, y);

  Scalar exact_nI = exact_rho * exact_YI / this->mI;
  Scalar exact_nE = exact_nI;
  Scalar exact_nA = (exact_rho - (this->mI + this->mE) * exact_nI) / this->mA;
  Scalar exact_YA = 1.0 - (this->mI + this->mE) *exact_nI / exact_rho;

  Scalar heatCapacity = this->CV_I * exact_nI + this->CV_E * exact_nE + this->CV_A * exact_nA;
  Scalar exact_rhoE = 0.5 * exact_rho * (exact_u * exact_u + exact_v * exact_v)
                      + heatCapacity * exact_T + exact_nI * this->formEnergy_I;

  state[0] = exact_rho;
  state[1] = exact_rho * exact_u;
  state[2] = exact_rho * exact_v;
  state[3] = exact_rhoE;
  state[4] = this->mI * exact_nI;
  // state[5] = mE * exact_nE;

  return;
}

/* ------------------------------------------------
 *
 *         2D 2T PERIODIC AMBIPOLAR TERNARY MIXTURE
 *
 *
 *
 * -----------------------------------------------
 */

template <typename Scalar>
MASA::ternary_2d_2t_periodic_ambipolar<Scalar>::ternary_2d_2t_periodic_ambipolar()
{
  this->mmsname = "ternary_2d_2t_periodic_ambipolar";

  this->register_var("TE0", &TE0);
  this->register_var("dTEx", &dTEx);
  this->register_var("dTEy", &dTEy);
  this->register_var("kTEx", &kTEx);
  this->register_var("kTEy", &kTEy);
  this->register_var("offset_TEx", &offset_TEx);
  this->register_var("offset_TEy", &offset_TEy);

  this->register_var("k_E", &k_E);

  this->register_var("nu_I", &nu_I);
  this->register_var("nu_A", &nu_A);

  // init defaults
  this->init_var();
}

template <typename Scalar>
int MASA::ternary_2d_2t_periodic_ambipolar<Scalar>::init_var()
{
  int err = 0;

  err += MASA::ternary_2d_periodic<Scalar>::init_var();

  err += this->set_var("TE0", 1.38);
  err += this->set_var("dTEx", 1.38);
  err += this->set_var("dTEy", 1.38);
  err += this->set_var("kTEx", 1.38);
  err += this->set_var("kTEy", 1.38);
  err += this->set_var("offset_TEx", 1.38);
  err += this->set_var("offset_TEy", 1.38);

  err += this->set_var("k_E", 1.38);

  err += this->set_var("nu_I", 1.38);
  err += this->set_var("nu_A", 1.38);

  return err;
}

template <typename Scalar>
void MASA::ternary_2d_2t_periodic_ambipolar<Scalar>::eval_q_state(Scalar x1,Scalar y1,std::vector<Scalar> &source)
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
  U[0] = this->eval_exact_u(x, y);
  U[1] = this->eval_exact_v(x, y);

  ADScalar rho = this->eval_exact_rho(x, y);
  // ADScalar YE = eval_exact_YE(x, y);
  ADScalar YI = this->eval_exact_YI(x, y);
  ADScalar T = this->eval_exact_T(x, y);
  ADScalar TE = eval_exact_TE(x, y);

  ADScalar nI = rho * YI / this->mI;
  ADScalar nE = nI;
  ADScalar YE = this->mE * nI / rho;
  ADScalar nA = (rho - (this->mI + this->mE) * nI) / this->mA;
  ADScalar YA = 1.0 - (this->mI + this->mE) * nI / rho;
  ADScalar nTotal = nA + nI + nE;
  ADScalar XI = nI / nTotal;
  ADScalar XE = nE / nTotal;
  ADScalar XA = nA / nTotal;

  ADScalar pe = nE * this->R * TE;
  ADScalar p = (nA + nI) * this->R * T + pe;
  ADScalar Ch = this->CV_I * nI + this->CV_A * nA;
  ADScalar Ue = this->CV_E * nE * TE;
  ADScalar rhoE = 0.5 * rho * (U.dot(U)) + Ch * T + Ue + nI * this->formEnergy_I;

  NumberVector<NDIM, ADScalar> rhoU;
  rhoU[0] = rho * U[0];
  rhoU[1] = rho * U[1];

  source[0] = raw_value(divergence(rhoU));

  NumberVector<NDIM, ADScalar> gradXA = XA.derivatives();
  NumberVector<NDIM, ADScalar> gradXI = XI.derivatives();
  NumberVector<NDIM, ADScalar> gradXE = XE.derivatives();
  NumberVector<NDIM, ADScalar> V_A1 = - this->D_A / XA * gradXA;
  NumberVector<NDIM, ADScalar> V_I2 = - this->D_I / XI * gradXI;
  NumberVector<NDIM, ADScalar> V_E2 = - this->D_E / XE * gradXE;

  ADScalar mob_I = this->qe / this->kB * this->ZI / T * this->D_I;
  ADScalar mob_E = this->qe / this->kB * this->ZE / TE * this->D_E;
  ADScalar mho = mob_I * nI * this->ZI + mob_E * nI * this->ZE;

  NumberVector<NDIM, ADScalar> ambE = - (V_I2 * this->ZI + V_E2 * this->ZE) * nI / mho;
  NumberVector<NDIM, ADScalar> V_I1 = V_I2 + mob_I * ambE;
  NumberVector<NDIM, ADScalar> V_E1 = V_E2 + mob_E * ambE;

  NumberVector<NDIM, ADScalar> Vc = YI * V_I1 + YE * V_E1 + YA * V_A1;
  NumberVector<NDIM, ADScalar> V_I = V_I1 - Vc;
  NumberVector<NDIM, ADScalar> V_E = V_E1 - Vc;
  NumberVector<NDIM, ADScalar> V_A = V_A1 - Vc;

  source[4] = raw_value(divergence(this->mI * nI * (U + V_I)));
  // source[5] = raw_value(divergence(mE * nE * (U + V_E)));

  // The shear strain tensor
  NumberVector<NDIM, typename ADScalar::derivatives_type> GradU = gradient(U);

  // The identity tensor I
  NumberVector<NDIM, NumberVector<NDIM, Scalar>> Identity = NumberVector<NDIM, Scalar>::identity();

  // The shear stress tensor
  NumberVector<NDIM, NumberVector<NDIM, ADScalar>> Tau = this->mu * (GradU + transpose(GradU))
                                                         + (this->muB - 2./3. * this->mu)*divergence(U)*Identity;

  NumberVector<NDIM, Scalar> source_rhoU = raw_value(divergence(rho*U.outerproduct(U) - Tau) + p.derivatives());
  source[1] = source_rhoU[0];
  source[2] = source_rhoU[1];

  // Electron heat flux
  NumberVector<NDIM, ADScalar> qe = - k_E * TE.derivatives()
                                    + nE * this->CP_E * TE * V_E;

  // Energy transfer by elastic collision
  Scalar m_EA = this->mE * this->mA / (this->mE + this->mA) / (this->mE + this->mA);
  Scalar m_EI = this->mE * this->mI / (this->mE + this->mI) / (this->mE + this->mI);
  ADScalar Wel = - 2.0 * nE * (m_EA * nu_A + m_EI * nu_I)
                            * 1.5 * this->R * (TE - T);

  source[5] = raw_value(divergence(Ue * U + qe) + pe * divergence(U) - Wel);

  // Temperature flux
  NumberVector<NDIM, ADScalar> q = - this->k_heat * T.derivatives()
                                   + nI * (this->CP_I * T + this->formEnergy_I) * V_I
                                   + nA * this->CP_A * T * V_A;

  source[3] = raw_value(divergence((rhoE + p) * U + q + qe - Tau.dot(U)));

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
void MASA::ternary_2d_2t_periodic_ambipolar<Scalar>::eval_exact_state(Scalar x,Scalar y,std::vector<Scalar> &state)
{
  state.resize(6);

  Scalar exact_u = this->eval_exact_u(x, y);
  Scalar exact_v = this->eval_exact_v(x, y);
  Scalar exact_rho = this->eval_exact_rho(x, y);
  // Scalar exact_YE = eval_exact_YE(x, y); // not used for ambipolar case.
  Scalar exact_YI = this->eval_exact_YI(x, y);
  Scalar exact_T = this->eval_exact_T(x, y);
  Scalar exact_TE = eval_exact_TE(x, y);

  Scalar exact_nI = exact_rho * exact_YI / this->mI;
  Scalar exact_nE = exact_nI;
  Scalar exact_nA = (exact_rho - (this->mI + this->mE) * exact_nI) / this->mA;
  Scalar exact_YA = 1.0 - (this->mI + this->mE) *exact_nI / exact_rho;

  Scalar exact_Ue = this->CV_E * exact_nE * exact_TE;

  Scalar Ch = this->CV_I * exact_nI + this->CV_A * exact_nA;
  Scalar exact_rhoE = 0.5 * exact_rho * (exact_u * exact_u + exact_v * exact_v)
                      + Ch * exact_T + exact_Ue + exact_nI * this->formEnergy_I;

  state[0] = exact_rho;
  state[1] = exact_rho * exact_u;
  state[2] = exact_rho * exact_v;
  state[3] = exact_rhoE;
  state[4] = this->mI * exact_nI;
  state[5] = exact_Ue;

  return;
}

template<typename Scalar> template<typename inputScalar>
inputScalar MASA::ternary_2d_2t_periodic_ambipolar<Scalar>::eval_exact_TE(inputScalar x, inputScalar y) {
  using std::cos;
  using std::sin;

  inputScalar exact_TE = TE0 + dTEx * cos(2.0 * this->pi * kTEx * (x / this->Lx - offset_TEx))
                             + dTEy * cos(2.0 * this->pi * kTEy * (y / this->Ly - offset_TEy));

  return exact_TE;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::ternary_2d_periodic);
MASA_INSTANTIATE_ALL(MASA::ternary_2d_periodic_ambipolar);
MASA_INSTANTIATE_ALL(MASA::ternary_2d_2t_periodic_ambipolar);

#endif // HAVE_METAPHYSICL
