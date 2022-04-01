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

const unsigned int NDIM = 2;

using namespace MASA;

/* ------------------------------------------------
 *
 *         EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */

template <typename Scalar>
MASA::periodic_argon_ternary_2d<Scalar>::periodic_argon_ternary_2d()
{
  this->mmsname = "periodic_argon_ternary_2d";
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
  this->register_var("dX", &dX);
  this->register_var("T0", &T0);
  this->register_var("dT", &dT);

  this->register_var("mA", &mA);
  this->register_var("mI", &mI);
  this->register_var("mE", &mE);

  this->register_var("CV_A", &CV_A);
  this->register_var("CV_I", &CV_I);
  this->register_var("CV_E", &CV_E);

  this->register_var("formEnergy_I", &formEnergy_I);

  this->register_var("Lx", &Lx);
  this->register_var("Ly", &Ly);
  this->register_var("kx", &kx);
  this->register_var("ky", &ky);
  this->register_var("offset_x", &offset_x);
  this->register_var("offset_y", &offset_y);

  this->register_var("kTx", &kTx);
  this->register_var("kTy", &kTy);
  this->register_var("offset_Tx", &offset_Tx);
  this->register_var("offset_Ty", &offset_Ty);

  // init defaults
  this->init_var();

}

template <typename Scalar>
int MASA::periodic_argon_ternary_2d<Scalar>::init_var()
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
  err += this->set_var("dX", 1.38);
  err += this->set_var("T0", 1.38);
  err += this->set_var("dT", 1.38);

  err += this->set_var("mA", 1.38);
  err += this->set_var("mI", 1.38);
  err += this->set_var("mE", 1.38);

  err += this->set_var("CV_A", 1.38);
  err += this->set_var("CV_I", 1.38);
  err += this->set_var("CV_E", 1.38);

  err += this->set_var("formEnergy_I", 1.38);

  err += this->set_var("Lx", 1.38);
  err += this->set_var("Ly", 1.38);
  err += this->set_var("kx", 1.38);
  err += this->set_var("ky", 1.38);
  err += this->set_var("offset_x", 1.38);
  err += this->set_var("offset_y", 1.38);

  err += this->set_var("kTx", 1.38);
  err += this->set_var("kTy", 1.38);
  err += this->set_var("offset_Tx", 1.38);
  err += this->set_var("offset_Ty", 1.38);

  return err;
}


template <typename Scalar>
Scalar MASA::periodic_argon_ternary_2d<Scalar>::eval_q_state(Scalar x,Scalar y,int eq)
{
  using std::cos;
  using std::sin;

  Scalar Q_u_t = y;

  return(Q_u_t);
}

/* ------------------------------------------------
 *
 *
 *   Analytical terms
 *
 * -----------------------------------------------
 */

template <typename Scalar>
Scalar MASA::periodic_argon_ternary_2d<Scalar>::eval_exact_state(Scalar x,Scalar y,int eq)
{
  using std::cos;
  using std::sin;

  Scalar exact_u = u0 + dux * sin(2.0 * pi * kux * (x / Lx - offset_ux))
                      + duy * sin(2.0 * pi * kuy * (y / Ly - offset_uy));
  Scalar exact_v = v0 + dvx * sin(2.0 * pi * kvx * (x / Lx - offset_vx))
                      + dvy * sin(2.0 * pi * kvy * (y / Ly - offset_vy));
  Scalar exact_ne = n0 * (X0 + dX * sin(2.0 * pi * kx * (x / Lx - offset_x)) * sin(2.0 * pi * ky * (y / Ly - offset_y)));
  Scalar exact_T = T0 + dT * sin(2.0 * pi * kTx * (x / Lx - offset_Tx)) * sin(2.0 * pi * kTy * (y / Ly - offset_Ty));
  Scalar exact_rho = (mI + mE) * exact_ne + mA * (n0 - 2.0 * exact_ne);

  Scalar heatCapacity = (CV_I + CV_E) * exact_ne + CV_A * (n0 - 2.0 * exact_ne);
  Scalar exact_rhoE = 0.5 * exact_rho * (exact_u * exact_u + exact_v * exact_v)
                      + heatCapacity * exact_T + exact_ne * formEnergy_I;

  switch (eq) {
    case 0: // rho
      return exact_rho;
    break;
    case 1: // rho * u
      return exact_rho * exact_u;
    break;
    case 2: // rho * v
      return exact_rho * exact_v;
    break;
    case 3: // rho * E
      return exact_rhoE;
    break;
    case 4: // rho * Y_I
      return mI * exact_ne;
    break;
    default:
      return -1.0;
    break;
  }
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::periodic_argon_ternary_2d);
