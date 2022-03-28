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
MASA::test_plasma_2d<Scalar>::test_plasma_2d()
{
  this->mmsname = "test_plasma_2d";
  this->dimension=2;
  this->register_var("test_var",&test_var);

  // init defaults
  this->init_var();

}

template <typename Scalar>
int MASA::test_plasma_2d<Scalar>::init_var()
{
  int err = 0;

  // randomly generated
  err += this->set_var("test_var",1.38);

  return err;
}


template <typename Scalar>
Scalar MASA::test_plasma_2d<Scalar>::eval_q_state(Scalar x,Scalar y,int eq)
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
Scalar MASA::test_plasma_2d<Scalar>::eval_exact_state(Scalar x,Scalar y,int eq)
{
  using std::cos;
  using std::sin;

  Scalar exact_u = x;

  return exact_u;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::test_plasma_2d);
