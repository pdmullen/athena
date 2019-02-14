//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX1(
//                         const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void CellCenteredBoundaryVariable::ReflectInnerX1(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<=nu_; ++n) {
    if (n==(IVX)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            var_cc(IVX,k,j,il-i) = -var_cc(IVX,k,j,(il+i-1));  // reflect 1-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            var_cc(n,k,j,il-i) = var_cc(n,k,j,(il+i-1));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX1(
//                         const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void CellCenteredBoundaryVariable::ReflectOuterX1(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<=nu_; ++n) {
    if (n==(IVX)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            var_cc(IVX,k,j,iu+i) = -var_cc(IVX,k,j,(iu-i+1));  // reflect 1-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            var_cc(n,k,j,iu+i) = var_cc(n,k,j,(iu-i+1));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflecInnerX2(
//                         const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void CellCenteredBoundaryVariable::ReflectInnerX2(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<=nu_; ++n) {
    if (n==(IVY)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(IVY,k,jl-j,i) = -var_cc(IVY,k,jl+j-1,i);  // reflect 2-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(n,k,jl-j,i) = var_cc(n,k,jl+j-1,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX2(
//                         const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void CellCenteredBoundaryVariable::ReflectOuterX2(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<=nu_; ++n) {
    if (n==(IVY)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(IVY,k,ju+j,i) = -var_cc(IVY,k,ju-j+1,i);  // reflect 2-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(n,k,ju+j,i) = var_cc(n,k,ju-j+1,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectInnerX3(
//                         const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void CellCenteredBoundaryVariable::ReflectInnerX3(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<=nu_; ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(IVZ,kl-k,j,i) = -var_cc(IVZ,kl+k-1,j,i);  // reflect 3-velocity
          }
        }
      }
    } else {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(n,kl-k,j,i) = var_cc(n,kl+k-1,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CellCenteredBoundaryVariable::ReflectOuterX3(
//                         const Real time, const Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void CellCenteredBoundaryVariable::ReflectOuterX3(
    MeshBlock *pmb, Coordinates *pco, Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<=nu_; ++n) {
    if (n==(IVZ)) {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(IVZ,ku+k,j,i) = -var_cc(IVZ,ku-k+1,j,i);  // reflect 3-velocity
          }
        }
      }
    } else {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            var_cc(n,ku+k,j,i) = var_cc(n,ku-k+1,j,i);
          }
        }
      }
    }
  }
  return;
}
