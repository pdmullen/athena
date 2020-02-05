//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle_gravity.cpp
//  \brief implements the members of the ParticleGravity class.

// C++ standard libraries
#include <cstring>  // strcmp()
#include <sstream>  // stringstream

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "particle_gravity.hpp"
#include "particles.hpp"

// Class variables
int ParticleGravity::iwx(-1), ParticleGravity::iwy(-1), ParticleGravity::iwz(-1);

//--------------------------------------------------------------------------------------
//! \fn ParticleGravity::ParticleGravity(Particles *ppar)
//  \brief constructs a new ParticleGravity instance.

ParticleGravity::ParticleGravity(Particles *ppar) {
  // Remember my parent Particles instance.
  pmy_par = ppar;

  // Remember the coordinates.
  pcoord = ppar->pmy_block->pcoord;

  // Remember the dimensions of my meshblock.
  RegionSize& block_size = ppar->pmy_block->block_size;
  nx1 = block_size.nx1;
  nx2 = block_size.nx2;
  nx3 = block_size.nx3;
  active1 = nx1 > 1;
  active2 = nx2 > 1;
  active3 = nx3 > 1;

  // Allocate space for gravitational force.
  gforce.NewAthenaArray(3, nx3, nx2, nx1);
}

//--------------------------------------------------------------------------------------
//! \fn ParticleGravity::~ParticleGravity(Particles *ppar)
//  \brief constructs a new ParticleGravity instance.

ParticleGravity::~ParticleGravity() {
  // Deallocate space for gravitational force.
  gforce.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleGravity::FindGravitationalForce(const AthenaArray<Real>& phi)
//  \brief computes the gravitational force from the potential phi.
// TODO(ccyang): currently only works with uniform cartesian.

void ParticleGravity::FindGravitationalForce(const AthenaArray<Real>& phi) {
  // Sanity check.
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleGravity::FindGravitationalForce]"
        << std::endl
        << "not implemented for non-Cartesian coordinates. " << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  if ((active1 && pcoord->dx1v(0) != pcoord->dx1v(1)) ||
      (active2 && pcoord->dx2v(0) != pcoord->dx2v(1)) ||
      (active3 && pcoord->dx3v(0) != pcoord->dx3v(1))) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleGravity::FindGravitationalForce]"
        << std::endl
        << "not implemented for non-uniform grid. " << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Set up the loop dimensions.
  const int is = active1 ? 1 : 0;
  const int ie = active1 ? nx1 - 2 : 0;
  const int js = active2 ? 1 : 0;
  const int je = active2 ? nx2 - 2 : 0;
  const int ks = active3 ? 1 : 0;
  const int ke = active3 ? nx3 - 2 : 0;

  // Get the grid spacing.
  Real a1(0.5 / pcoord->dx1v(0)), a2(0.5 / pcoord->dx2v(0)), a3(0.5 / pcoord->dx3v(0));

  // Compute the force.
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        gforce(0,k,j,i) = active1 ? a1 * (phi(k,j,i-1) - phi(k,j,i+1)) : 0.0;
        gforce(1,k,j,i) = active2 ? a2 * (phi(k,j-1,i) - phi(k,j+1,i)) : 0.0;
        gforce(2,k,j,i) = active3 ? a3 * (phi(k-1,j,i) - phi(k+1,j,i)) : 0.0;
      }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleGravity::Initialize()
//  \brief initializes the class.

void ParticleGravity::Initialize() {
  // Get the indices to working arrays for particles.
  iwx = Particles::AddWorkingArray();
  iwy = Particles::AddWorkingArray();
  iwz = Particles::AddWorkingArray();
}
