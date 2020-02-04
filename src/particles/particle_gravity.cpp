//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle_gravity.cpp
//  \brief implements the members of the ParticleGravity class.

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "particle_gravity.hpp"
#include "particles.hpp"

//--------------------------------------------------------------------------------------
//! \fn ParticleGravity::ParticleGravity(Particles *ppar)
//  \brief constructs a new ParticleGravity instance.

ParticleGravity::ParticleGravity(Particles *ppar) {
  // Remember my parent Particles instance.
  pmy_par = ppar;

  // Remember the dimensions of my meshblock.
  RegionSize& block_size = ppar->pmy_block->block_size;
  nx1 = block_size.nx1;
  nx2 = block_size.nx2;
  nx3 = block_size.nx3;

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
