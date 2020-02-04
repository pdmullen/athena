#ifndef PARTICLES_PARTICLE_GRAVITY_HPP_
#define PARTICLES_PARTICLE_GRAVITY_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle_gravity.hpp
//  \brief provides the header for the ParticleGravity class.

// Forward definition
class Particles;

//--------------------------------------------------------------------------------------
//! \class ParticleGravity
//  \brief defines the class for managing the gravity on particles.

class ParticleGravity {
 public:
   // Constructor
   ParticleGravity(Particles *ppar);

   // Destructor
   ~ParticleGravity();

 private:
   // Attributes
   AthenaArray<Real> gforce;  // gravitational force
   Particles *pmy_par;  // pointer to parent Particles instance
   int nx1, nx2, nx3;   // block dimensions
};

#endif  // PARTICLES_PARTICLE_GRAVITY_HPP_
