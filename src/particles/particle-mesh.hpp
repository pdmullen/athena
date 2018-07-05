//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-mesh.hpp
//  \brief defines ParticleMesh class used for communication between meshblocks needed 
//         by particle-mesh methods.

// Athena++ classes headers
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../mesh/mesh.hpp"

// Particle-mesh constants.
const Real RINF = 1;  // radius of influence
const int NGPM = 1;   // number of ghost cells needed.

//--------------------------------------------------------------------------------------
//! \class ParticleMesh
//  \brief defines the class for particle-mesh methods

class ParticleMesh {

public:
  // Constructor and destructor
  ParticleMesh(int nmeshaux, MeshBlock *pmb);
  ~ParticleMesh();

private:
  // Instance Variables
  AthenaArray<Real> meshaux;         // auxiliaries to the meshblock
  int nmeshaux_;                     // number of auxiliaries to the meshblock
  int is_, ie_, js_, je_, ks_, ke_;  // beginning and ending indices

  MeshBlock *pmb_;         // ptr to my meshblock
  BoundaryValues *pbval_;  // ptr to my BoundaryValues
  BoundaryData bd_;        // boundary data

};