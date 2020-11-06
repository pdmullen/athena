//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particle_buffer.cpp
//! \brief implements ParticleBuffer class for communication of particles.
//======================================================================================

// C++ standard libraries
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "particle_buffer.hpp"

// Class variable initialization
int ParticleBuffer::nint = 0, ParticleBuffer::nreal = 0;

//--------------------------------------------------------------------------------------
//! \fn void ParticleBuffer::SetNumberOfProperties(int nint0, int nreal0)
//! \brief sets the number of integer and real propertes.

void ParticleBuffer::SetNumberOfProperties(int nint0, int nreal0) {
  // Sanity check
  if (nint0 < 0 || nreal0 < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::SetNumberOfProperties]"
        << std::endl
        << "Invalid nint0 = " << nint0 << " or nreal0 = " << nreal0 << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  // Set the numbers.
  nint = nint0;
  nreal = nreal0;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleBuffer::ParticleBuffer()
//! \brief initiates a default instance of ParticleBuffer.

ParticleBuffer::ParticleBuffer() {
  ibuf = NULL;
  rbuf = NULL;
  nparmax = npar = 0;
#ifdef MPI_PARALLEL
  reqi = reqr = MPI_REQUEST_NULL;
  flagn = flagi = flagr = 0;
  tag = -1;
#endif
}

//--------------------------------------------------------------------------------------
//! \fn ParticleBuffer::ParticleBuffer(int nparmax0)
//! \brief initiates a new instance of ParticleBuffer with nparmax = nparmax0.

ParticleBuffer::ParticleBuffer(int nparmax0) {
  // Sanity check
  if (nparmax0 <= 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::ParticleBuffer]" << std::endl
        << "Invalid nparmax0 = " << nparmax0 << std::endl;
    ATHENA_ERROR(msg);

    ibuf = NULL;
    rbuf = NULL;
    nparmax = npar = 0;
    return;
  }

  // Initialize the instance variables.
  nparmax = nparmax0;
  ibuf = new int[nint * nparmax];
  rbuf = new Real[nreal * nparmax];
  npar = 0;
#ifdef MPI_PARALLEL
  reqi = reqr = MPI_REQUEST_NULL;
  flagn = flagi = flagr = 0;
  tag = -1;
#endif
}

//--------------------------------------------------------------------------------------
//! \fn ParticleBuffer::ParticleBuffer(int nparmax0)
//! \brief destroys an instance of ParticleBuffer.

ParticleBuffer::~ParticleBuffer() {
  if (ibuf != NULL) delete [] ibuf;
  if (rbuf != NULL) delete [] rbuf;
#ifdef MPI_PARALLEL
  if (reqi != MPI_REQUEST_NULL) MPI_Request_free(&reqi);
  if (reqr != MPI_REQUEST_NULL) MPI_Request_free(&reqr);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleBuffer::Reallocate(int new_nparmax)
//! \brief reallocates the buffers; the old content is preserved.

void ParticleBuffer::Reallocate(int new_nparmax) {
  // Sanity check
  if (new_nparmax <= 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::Reallocate]" << std::endl
        << "Invalid new_nparmax = " << new_nparmax << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
  if (new_nparmax < npar) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::Reallocate]" << std::endl
        << "new_nparmax = " << new_nparmax << " < npar = " << npar << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
#ifdef MPI_PARALLEL
  if (reqi != MPI_REQUEST_NULL || reqr != MPI_REQUEST_NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::Reallocate]" << std::endl
        << "MPI requests are active. " << std::endl;
    ATHENA_ERROR(msg);
    return;
  }
#endif

  // Allocate new space.
  nparmax = new_nparmax;
  int *ibuf_new = new int[nint * nparmax];
  Real *rbuf_new = new Real[nreal * nparmax];

  // Move existing data.
  if (npar > 0) {
    std::memcpy(ibuf_new, ibuf, nint * npar * sizeof(int));
    std::memcpy(rbuf_new, rbuf, nreal * npar * sizeof(Real));
  }

  // Delete old space.
  if (ibuf != NULL) delete [] ibuf;
  if (rbuf != NULL) delete [] rbuf;
  ibuf = ibuf_new;
  rbuf = rbuf_new;
}
