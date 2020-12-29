//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to self-gravity

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../gravity/gravity.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::SelfGravityMomentum
//  \brief Adds source terms for self-gravitational acceleration to conserved variables

void HydroSourceTerms::SelfGravityMomentum(const Real dt,const AthenaArray<Real> *flux,
                                   const AthenaArray<Real> &prim,
                                   AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  if (SELF_GRAVITY_ENABLED && !GRAVITY_FLUX_ENABLED) {
    Gravity *pgrav = pmb->pgrav;
    //Real four_pi_G = pgrav->four_pi_G;
    //Real grav_mean_rho = pgrav->grav_mean_rho;
    // acceleration in 1-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx1 = pmb->pcoord->dx1v(i);
          //Real dx2 = pmb->pcoord->dx2v(j);
          //Real dx3 = pmb->pcoord->dx3v(k);
          Real dtodx1 = dt/dx1;
          Real phil = 0.5*(pgrav->phi(k,j,i-1)+pgrav->phi(k,j,i  ));
          Real phir = 0.5*(pgrav->phi(k,j,i  )+pgrav->phi(k,j,i+1));
          // Update momenta and energy with d/dx1 terms
          cons(IM1,k,j,i) -= dtodx1*prim(IDN,k,j,i)*(phir-phil);
        }
      }
    }

    if (pmb->block_size.nx2 > 1) {
      // acceleration in 2-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            //Real dx1 = pmb->pcoord->dx1v(i);
            Real dx2 = pmb->pcoord->dx2v(j);
            //Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx2 = dt/dx2;
            Real phil = 0.5*(pgrav->phi(k,j-1,i)+pgrav->phi(k,j  ,i));
            Real phir = 0.5*(pgrav->phi(k,j  ,i)+pgrav->phi(k,j+1,i));
            cons(IM2,k,j,i) -= dtodx2*prim(IDN,k,j,i)*(phir-phil);
          }
        }
      }
    }

    if (pmb->block_size.nx3 > 1) {
      // acceleration in 3-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            //Real dx1 = pmb->pcoord->dx1v(i);
            //Real dx2 = pmb->pcoord->dx2v(j);
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx3 = dt/dx3;
            Real phil = 0.5*(pgrav->phi(k-1,j,i)+pgrav->phi(k  ,j,i));
            Real phir = 0.5*(pgrav->phi(k  ,j,i)+pgrav->phi(k+1,j,i));
            cons(IM3,k,j,i) -= dtodx3*prim(IDN,k,j,i)*(phir-phil);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::SelfGravityEnergy
//  \brief Adds source terms for self-gravitational energy releases to conserved variables

void HydroSourceTerms::SelfGravityEnergy(const Real dt,Real stage,const Real wght[2],
                                   const AthenaArray<Real> *fl,
                                   const AthenaArray<Real> *fl0,
                                   const AthenaArray<Real> *fl1,
                                   AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  if (SELF_GRAVITY_ENABLED && NON_BAROTROPIC_EOS) {
    Gravity *pgrav = pmb->pgrav;
    //Real four_pi_G = pgrav->four_pi_G;
    //Real grav_mean_rho = pgrav->grav_mean_rho;
    // acceleration in 1-direction
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real dx1 = pmb->pcoord->dx1v(i);
          //Real dx2 = pmb->pcoord->dx2v(j);
          //Real dx3 = pmb->pcoord->dx3v(k);
          Real dtodx1 = wght[1]*dt/dx1;
          Real phic = 0.5*  (pgrav->phi (k,j,i  )+pgrav->phi0(k,j,i  ));
          Real phil = 0.25*((pgrav->phi (k,j,i-1)+pgrav->phi (k,j,i  ))
                      +(pgrav->phi0(k,j,i-1)+pgrav->phi0(k,j,i  )));
          Real phir = 0.25*((pgrav->phi (k,j,i  )+pgrav->phi (k,j,i+1))
                      +(pgrav->phi0(k,j,i  )+pgrav->phi0(k,j,i+1)));
          Real fluxl = fl[X1DIR](IDN,k,j,i  );
          Real fluxr = fl[X1DIR](IDN,k,j,i+1);
          if (pmb->pmy_mesh->integrator != "vl2") {
            if (stage==2) {
              fluxl = (fluxl + fl0[X1DIR](k,j,i  ))/2.0;
              fluxr = (fluxr + fl0[X1DIR](k,j,i+1))/2.0;
            } else if (stage==3) {
              fluxl = (4.0*fluxl+fl0[X1DIR](k,j,i  )+fl1[X1DIR](k,j,i  ))/6.0;
              fluxr = (4.0*fluxr+fl0[X1DIR](k,j,i+1)+fl1[X1DIR](k,j,i+1))/6.0;
            }
          }
          // Update momenta and energy with d/dx1 terms
          cons(IEN,k,j,i) -= dtodx1*(fluxl*(phic - phil) +
                                     fluxr*(phir - phic));
          // add back contributions from previous stage (if necessary)
          if (pmb->pmy_mesh->integrator != "vl2" && stage >= 2) {
            dtodx1 = wght[0]*dt/dx1;
            phic = 0.5 *( pgrav->phi1(k,j,i  )+pgrav->phi0(k,j,i  ));
            phil = 0.25*((pgrav->phi1(k,j,i-1)+pgrav->phi1(k,j,i  ))
                        +(pgrav->phi0(k,j,i-1)+pgrav->phi0(k,j,i  )));
            phir = 0.25*((pgrav->phi1(k,j,i  )+pgrav->phi1(k,j,i+1))
                        +(pgrav->phi0(k,j,i  )+pgrav->phi0(k,j,i+1)));
            if (stage==2) {
              fluxl = fl0[X1DIR](k,j,i  );
              fluxr = fl0[X1DIR](k,j,i+1);
              cons(IEN,k,j,i) += dtodx1*(fluxl*(phic - phil) +
                                         fluxr*(phir - phic));
            } else if (stage==3) {
              fluxl = 0.5*(fl0[X1DIR](k,j,i  ) + fl1[X1DIR](k,j,i  ));
              fluxr = 0.5*(fl0[X1DIR](k,j,i+1) + fl1[X1DIR](k,j,i+1));
              cons(IEN,k,j,i) += dtodx1*(fluxl*(phic - phil) +
                                         fluxr*(phir - phic));
            }
          }
        }
      }
    }

    if (pmb->block_size.nx2 > 1) {
      // acceleration in 2-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            //Real dx1 = pmb->pcoord->dx1v(i);
            Real dx2 = pmb->pcoord->dx2v(j);
            //Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx2 = wght[1]*dt/dx2;
            Real phic = 0.5*  (pgrav->phi (k,j  ,i)+pgrav->phi0(k,j  ,i));
            Real phil = 0.25*((pgrav->phi (k,j-1,i)+pgrav->phi (k,j  ,i))
                        +(pgrav->phi0(k,j-1,i)+pgrav->phi0(k,j  ,i)));
            Real phir = 0.25*((pgrav->phi (k,j  ,i)+pgrav->phi (k,j+1,i))
                        +(pgrav->phi0(k,j  ,i)+pgrav->phi0(k,j+1,i)));
            Real fluxl = fl[X2DIR](IDN,k,j  ,i);
            Real fluxr = fl[X2DIR](IDN,k,j+1,i);
            if (pmb->pmy_mesh->integrator != "vl2") {
              if (stage==2) {
                fluxl = (fluxl + fl0[X2DIR](k,j  ,i))/2.0;
                fluxr = (fluxr + fl0[X2DIR](k,j+1,i))/2.0;
              } else if (stage==3) {
                fluxl = (4.0*fluxl+fl0[X2DIR](k,j  ,i)+fl1[X2DIR](k,j  ,i))/6.0;
                fluxr = (4.0*fluxr+fl0[X2DIR](k,j+1,i)+fl1[X2DIR](k,j+1,i))/6.0;
              }
            }
            // Update momenta and energy with d/dx1 terms
            cons(IEN,k,j,i) -= dtodx2*(fluxl*(phic - phil) +
                                       fluxr*(phir - phic));
            // add back contributions from previous stage (if necessary)
            if (pmb->pmy_mesh->integrator != "vl2" && stage > 1) {
              dtodx2 = wght[0]*dt/dx2;
              phic = 0.5 *( pgrav->phi1(k,j  ,i)+pgrav->phi0(k,j  ,i));
              phil = 0.25*((pgrav->phi1(k,j-1,i)+pgrav->phi1(k,j  ,i))
                          +(pgrav->phi0(k,j-1,i)+pgrav->phi0(k,j  ,i)));
              phir = 0.25*((pgrav->phi1(k,j  ,i)+pgrav->phi1(k,j+1,i))
                          +(pgrav->phi0(k,j  ,i)+pgrav->phi0(k,j+1,i)));
              if (stage==2) {
                fluxl = fl0[X2DIR](k,j  ,i);
                fluxr = fl0[X2DIR](k,j+1,i);
              } else { // stage==3
                fluxl = 0.5*(fl0[X2DIR](k,j  ,i)+fl1[X2DIR](k,j  ,i));
                fluxr = 0.5*(fl0[X2DIR](k,j+1,i)+fl1[X2DIR](k,j+1,i));
              }
              cons(IEN,k,j,i) += dtodx2*(fluxl*(phic - phil) +
                                         fluxr*(phir - phic));
            }
          }
        }
      }
    }

    if (pmb->block_size.nx3 > 1) {
      // acceleration in 3-direction
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            //Real dx1 = pmb->pcoord->dx1v(i);
            //Real dx2 = pmb->pcoord->dx2v(j);
            Real dx3 = pmb->pcoord->dx3v(k);
            Real dtodx3 = wght[1]*dt/dx3;
            Real phic = 0.5*  (pgrav->phi (k  ,j,i)+pgrav->phi0(k  ,j,i));
            Real phil = 0.25*((pgrav->phi (k-1,j,i)+pgrav->phi (k  ,j,i))
                        +(pgrav->phi0(k-1,j,i)+pgrav->phi0(k  ,j,i)));
            Real phir = 0.25*((pgrav->phi (k  ,j,i)+pgrav->phi (k+1,j,i))
                        +(pgrav->phi0(k  ,j,i)+pgrav->phi0(k+1,j,i)));
            Real fluxl = fl[X3DIR](IDN,k  ,j,i);
            Real fluxr = fl[X3DIR](IDN,k+1,j,i);
            if (pmb->pmy_mesh->integrator != "vl2") {
              if (stage==2) {
                fluxl = (fluxl + fl0[X3DIR](k  ,j,i))/2.0;
                fluxr = (fluxr + fl0[X3DIR](k+1,j,i))/2.0;
              } else if (stage==3) {
                fluxl = (4.0*fluxl+fl0[X3DIR](k  ,j,i)+fl1[X3DIR](k  ,j,i))/6.0;
                fluxr = (4.0*fluxr+fl0[X3DIR](k+1,j,i)+fl1[X3DIR](k+1,j,i))/6.0;
              }
            }
            // Update momenta and energy with d/dx1 terms
            cons(IEN,k,j,i) -= dtodx3*(fluxl*(phic - phil) +
                                       fluxr*(phir - phic));
            // add back contributions from previous stage (if necessary)
            if (pmb->pmy_mesh->integrator != "vl2" && stage > 1) {
              dtodx3 = wght[0]*dt/dx3;
              phic = 0.5 *( pgrav->phi1(k  ,j,i)+pgrav->phi0(k  ,j,i));
              phil = 0.25*((pgrav->phi1(k-1,j,i)+pgrav->phi1(k  ,j,i))
                          +(pgrav->phi0(k-1,j,i)+pgrav->phi0(k  ,j,i)));
              phir = 0.25*((pgrav->phi1(k  ,j,i)+pgrav->phi1(k+1,j,i))
                          +(pgrav->phi0(k  ,j,i)+pgrav->phi0(k+1,j,i)));
              if (stage==2) {
                fluxl = fl0[X3DIR](k  ,j,i);
                fluxr = fl0[X3DIR](k+1,j,i);
              } else { // stage==3
                fluxl = 0.5*(fl0[X3DIR](k  ,j,i)+fl1[X3DIR](k  ,j,i));
                fluxr = 0.5*(fl0[X3DIR](k+1,j,i)+fl1[X3DIR](k+1,j,i));
              }
              cons(IEN,k,j,i) += dtodx3*(fluxl*(phic - phil) +
                                         fluxr*(phir - phic));
            }
          }
        }
      }
    }
  }
  return;
}
