//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grav_src_task_list.cpp
//! \brief Self-Gravity Energy Source Term
//!
//! REFERENCE:
//! Mullen, P. D., Hanawa, T., & Gammie, C. F. 2020, ApJS

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//  GravitySourceTaskList constructor

GravitySourceTaskList::GravitySourceTaskList(
  ParameterInput *pin, Mesh *pm, TimeIntegratorTaskList *ptlist) :
  ptlist_(ptlist) {
  nstages = ptlist_->nstages;

  if (ptlist_->integrator != "vl2" && ptlist_->integrator != "rk2"
      && ptlist_->integrator != "rk3") {
    std::stringstream msg;
    msg << "### FATAL ERROR in GravitySourceTaskList constructor" << std::endl
        << "integrator=" << ptlist_->integrator << " not yet compatible with "
        << "non-isothermal self-gravitating hydrodynamics" << std::endl;
    ATHENA_ERROR(msg);
  }

  // Applies self-gravity source terms
  // Now assemble list of tasks for each stage of time integrator
  {using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
    // calculate hydro/field diffusive fluxes
    AddTask(SRCTERM_HYD,NONE);
    AddTask(SEND_HYD,SRCTERM_HYD);
    AddTask(RECV_HYD,NONE);
    AddTask(SETB_HYD,(RECV_HYD|SRCTERM_HYD));
    AddTask(CONS2PRIM,SETB_HYD);
    AddTask(PHY_BVAL,CONS2PRIM);
    if (!STS_ENABLED || pm->sts_integrator == "rkl1") {
      AddTask(USERWORK,PHY_BVAL);
      AddTask(NEW_DT,USERWORK);
      if (pm->adaptive) {
        AddTask(FLAG_AMR,USERWORK);
        AddTask(CLEAR_ALLBND,FLAG_AMR);
      } else {
        AddTask(CLEAR_ALLBND,NEW_DT);
      }
    } else {
      AddTask(CLEAR_ALLBND,PHY_BVAL);
    }
  } // end of using namespace block
}

//---------------------------------------------------------------------------------------
//  Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//  ntask.

void GravitySourceTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;

  using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)

  if (id == CLEAR_ALLBND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySourceTaskList::ClearAllBoundary_GSRC);
    task_list_[ntasks].lb_time = false;
  } else if (id == SRCTERM_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&GravitySourceTaskList::AddSourceTermsHydro_GSRC);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendField);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydro);
    task_list_[ntasks].lb_time = false;
  } else if (id == RECV_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveField);
    task_list_[ntasks].lb_time = false;
  } else if (id == SETB_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SETB_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesField);
    task_list_[ntasks].lb_time = true;
  } else if (id == CONS2PRIM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
    task_list_[ntasks].lb_time = true;
  } else if (id == PHY_BVAL) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
    task_list_[ntasks].lb_time = true;
  } else if (id == USERWORK) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
    task_list_[ntasks].lb_time = true;
  } else if (id == NEW_DT) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
    task_list_[ntasks].lb_time = true;
  } else if (id == FLAG_AMR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
    task_list_[ntasks].lb_time = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
// Functions to end MPI communication

TaskStatus GravitySourceTaskList::ClearAllBoundary_GSRC(MeshBlock *pmb, int stage) {
  pmb->pbval->ClearBoundarySubset(BoundaryCommSubset::all,
                                  pmb->pbval->bvars_gsrc);
  return TaskStatus::success;
}


void GravitySourceTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  pmb->pbval->StartReceivingSubset(BoundaryCommSubset::all, pmb->pbval->bvars_gsrc);
  return;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

TaskStatus GravitySourceTaskList::AddSourceTermsHydro_GSRC(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  if (stage <= ptlist_->nstages) {
    // Evaluate the self-gravity energy source term
    Real wght[2];
    wght[0] = 0.0;
    wght[1] = ptlist_->stage_wghts[(stage-1)].beta;
    if (ptlist_->integrator != "vl2") {
      if (stage==2) {
        wght[0] = wght[1];
        wght[1] = 2.0 * wght[0];
      } else if (stage==3) {
        wght[0] = 2.0 * ptlist_->stage_wghts[(stage-2)].beta * wght[1];
        wght[1] = wght[0] + ptlist_->stage_wghts[(stage-1)].beta;
      }
    }
    ph->hsrc.SelfGravityEnergy(pmb->pmy_mesh->dt, stage, wght,
                               ph->flux, ph->fl0, ph->fl1,
                               ph->u);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::next;
}
