//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file kida.cpp
//  \brief implementation of functions in class ChemNetwork, using the simple
//  network for kida style network files.
//======================================================================================

// this class header
#include "kida.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../../radiation/radiation.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../../utils/string_utils.hpp"
#include "../utils/thermo.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf
#include <fstream>   //file()
#include <stdio.h>    // c style file
#include <algorithm>    // std::find()
#include <iterator>     // std::distance()

#ifdef DEBUG
static bool output_rates = true;
static bool output_thermo = true;
#endif

ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) :  
  pmy_spec_(pmb->pscalars), pmy_mb_(pmb), id7max_(0), n_cr_(0),
  icr_H_(-1), icr_H2_(-1), icr_He_(-1), n_crp_(0), n_ph_(0), iph_H2_(-1), 
  n_2body_(0), i2body_H2_H_(-1), i2body_H2_H2_(-1), i2body_H_e_(-1),
  n_2bodytr_(0), n_gr_(0), igr_H_(-1), n_sr_(0),
  n_freq_(0), index_gpe_(0), index_cr_(0), gradv_(0.) {

	//set the parameters from input file
	zdg_ = pin->GetOrAddReal("chemistry", "Zdg", 1.);//dust and gas metallicity
  o2pH2_ = pin->GetOrAddReal("chemistry", "o2pH2", 3.);//ortho to para H2 ratio
  //units
	unit_density_in_nH_ = pin->GetReal("chemistry", "unit_density_in_nH");
	unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
	unit_vel_in_cms_ = pin->GetReal("chemistry", "unit_vel_in_cms");
	unit_radiation_in_draine1987_ = pin->GetReal(
                                "chemistry", "unit_radiation_in_draine1987");
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
  //temperature
  if (NON_BAROTROPIC_EOS) {
    temperature_ = 0.;
  } else {
    //isothermal
    temperature_ = pin->GetReal("chemistry", "temperature");
  }
  //whether to cap temperature if the reaction is outside of the temperature range
  //only for 2 body reactions
  is_Tcap_2body_ = pin->GetOrAddBoolean("chemistry", "is_Tcap_2body", false);
	//minimum temperature for reaction rates, also applied to energy equation
	temp_min_rates_ = pin->GetOrAddReal("chemistry", "temp_min_rates", 1.);
  //minimum temperature below which cooling is turned off
	temp_min_cool_ = pin->GetOrAddReal("chemistry", "temp_min_cool", 1.);
  //dust temperature for dust thermo cooling of the gas at high densities
  temp_dust_thermo_ = pin->GetOrAddReal("chemistry", "temp_dust_thermo", 10.);
  //folder of the network
  network_dir_ = pin->GetString("chemistry", "network_dir");
	//CO cooling parameters
	//Maximum CO cooling length in cm. default 100pc.
	Leff_CO_max_ = pin->GetOrAddReal("chemistry", "Leff_CO_max", 3.0e20);

  //read in the species
  std::string species_file_name = network_dir_ + "/species.dat";
  std::ifstream species_file(species_file_name);
  if (!species_file) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "Cannot open file" << species_file_name << std::endl;
    ATHENA_ERROR(msg);
  }
  std::string line;
  int nline = 0;
  while (getline(species_file, line)) {
    //trim white spaces
    StringUtils::trim(line);
    //skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
        continue;
    }
    KidaSpecies si(line, nline);
    species_.push_back(si);
    species_names[nline] = si.name;
    ispec_map_[si.name] = nline;
    nline++;
  }
  if (nline != NSCALARS) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "number of species in species.dat does not match the number of scalars (" 
      << NSCALARS << ")" << std::endl;
    ATHENA_ERROR(msg);
  }

  //read in the reactions
  std::string reactions_file_name = network_dir_ + "/reactions.dat";
  std::ifstream reactions_file(reactions_file_name);
  std::vector<int> rids; //array of id for reactions
  if (!reactions_file) {
    std::stringstream msg; //error message
    msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
      << "Cannot open file" << reactions_file_name << std::endl;
    ATHENA_ERROR(msg);
  }
  nline = 0;
  while (getline(reactions_file, line)) {
    //trim white spaces
    StringUtils::trim(line);
    //skip blank lines and comments
    if(line.empty() || (line.find("!") == 0)) {
        continue;
    }
    KidaReaction ri(line);
    if (ri.id_ <= 0) {
      std::stringstream msg; //error message
      msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]" << std::endl
        << "reaction ID ( " << ri.id_ << ") not positive" << std::endl;
      ATHENA_ERROR(msg);
    }
    auto itr = std::find(rids.rbegin(), rids.rend(), ri.id_); 
    const int ifind = rids.rend() - itr - 1;
    if (ifind < 0) {
      rids.push_back(ri.id_);
      if (ri.id_ > id7max_ && ri.formula_ == 7) {
        id7max_ = ri.id_;
      }
      reactions_.push_back(ri);
    } else {
      if (reactions_[ifind].reactants_ == ri.reactants_ 
          && reactions_[ifind].products_ == ri.products_
          && reactions_[ifind].itype_ == ri.itype_ 
          && reactions_[ifind].Tmax_ < ri.Tmin_
          && reactions_[ifind].formula_ != 7 && ri.formula_ != 7
          && ri.itype_ >= 4 && ri.itype_<= 8 ) {
        if (ifind != nline - 1) {
          std::stringstream msg; //error message
          msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]"
            << std::endl << "reactions ID ( " << ri.id_ << " ) of different" 
            << " temperature ranges are not arranged next to each other."
            << std::endl;
          ATHENA_ERROR(msg);
        }
        rids.push_back(ri.id_);
        reactions_.push_back(ri);
        const int nprev = std::count(id_2bodytr_.begin(), id_2bodytr_.end(), ri.id_);
        if (nprev == 0) {
          id_2bodytr_.push_back(ri.id_);
          n_2bodytr_++;
        } else if (nprev < n_range_ - 1) {
          id_2bodytr_.push_back(ri.id_);
        } else {
          std::stringstream msg; //error message
          msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]"
            << std::endl << "reaction ID ( " << ri.id_ << " )" 
            << "too many temperature ranges" << std::endl;
          ATHENA_ERROR(msg);
        }
      } else {
        std::stringstream msg; //error message
        msg << "### FATAL ERROR in ChemNetwork constructor [ChemNetwork]"
          << std::endl << "reaction ID ( " << ri.id_ << ") not unique" << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    nline++;
  }
  nr_ = nline; //number of reactions
  if (id7max_ > 0) {
    id7map_.NewAthenaArray(id7max_+1);
    id7type_.NewAthenaArray(id7max_+1);
    for (int i=0; i<id7max_+1; i++) {
      id7map_(i) = -1;
      id7type_(i) = ReactionType::none;
    }
  }

  //initialize coefficients of reactions
  InitializeReactions();

  //radiation related variables
  const int nfreq = pin->GetOrAddInteger("radiation", "n_frequency", 1);
  n_freq_ = n_ph_ + 2;
  std::stringstream msg;
  //check whether number of frequencies equal to the input file specification
  if (nfreq != n_freq_) {
    msg << "### FATAL ERROR in ChemNetwork constructor" << std::endl
      << "number of frequencies in radiation: " << nfreq 
      << " not equal to that in chemistry: " << n_freq_  << std::endl;
    ATHENA_ERROR(msg);
  }
  index_gpe_ = n_ph_;
  index_cr_ = n_ph_ + 1;
  rad_.NewAthenaArray(n_freq_);

#ifdef DEBUG
  PrintProperties();
#endif
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::InitializeReactions() {
  KidaReaction *pr = NULL;
  //error message
  bool error=false;
  ReactionType rtype;
  //count reactions
  for (int ir=0; ir<nr_; ir++) {
    CheckReaction(reactions_[ir]);
    pr = &reactions_[ir];
    rtype = SortReaction(pr);
    switch(rtype) {
      case ReactionType::cr: n_cr_++; break;
      case ReactionType::crp: n_crp_++; break;
      case ReactionType::photo: n_ph_++; break;
      case ReactionType::twobody: n_2body_++; break;
      case ReactionType::twobodytr: break;
      case ReactionType::grain: n_gr_++; break;
      case ReactionType::special: n_sr_++; break;
      default: std::stringstream msg; 
               msg << "### FATAL ERROR in ChemNetwork InitializeReactions()"
                 << " [ChemNetwork]: reaction type not recognized." << std::endl;
               ATHENA_ERROR(msg);
               break;
    }
  }
  //create arrays
  if (n_cr_ > 0) {
    incr_.NewAthenaArray(n_cr_);
    outcr1_.NewAthenaArray(n_cr_);
    outcr2_.NewAthenaArray(n_cr_);
    kcr_base_.NewAthenaArray(n_cr_);
    kcr_.NewAthenaArray(n_cr_);
  }
  if (n_crp_ > 0) {
    incrp_.NewAthenaArray(n_crp_);
    outcrp1_.NewAthenaArray(n_crp_);
    outcrp2_.NewAthenaArray(n_crp_);
    kcrp_base_.NewAthenaArray(n_crp_);
    kcrp_.NewAthenaArray(n_crp_);
  }
  if (n_ph_ > 0) {
    inph_.NewAthenaArray(n_ph_);
    outph1_.NewAthenaArray(n_ph_);
    outph2_.NewAthenaArray(n_ph_);
    kph_base_.NewAthenaArray(n_ph_);
    kph_avfac_.NewAthenaArray(n_ph_);
    kph_.NewAthenaArray(n_ph_);
  }
  if (n_2body_ > 0) {
    in2body1_.NewAthenaArray(n_2body_);
    in2body2_.NewAthenaArray(n_2body_);
    out2body1_.NewAthenaArray(n_2body_);
    out2body2_.NewAthenaArray(n_2body_);
    out2body3_.NewAthenaArray(n_2body_);
    out2body4_.NewAthenaArray(n_2body_);
    frml_2body_.NewAthenaArray(n_2body_);
    a2body_.NewAthenaArray(n_2body_);
    b2body_.NewAthenaArray(n_2body_);
    c2body_.NewAthenaArray(n_2body_);
    Tmin_2body_.NewAthenaArray(n_2body_);
    Tmax_2body_.NewAthenaArray(n_2body_);
    k2body_.NewAthenaArray(n_2body_);
  }
  if (n_2bodytr_ > 0) {
    in2bodytr1_.NewAthenaArray(n_2bodytr_);
    in2bodytr2_.NewAthenaArray(n_2bodytr_);
    out2bodytr1_.NewAthenaArray(n_2bodytr_);
    out2bodytr2_.NewAthenaArray(n_2bodytr_);
    out2bodytr3_.NewAthenaArray(n_2bodytr_);
    nr_2bodytr_.NewAthenaArray(n_2bodytr_);
    frml_2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    a2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    b2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    c2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    Tmin_2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    Tmax_2bodytr_.NewAthenaArray(n_2bodytr_, n_range_);
    k2bodytr_.NewAthenaArray(n_2bodytr_);
  }
  if (n_gr_ > 0) {
    ingr1_.NewAthenaArray(n_gr_);
    ingr2_.NewAthenaArray(n_gr_);
    outgr_.NewAthenaArray(n_gr_);
    kgr_.NewAthenaArray(n_gr_);
  }
  if (n_sr_ > 0) {
    insr_.NewAthenaArray(n_sr_, n_insr_);
    outsr_.NewAthenaArray(n_sr_, n_outsr_);
    ksr_.NewAthenaArray(n_sr_);
  }
  
  int icr=0, icrp=0, iph=0, i2body=0, i2bodytr=0, igr=0, isr=0;
  std::vector<int> idtr;
  for (int ir=0; ir<nr_; ir++) {
    pr = &reactions_[ir];
    rtype = SortReaction(pr);
    //---------------- cr - direct cosmic-ray ionization --------------
    if (rtype == ReactionType::cr) {
      std::string in_spec; //input species
      if (pr->reactants_[0] == "CR") {
        in_spec = pr->reactants_[1];
      } else {
        in_spec = pr->reactants_[0];
      }
      if (in_spec == "H") {
        icr_H_ = icr;
      } else if (in_spec == "H2") {
        icr_H2_ = icr;
      } else if (in_spec == "He") {
        icr_He_ = icr;
      }
      if (pr->formula_ == 1) {
        incr_(icr) = ispec_map_[in_spec];
        outcr1_(icr) = ispec_map_[ pr->products_[0]];
        outcr2_(icr) = ispec_map_[ pr->products_[1]];
        kcr_base_(icr) = pr->alpha_;
        kcr_(icr) = 0.;
        icr++;
      } else if (pr->formula_ == 7) {
        incr_(icr) = ispec_map_[in_spec];
        outcr1_(icr) = ispec_map_[ pr->products_[0]];
        outcr2_(icr) = ispec_map_[ pr->products_[1]];
        id7map_(pr->id_) = icr;
        id7type_(pr->id_) = ReactionType::cr;
        kcr_base_(icr) = 0.;
        kcr_(icr) = 0.;
        icr++;
      } else {
        error = true;
      }

    //---------------- crp - cosmic-ray induced photo ionization --------
    } else if (rtype == ReactionType::crp) {
      std::string in_spec; //input species
      if (pr->reactants_[0] == "CRP") {
        in_spec = pr->reactants_[1];
      } else {
        in_spec = pr->reactants_[0];
      }
      if (pr->formula_ == 1) {
        incrp_(icrp) = ispec_map_[in_spec];
        outcrp1_(icrp) = ispec_map_[ pr->products_[0]];
        outcrp2_(icrp) = ispec_map_[ pr->products_[1]];
        kcrp_base_(icrp) = pr->alpha_;
        kcrp_(icrp) = 0.;
        icrp++;
      } else if (pr->formula_ == 7) {
        incrp_(icrp) = ispec_map_[in_spec];
        outcrp1_(icrp) = ispec_map_[ pr->products_[0]];
        outcrp2_(icrp) = ispec_map_[ pr->products_[1]];
        id7map_(pr->id_) = icrp;
        id7type_(pr->id_) = ReactionType::crp;
        kcrp_base_(icrp) = 0.;
        kcrp_(icrp) = 0.;
        icrp++;
      } else {
        error = true;
      }

    //---------------- photo - FUV ionization/dissociation ----------------
    } else if (rtype == ReactionType::photo) {
      std::string in_spec; //input species
      if (pr->reactants_[0] == "Photon") {
        in_spec = pr->reactants_[1];
      } else {
        in_spec = pr->reactants_[0];
      }
      if (in_spec == "H2" && pr->products_[0] == "H" && pr->products_[1] == "H") {
        iph_H2_ = iph;
      }
      if (pr->formula_ == 2) {
        inph_(iph) = ispec_map_[in_spec];
        outph1_(iph) = ispec_map_[ pr->products_[0]];
        outph2_(iph) = ispec_map_[ pr->products_[1]];
        kph_base_(iph) = pr->alpha_;
        kph_avfac_(iph) = pr->gamma_;
        smap_ph_[in_spec] = iph;
        kph_(iph) = 0.;
        iph++;
      } else {
        error = true;
      }

    //---------------- twobody - 2body reaction ---------------------------
    } else if (rtype == ReactionType::twobody) {
      if (pr->reactants_[0] == "H2" && pr->reactants_[1] == "H") {
        i2body_H2_H_ = i2body;
      }
      if (pr->reactants_[0] == "H2" && pr->reactants_[1] == "H2") {
        i2body_H2_H2_ = i2body;
      }
      if (pr->reactants_[0] == "H" && pr->reactants_[1] == "e-") {
        i2body_H_e_ = i2body;
      }
      if (pr->formula_ == 3 || pr->formula_ == 4 || pr->formula_ == 5) {
        in2body1_(i2body) = ispec_map_[ pr->reactants_[0]];
        in2body2_(i2body) = ispec_map_[ pr->reactants_[1]];
        out2body1_(i2body) = ispec_map_[ pr->products_[0]];
        if (pr->products_.size() >= 2) {
          out2body2_(i2body) = ispec_map_[ pr->products_[1]];
        } else {
          out2body2_(i2body) = -1;
        }
        if (pr->products_.size() >= 3) {
          out2body3_(i2body) = ispec_map_[ pr->products_[2]];
        } else {
          out2body3_(i2body) = -1;
        }
        if (pr->products_.size() == 4) {
          out2body4_(i2body) = ispec_map_[ pr->products_[3]];
        } else {
          out2body4_(i2body) = -1;
        }
        frml_2body_(i2body) = pr->formula_;
        a2body_(i2body) = pr->alpha_;
        b2body_(i2body) = pr->beta_;
        c2body_(i2body) = pr->gamma_;
        Tmin_2body_(i2body) = pr->Tmin_;
        Tmax_2body_(i2body) = pr->Tmax_;
        k2body_(i2body) = 0.;
        i2body++;
      } else if (pr->formula_ == 7) {
        in2body1_(i2body) = ispec_map_[ pr->reactants_[0]];
        in2body2_(i2body) = ispec_map_[ pr->reactants_[1]];
        out2body1_(i2body) = ispec_map_[ pr->products_[0]];
        if (pr->products_.size() >= 2) {
          out2body2_(i2body) = ispec_map_[ pr->products_[1]];
        } else {
          out2body2_(i2body) = -1;
        }
        if (pr->products_.size() >= 3) {
          out2body3_(i2body) = ispec_map_[ pr->products_[2]];
        } else {
          out2body3_(i2body) = -1;
        }
        if (pr->products_.size() == 4) {
          out2body4_(i2body) = ispec_map_[ pr->products_[3]];
        } else {
          out2body4_(i2body) = -1;
        }
        frml_2body_(i2body) = pr->formula_;
        a2body_(i2body) = 0.;
        b2body_(i2body) = 0.;
        c2body_(i2body) = 0.;
        Tmin_2body_(i2body) = pr->Tmin_;
        Tmax_2body_(i2body) = pr->Tmax_;
        k2body_(i2body) = 0.;
        id7map_(pr->id_) = i2body;
        id7type_(pr->id_) = ReactionType::twobody;
        i2body++;
      } else {
        error = true;
      }
    //------- twobodytr - 2body reaction with temperature range ----------
    } else if (rtype == ReactionType::twobodytr) {
      const int nprev = std::count(idtr.begin(), idtr.end(), pr->id_);
      if (nprev == 0) {
        in2bodytr1_(i2bodytr) = ispec_map_[ pr->reactants_[0]];
        in2bodytr2_(i2bodytr) = ispec_map_[ pr->reactants_[1]];
        out2bodytr1_(i2bodytr) = ispec_map_[ pr->products_[0]];
        if (pr->products_.size() >= 2) {
          out2bodytr2_(i2bodytr) = ispec_map_[ pr->products_[1]];
        } else {
          out2bodytr2_(i2bodytr) = -1;
        }
        if (pr->products_.size() >= 3) {
          out2bodytr3_(i2bodytr) = ispec_map_[ pr->products_[2]];
        } else {
          out2bodytr3_(i2bodytr) = -1;
        }
        k2bodytr_(i2bodytr) = 0.;
        i2bodytr++;
      } 
      nr_2bodytr_(i2bodytr-1) = nprev + 1;
      frml_2bodytr_(i2bodytr-1, nprev) = pr->formula_;
      a2bodytr_(i2bodytr-1, nprev) = pr->alpha_;
      b2bodytr_(i2bodytr-1, nprev) = pr->beta_;
      c2bodytr_(i2bodytr-1, nprev) = pr->gamma_;
      Tmin_2bodytr_(i2bodytr-1, nprev) = pr->Tmin_;
      Tmax_2bodytr_(i2bodytr-1, nprev) = pr->Tmax_;
      idtr.push_back(pr->id_);

    //-------------------- grain - grain assisted reaction ----------------
    } else if (rtype == ReactionType::grain) {
      if (pr->reactants_[0] == "H" && pr->reactants_[1] == "H") {
        igr_H_ = igr;
      }
      if (pr->formula_ == 7) {
        ingr1_(igr) = ispec_map_[pr->reactants_[0]];
        ingr2_(igr) = ispec_map_[pr->reactants_[1]];
        outgr_(igr) = ispec_map_[pr->products_[0]];
        id7map_(pr->id_) = igr;
        id7type_(pr->id_) = ReactionType::grain;
        kgr_(igr) = 0.;
        igr++;
      } else{
        error = true;
      }

    //------------------ special - special reactions -----------------------
    } else if (rtype == ReactionType::special) {
      if (pr->formula_ == 7) {
        id7map_(pr->id_) = isr;
        id7type_(pr->id_) = ReactionType::special;
        ksr_(isr) = 0.;
        for (int jin=0; jin<n_insr_; jin++) {
          if (jin < pr->reactants_.size()) {
            insr_(isr, jin) = ispec_map_[pr->reactants_[jin]];
          } else {
            insr_(isr, jin) = -1;
          }
        }
        for (int jout=0; jout<n_outsr_; jout++) {
          if (jout < pr->products_.size()) {
            outsr_(isr, jout) = ispec_map_[pr->products_[jout]];
          } else {
            outsr_(isr, jout) = -1;
          }
        }
        isr++;
      } else {
        error = true;
      }

    //------------------ formula not recogonized -------------------------
    } else{
      error = true;
    }
    if (error) {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
          << std::endl
          << "reaction ID=" << pr->id_ << ", itype=" << pr->itype_ 
          << " and forumla=" << pr->formula_ << " undefined." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  //sanity check
  if (icr != n_cr_ || icrp != n_crp_ || iph != n_ph_ || i2body != n_2body_ 
      || i2bodytr != n_2bodytr_ || igr != n_gr_ || isr != n_sr_) {
    std::stringstream msg; 
    msg << "### FATAL ERROR in ChemNetwork InitializeReactions() [ChemNetwork]"
      << ": counts of reactions does not match." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

void ChemNetwork::UpdateRates(const Real y[NSCALARS], const Real E) {
  const Real y_H2 = y[ispec_map_["H2"]];
  const Real y_H = y[ispec_map_["H"]];
  Real y_e;
  if (ispec_map_.find("e-") != ispec_map_.end()) {
    y_e = y[ispec_map_["e-"]];
  } else {
    y_e = 0.;
  }

  Real T, Tcap;
  if (NON_BAROTROPIC_EOS) {
    T = E / Thermo::CvCold(y_H2, xHe_, y_e);
  } else {
    //isohermal EOS
    T = temperature_;
  }
	//cap T above some minimum temperature
	if (T < temp_min_rates_) {
		T = temp_min_rates_;
	} 
	//cosmic ray reactions
	for (int i=0; i<n_cr_; i++) {
		kcr_(i) = kcr_base_(i) * rad_(index_cr_);
	}

	//cosmic ray induced photo reactions
	for (int i=0; i<n_crp_; i++) {
		kcrp_(i) = kcrp_base_(i) * rad_(index_cr_) * 2*y_H2;
	}

	//FUV reactions
	for (int i=0; i<n_ph_; i++) {
    kph_(i) = kph_base_(i) * rad_(i);
	}

  //2body reactions
  if (is_Tcap_2body_) {
    for (int i=0; i<n_2body_; i++) {
      if (T < Tmin_2body_(i)) {
        Tcap = Tmin_2body_(i);
      } else if (T > Tmax_2body_(i)) {
        Tcap = Tmax_2body_(i);
      } else {
        Tcap = T;
      }
      if (frml_2body_(i) == 3) {
        k2body_(i) = a2body_(i)*pow(Tcap/300., b2body_(i))*exp(-c2body_(i)/Tcap) * nH_;
      } else if (frml_2body_(i) == 4) {
        k2body_(i) = a2body_(i)*b2body_(i)*( 0.62 
                                          + 0.4767*c2body_(i)*sqrt(300./Tcap) ) * nH_;
      } else if (frml_2body_(i) == 5) {
        k2body_(i) = a2body_(i)*b2body_(i)*( 1 + 0.0967*c2body_(i)*sqrt(300./Tcap) 
                                           + 28.501*c2body_(i)*c2body_(i)/Tcap ) * nH_;
      }
    }
  } else {
    for (int i=0; i<n_2body_; i++) {
      if (frml_2body_(i) == 3) {
        k2body_(i) = a2body_(i)*pow(T/300., b2body_(i))*exp(-c2body_(i)/T) * nH_;
      } else if (frml_2body_(i) == 4) {
        k2body_(i) = a2body_(i)*b2body_(i)*( 0.62 
                                            + 0.4767*c2body_(i)*sqrt(300./T) ) * nH_;
      } else if (frml_2body_(i) == 5) {
        k2body_(i) = a2body_(i)*b2body_(i)*( 1 + 0.0967*c2body_(i)*sqrt(300./T) 
                                             + 28.501*c2body_(i)*c2body_(i)/T ) * nH_;
      }
    }
  }

  //2bodytr reactions
  if (is_Tcap_2body_) {
    for (int i=0; i<n_2bodytr_; i++) {
      int nr = nr_2bodytr_(i);
      int irange1 = 0;
      int irange2 = 0;
      Real rate1 = 0.;
      Real rate2 = 0.;
      if ( T < Tmin_2bodytr_(i,0) ) {
        Tcap = Tmin_2bodytr_(i,0);
      } else if ( T > Tmax_2bodytr_(i,nr-1) ) {
        Tcap = Tmax_2bodytr_(i,nr-1);
      } else {
        Tcap = T;
      }
      //select which temperature range to use
      if ( Tcap <= Tmax_2bodytr_(i,0) ) {
        irange1 = 0;
        irange2 = 0;
      } else if ( Tcap <= Tmin_2bodytr_(i,1) ) {
        irange1 = 0;
        irange2 = 1;
      } else if ( Tcap <= Tmax_2bodytr_(i,1) ) {
        irange1 = 1;
        irange2 = 1;
      } else {
        if (nr == 2) {
          irange1 = 1;
          irange2 = 1;
        } else if (nr == 3) {
          if ( Tcap <= Tmin_2bodytr_(i,2) ) {
            irange1 = 1;
            irange2 = 2;
          } else {
            irange1 = 2;
            irange2 = 2;
          }
        } else {
          std::stringstream msg; 
          msg << "### fatal error in chemnetwork UpdateRates() [chemnetwork]: "
            << "2bodytr reaction with more than 3 temperature ranges not implemented."
            << std::endl; 
          ATHENA_ERROR(msg);
        }
      }
      //calculate rates
      if (frml_2bodytr_(i,irange1) == 3) {
        rate1 = a2bodytr_(i,irange1)*pow(Tcap/300., b2bodytr_(i,irange1))
                    *exp(-c2bodytr_(i,irange1)/Tcap) * nH_;
      } else if (frml_2bodytr_(i,irange1) == 4) {
        rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*( 0.62 
                               + 0.4767*c2bodytr_(i,irange1)*sqrt(300./Tcap) ) * nH_;
      } else if (frml_2bodytr_(i,irange1) == 5) {
        rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*(
            1 + 0.0967*c2bodytr_(i,irange1)*sqrt(300./Tcap) 
              + 28.501*c2bodytr_(i,irange1)*c2bodytr_(i,irange1)/Tcap ) * nH_;
      }
      if (irange1 == irange2) {
        rate2 = rate1;
      } else {
        if (frml_2bodytr_(i,irange2) == 3) {
          rate2 = a2bodytr_(i,irange2)*pow(Tcap/300., b2bodytr_(i,irange2))
                      *exp(-c2bodytr_(i,irange2)/Tcap) * nH_;
        } else if (frml_2bodytr_(i,irange2) == 4) {
          rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)*( 0.62 
                                 + 0.4767*c2bodytr_(i,irange2)*sqrt(300./Tcap) ) * nH_;
        } else if (frml_2bodytr_(i,irange2) == 5) {
          rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)*(
              1 + 0.0967*c2bodytr_(i,irange2)*sqrt(300./Tcap) 
                + 28.501*c2bodytr_(i,irange2)*c2bodytr_(i,irange2)/Tcap ) * nH_;
        }
      }
      //assign reaction rate
      k2bodytr_(i) = (rate1 + rate2) * 0.5;
    }
  } else {
    for (int i=0; i<n_2bodytr_; i++) {
      int nr = nr_2bodytr_(i);
      int irange1 = 0;
      int irange2 = 0;
      Real rate1 = 0.;
      Real rate2 = 0.;
      //select which temperature range to use
      if ( T <= Tmax_2bodytr_(i,0) ) {
        irange1 = 0;
        irange2 = 0;
      } else if ( T <= Tmin_2bodytr_(i,1) ) {
        irange1 = 0;
        irange2 = 1;
      } else if ( T <= Tmax_2bodytr_(i,1) ) {
        irange1 = 1;
        irange2 = 1;
      } else {
        if (nr == 2) {
          irange1 = 1;
          irange2 = 1;
        } else if (nr == 3) {
          if ( T <= Tmin_2bodytr_(i,2) ) {
            irange1 = 1;
            irange2 = 2;
          } else {
            irange1 = 2;
            irange2 = 2;
          }
        } else {
          std::stringstream msg; 
          msg << "### fatal error in chemnetwork UpdateRates() [chemnetwork]: "
            << "2bodytr reaction with more than 3 temperature ranges not implemented."
            << std::endl; 
          ATHENA_ERROR(msg);
        }
      }
      //calculate rates
      if (frml_2bodytr_(i,irange1) == 3) {
        rate1 = a2bodytr_(i,irange1)*pow(T/300., b2bodytr_(i,irange1))
                    *exp(-c2bodytr_(i,irange1)/T) * nH_;
      } else if (frml_2bodytr_(i,irange1) == 4) {
        rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*( 0.62 
                               + 0.4767*c2bodytr_(i,irange1)*sqrt(300./T) ) * nH_;
      } else if (frml_2bodytr_(i,irange1) == 5) {
        rate1 = a2bodytr_(i,irange1)*b2bodytr_(i,irange1)*(
            1 + 0.0967*c2bodytr_(i,irange1)*sqrt(300./T) 
              + 28.501*c2bodytr_(i,irange1)*c2bodytr_(i,irange1)/T ) * nH_;
      }
      if (irange1 == irange2) {
        rate2 = rate1;
      } else {
        if (frml_2bodytr_(i,irange2) == 3) {
          rate2 = a2bodytr_(i,irange2)*pow(T/300., b2bodytr_(i,irange2))
                      *exp(-c2bodytr_(i,irange2)/T) * nH_;
        } else if (frml_2bodytr_(i,irange2) == 4) {
          rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)*( 0.62 
                                 + 0.4767*c2bodytr_(i,irange2)*sqrt(300./T) ) * nH_;
        } else if (frml_2bodytr_(i,irange2) == 5) {
          rate2 = a2bodytr_(i,irange2)*b2bodytr_(i,irange2)*(
              1 + 0.0967*c2bodytr_(i,irange2)*sqrt(300./T) 
                + 28.501*c2bodytr_(i,irange2)*c2bodytr_(i,irange2)/T ) * nH_;
        }
      }
      //assign reaction rate
      k2bodytr_(i) = (rate1 + rate2) * 0.5;
    }
  }

  //special rates and grain assisted reactions
  UpdateRatesSpecial(y, E);
  return;
}

//sort the type of the reaction, check format
ReactionType ChemNetwork::SortReaction(KidaReaction* pr) const {
  //---------------- 1 - direct cosmic-ray ionization --------------
  if (pr->itype_ == 1) {
    //check format of reaction 
    if (pr->reactants_.size() == 2 && pr->products_.size() == 2
        && (pr->reactants_[0] == "CR" || pr->reactants_[1] == "CR") ) {
    } else {
      std::stringstream msg; 
      msg << "### fatal error in chemnetwork sortreaction() [chemnetwork]"
          << std::endl << "wrong format in cr reaction id=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::cr;

  //---------------- 2 - cosmic-ray induced photo ionization --------
  } else if (pr->itype_ == 2) {
    //check format of reaction 
    if (pr->reactants_.size() == 2 && pr->products_.size() == 2
        && (pr->reactants_[0] == "CRP" || pr->reactants_[1] == "CRP") ) {
    } else {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
         << std::endl << "Wrong format in CRP reaction ID=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::crp;

  //---------------- 3 - FUV ionization/dissociation ----------------
  } else if (pr->itype_ == 3) {
    //check format of reaction 
    if (pr->reactants_.size() == 2 && pr->products_.size() == 2
        && (pr->reactants_[0] == "Photon" || pr->reactants_[1] == "Photon") ) {
    } else {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
         << std::endl << "Wrong format in FUV reaction ID=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::photo;

  //---------------- 4-8 - 2body reaction ---------------------------
  } else if (pr->itype_ >= 4 && pr->itype_ <= 8) {
    //check format
    if (std::find(id_2bodytr_.begin(), id_2bodytr_.end(), pr->id_)
                  == id_2bodytr_.end()) {
      if (pr->reactants_.size() != 2 || 
          (pr->products_.size() != 1 && pr->products_.size() != 2
           && pr->products_.size() != 3 && pr->products_.size() != 4)) {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
            << std::endl << "Wrong format in 2body reaction ID=" << pr->id_
            << std::endl;
        ATHENA_ERROR(msg);
      }
      return ReactionType::twobody;
    } else { //2 body reaction with temperature range
      if (pr->reactants_.size() != 2 || 
          (pr->products_.size() != 1 && pr->products_.size() != 2
           && pr->products_.size() != 3)) {
        std::stringstream msg; 
        msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
            << std::endl << "Wrong format in 2bodytr reaction ID=" << pr->id_
            << std::endl;
        ATHENA_ERROR(msg);
      }
      return ReactionType::twobodytr;
    }

  //-------------------- 9 - grain assisted reaction ----------------
  } else if (pr->itype_ == 9) {
    //check format
    if (pr->reactants_.size() != 2 || pr->products_.size() != 1) {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in gr reaction ID=" << pr->id_ << std::endl;
      ATHENA_ERROR(msg);
    }
    if (pr->reactants_[1] != "e-" && pr->reactants_[1] != "H") {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "second reactant must be H or e-." << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::grain;

  //------------------ 10 - special reactions -----------------------
  } else if (pr->itype_ == 10) {
    if (pr->reactants_.size() > n_insr_ || pr->products_.size() > n_outsr_) {
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
          << std::endl << "Wrong format in special reaction ID=" << pr->id_ 
          << std::endl;
      ATHENA_ERROR(msg);
    }
    return ReactionType::special;

  //------------------ type not recogonized -------------------------
  } else{
    std::stringstream msg; 
    msg << "### FATAL ERROR in ChemNetwork SortReaction() [ChemNetwork]"
      << std::endl << "Wrong format in special reaction ID=" << pr->id_ 
      << std::endl;
    ATHENA_ERROR(msg);
  }
}

void ChemNetwork::CheckReaction(KidaReaction reaction) {
  int atom_count_in[KidaSpecies::natom_];
  int atom_count_out[KidaSpecies::natom_];
  int charge_in = 0;
  int charge_out = 0;
  for (int ia=0; ia<KidaSpecies::natom_; ia++) {
    atom_count_in[ia] = 0;
    atom_count_out[ia] = 0;
  }
  for (int i=0; i<reaction.reactants_.size(); i++) {
    if (reaction.reactants_[i] == "CR" || reaction.reactants_[i] == "CRP" 
        || reaction.reactants_[i] == "Photon") {
      continue;
    }
    for (int ia=0; ia<KidaSpecies::natom_; ia++) {
      atom_count_in[ia] +=
        species_[ispec_map_[reaction.reactants_[i]]].atom_count_[ia];
    }
    charge_in += species_[ispec_map_[reaction.reactants_[i]]].charge_;
  }

  for (int i=0; i<reaction.products_.size(); i++) {
    for (int ia=0; ia<KidaSpecies::natom_; ia++) {
      atom_count_out[ia] +=
        species_[ispec_map_[reaction.products_[i]]].atom_count_[ia];
    }
    charge_out += species_[ispec_map_[reaction.products_[i]]].charge_;
  }

  if (charge_in != charge_out) {
    reaction.Print();
    std::stringstream msg; 
    msg << "### FATAL ERROR in ChemNetwork CheckReaction() [ChemNetwork] :"
        << "charge not conserved." << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int ia=0; ia<KidaSpecies::natom_; ia++) {
    if (atom_count_in[ia] != atom_count_out[ia]) {
      reaction.Print();
      std::stringstream msg; 
      msg << "### FATAL ERROR in ChemNetwork CheckReaction() [ChemNetwork] :"
          << "atoms not conserved." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  return;
}


void ChemNetwork::PrintProperties() const {
  //print each species.
  for (int i=0; i<NSCALARS; i++) {
    std::cout << "species: " << i << std::endl;
    std::cout << "name=" << species_[i].name << ", index=" << species_[i].index 
      << ", charge=" << species_[i].charge_ << std::endl;
    std::cout << "atom_count_ = ";
    for (int j=0; j < species_[i].natom_; j++) {
      std::cout << species_[i].atom_count_[j] << " ";
    }
    std::cout << std::endl;
  }

  //print each reactions.
  std::cout << "number of reactions: " << nr_ << std::endl;
  for (int i=0; i<reactions_.size(); i++) {
    reactions_[i].Print();
    std::cout << "alpha=" << reactions_[i].alpha_ << "," 
              << "beta=" << reactions_[i]. beta_ << "," 
              << "gamma=" << reactions_[i].gamma_ << "," 
              << "Tmin=" << reactions_[i].Tmin_ << "," 
              << "Tmax=" << reactions_[i].Tmax_ << "," 
              << "itype=" << reactions_[i].itype_ << ","
              << "forumla=" << reactions_[i].formula_ << std::endl;
  }

  //print reaction coefficients
  //cosmic-ray reactions
  std::cout << "CR reations:" << std::endl;
  for (int i=0; i<n_cr_; i++) {
    std::cout<< species_names[incr_(i)] << " + CR -> "
      << species_names[outcr1_(i)] << " + " << species_names[outcr2_(i)] << ", "
      << "kcr_base_=" << kcr_base_(i) << std::endl;
  }

  //cosmic-ray induced photo reactions
  std::cout << "CRP reations:" << std::endl;
  for (int i=0; i<n_crp_; i++) {
    std::cout<< species_names[incrp_(i)] << " + CRP -> "
      << species_names[outcrp1_(i)] << " + " << species_names[outcrp2_(i)] << ", "
      << "kcrp_base_=" << kcrp_base_(i) << std::endl;
  }

  //FUV reactions
  std::cout << "FUV photo- ionization/dissociation:" << std::endl;
  for (int i=0; i<n_ph_; i++) {
    std::cout<< species_names[inph_(i)] << " + Photon -> "
      << species_names[outph1_(i)] << " + " << species_names[outph2_(i)] << ", "
      << "kph_base_=" << kph_base_(i) << ", kph_avfac_=" << kph_avfac_(i)
      << std::endl;
  }
  std::cout << "smap_ph_: " << std::endl;
  for (std::map<std::string,int>::const_iterator it=smap_ph_.begin();
       it!=smap_ph_.end(); it++) {
    std::cout << it->first << " => " << it->second << std::endl;
  }

  //2body reactions
  std::cout << "2body reactions:" << std::endl;
  for (int i=0; i<n_2body_; i++) {
    std::cout<< species_names[in2body1_(i)] << " + "
      << species_names[in2body2_(i)]<< " -> " << species_names[out2body1_(i)];
    if (out2body2_(i) >= 0) {
      std::cout<< " + " << species_names[out2body2_(i)];
    }
    if (out2body3_(i) >= 0) {
      std::cout<< " + " << species_names[out2body3_(i)];
    }
    if (out2body4_(i) >= 0) {
      std::cout<< " + " << species_names[out2body4_(i)];
    }
    std::cout<< ", " << "alpha=" << a2body_(i) << ", beta=" << b2body_(i)
      << ", gamma=" << c2body_(i) << ", Trange=[" << Tmin_2body_(i) << "," 
      << Tmax_2body_(i) << "]" << std::endl;
  }

  //2body reactions with temperature ranges
  std::cout << "2bodytr reactions:" << std::endl;
  for (int i=0; i<n_2bodytr_; i++) {
    std::cout<< species_names[in2bodytr1_(i)] << " + "
      << species_names[in2bodytr2_(i)]<< " -> " << species_names[out2bodytr1_(i)];
    if (out2bodytr2_(i) >= 0) {
      std::cout<< " + " << species_names[out2bodytr2_(i)];
    }
    if (out2bodytr3_(i) >= 0) {
      std::cout<< " + " << species_names[out2bodytr3_(i)];
    }
    std::cout << "    ,nr_2bodytr_=" << nr_2bodytr_(i) << std::endl;
    for (int j=0; j<nr_2bodytr_(i); j++) {
      std::cout<< "alpha=" << a2bodytr_(i, j) << ", beta=" 
        << b2bodytr_(i, j) << ", gamma=" << c2bodytr_(i, j) 
        << ", Trange=[" << Tmin_2bodytr_(i, j) << "," << Tmax_2bodytr_(i, j) 
        << "], formula=" << frml_2bodytr_(i, j) << std::endl;
    }
  }

  //grain assisted reactions
  std::cout << "gr reations:" << std::endl;
  for (int i=0; i<n_gr_; i++) {
    std::cout<< species_names[ingr1_(i)] << " + " << species_names[ingr2_(i)] 
      <<" -> " << species_names[outgr_(i)] << std::endl;
  }

  //special reactions
  std::cout << "special reations:" << std::endl;
  for (int i=0; i<n_sr_; i++) {
    for (int jin=0; jin<n_insr_; jin++) {
      if (insr_(i, jin) >= 0) {
        std::cout << species_names[insr_(i, jin)];
        if (jin < n_insr_-1 && insr_(i, jin+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << " -> ";
    for (int jout=0; jout<n_outsr_; jout++) {
      if (outsr_(i, jout) >= 0) {
        std::cout << species_names[outsr_(i, jout)];
        if (jout < n_outsr_-1 && outsr_(i, jout+1) >= 0) {
          std::cout << " + ";
        }
      }
    }
    std::cout << std::endl;
  }
  
  for (int i=0; i<id7max_+1; i++) {
    if (id7map_(i) >= 0) {
      std::cout << i << " => " << id7map_(i) << ", ";
      switch (id7type_(i)){
        case ReactionType::cr: std::cout << "cr"; break;
        case ReactionType::crp: std::cout << "crp"; break;
        case ReactionType::twobody: std::cout << "2body"; break;
        case ReactionType::grain: std::cout << "gr"; break;
        case ReactionType::special: std::cout << "sr"; break;
        default: std::stringstream msg; 
                 msg << "### FATAL ERROR in ChemNetwork PrintPropeties() "
                   << "[ChemNetwork]: reaction type not recognized for special rate."
                   << std::endl;
                 ATHENA_ERROR(msg);
                 break;
      }
      std::cout << std::endl;
    }
  }

  return;
}

void ChemNetwork::OutputRates(FILE *pf) const {
  //output the reactions and rates
	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "%4s + CR -> %4s + %4s,     kcr = %.2e\n", 
		 species_names[incr_(i)].c_str(), species_names[outcr1_(i)].c_str(),
     species_names[outcr2_(i)].c_str(), kcr_(i));
	}
	for (int i=0; i<n_crp_; i++) {
		fprintf(pf, "%4s + CRP -> %4s + %4s,     kcrp = %.2e\n", 
		 species_names[incrp_(i)].c_str(), species_names[outcrp1_(i)].c_str(),
     species_names[outcrp2_(i)].c_str(), kcrp_(i));
	}
	for (int i=0; i<n_ph_; i++) {
		fprintf(pf, "%4s + Photon -> %4s + %4s,     kph = %.2e\n", 
		 species_names[inph_(i)].c_str(), species_names[outph1_(i)].c_str(),
     species_names[outph2_(i)].c_str(), kph_(i));
	}
	for (int i=0; i<n_2body_; i++) {
    fprintf(pf, "%4s + %4s -> %4s",
        species_names[in2body1_(i)].c_str(),
        species_names[in2body2_(i)].c_str(),
        species_names[out2body1_(i)].c_str());
    if (out2body2_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2body2_(i)].c_str());
    }
    if (out2body3_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2body3_(i)].c_str());
    }
    if (out2body4_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2body4_(i)].c_str());
    }
    fprintf(pf,   ",     k2body = %.2e\n", k2body_(i));
	}
	for (int i=0; i<n_2bodytr_; i++) {
    fprintf(pf, "%4s + %4s -> %4s",
        species_names[in2bodytr1_(i)].c_str(),
        species_names[in2bodytr2_(i)].c_str(),
        species_names[out2bodytr1_(i)].c_str());
    if (out2bodytr2_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2bodytr2_(i)].c_str());
    }
    if (out2bodytr3_(i) >= 0) {
      fprintf(pf, " + %4s", species_names[out2bodytr3_(i)].c_str());
    }
    fprintf(pf,   ",     k2bodytr = %.2e\n", k2bodytr_(i));
	}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "%4s + %4s (+ gr) -> %4s (+ gr),       kgr = %.2e\n", 
		 species_names[ingr1_(i)].c_str(), species_names[ingr2_(i)].c_str(),
     species_names[outgr_(i)].c_str(), kgr_(i));
	}
  for (int i=0; i<n_sr_; i++) {
    for (int jin=0; jin<n_insr_; jin++) {
      if (insr_(i, jin) >= 0) {
        fprintf(pf, "%4s", species_names[insr_(i, jin)].c_str());
        if (jin < n_insr_-1 && insr_(i, jin+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
    fprintf(pf, " -> ");
    for (int jout=0; jout<n_outsr_; jout++) {
      if (outsr_(i, jout) >= 0) {
        fprintf(pf, "%4s", species_names[outsr_(i, jout)].c_str());
        if (jout < n_outsr_-1 && outsr_(i, jout+1) >= 0) {
          fprintf(pf, " + ");
        }
      }
    }
		fprintf(pf, ",       ksr = %.2e\n", ksr_(i));
  }
  return;
}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  Real rho, rho_floor, rad_sum;
  const int nang = pmy_mb_->prad->nang;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  //hydrogen atom number density
  nH_ =  rho * unit_density_in_nH_;
  //average radiation field of all angles
  for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
    rad_sum = 0;
    //radiation
    for (int iang=0; iang < nang; ++iang) {
      rad_sum += pmy_mb_->prad->ir(k, j, i, ifreq * nang + iang);
    }
    if (ifreq == index_cr_) {
      rad_(index_cr_) = rad_sum / float(nang);
    } else {
      rad_(ifreq) = rad_sum * unit_radiation_in_draine1987_ / float(nang) ;
    }
#ifdef DEBUG
    if (isnan(rad_(ifreq))) {
      printf("InitializeNextStep: ");
      printf("ifreq=%d, nang=%d, rad_sum=%.2e\n", ifreq, nang, rad_sum);
      OutputRates(stdout);
    }
#endif
  }
  //CO cooling paramters
  SetGrad_v(k, j, i);
  return;
}

void ChemNetwork::RHS(const Real t, const Real y[NSCALARS], const Real ED,
                      Real ydot[NSCALARS]){
	Real rate = 0;
  Real E_ergs = ED * unit_E_in_cgs_ / nH_; //ernergy per hydrogen atom
	//store previous y includeing negative abundance correction
	Real y0[NSCALARS];//correct negative abundance, only for UpdateRates()
	Real ydotg[NSCALARS];

	for(int i=0; i<NSCALARS; i++) {
		ydotg[i] = 0.0;
  }

  //correct negative abundance to zero, used in rate update
  for (int i=0; i<NSCALARS; i++) {
    if (y[i] < 0) {
      y0[i] = 0;
    } else {
      y0[i] = y[i];
    }
    //throw error if nan, or inf, or large negative value occurs
    if ( isnan(y[i]) || isinf(y[i]) ) {
      printf("RHS: ");
      for (int j=0; j<NSCALARS; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), y[j]);
      }
      printf("\n");
      OutputRates(stdout);
      printf("rad_ = ");
      for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
        printf("%.2e  ", rad_(ifreq));
      }
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      std::stringstream msg;
      msg << "ChemNetwork (kida): RHS(y): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  UpdateRates(y0, E_ergs);

#ifdef DEBUG
  if (output_rates) {
    FILE *pf = fopen("chem_network.dat", "w");
    OutputRates(pf);
    fclose(pf);
    output_rates = false;
  }
#endif

  //cosmic ray reactions
  for (int i=0; i<n_cr_; i++) {
    rate = kcr_(i) * y[incr_(i)];
    ydotg[incr_(i)] -= rate;
    ydotg[outcr1_(i)] += rate;
    ydotg[outcr2_(i)] += rate;
  }

  //cosmic ray induced photo reactions
  for (int i=0; i<n_crp_; i++) {
    rate = kcrp_(i) * y[incrp_(i)];
    ydotg[incrp_(i)] -= rate;
    ydotg[outcrp1_(i)] += rate;
    ydotg[outcrp2_(i)] += rate;
  }

  //cosmic ray induced photo reactions
  for (int i=0; i<n_ph_; i++) {
    rate = kph_(i) * y[inph_(i)];
    ydotg[inph_(i)] -= rate;
    ydotg[outph1_(i)] += rate;
    ydotg[outph2_(i)] += rate;
  }

  //2body reactions
  for (int i=0; i<n_2body_; i++) {
    rate =  k2body_(i) * y[in2body1_(i)] * y[in2body2_(i)];
    if (y[in2body1_(i)] < 0 && y[in2body2_(i)] < 0) {
      rate *= -1.;
    }
    ydotg[in2body1_(i)] -= rate;
    ydotg[in2body2_(i)] -= rate;
    ydotg[out2body1_(i)] += rate;
    if (out2body2_(i) >= 0) {
      ydotg[out2body2_(i)] += rate;
    }
    if (out2body3_(i) >= 0) {
      ydotg[out2body3_(i)] += rate;
    }
    if (out2body4_(i) >= 0) {
      ydotg[out2body4_(i)] += rate;
    }
  }

  //2bodytr reactions
  for (int i=0; i<n_2bodytr_; i++) {
    rate =  k2bodytr_(i) * y[in2bodytr1_(i)] * y[in2bodytr2_(i)];
    if (y[in2bodytr1_(i)] < 0 && y[in2bodytr2_(i)] < 0) {
      rate *= -1.;
    }
    ydotg[in2bodytr1_(i)] -= rate;
    ydotg[in2bodytr2_(i)] -= rate;
    ydotg[out2bodytr1_(i)] += rate;
    if (out2bodytr2_(i) >= 0) {
      ydotg[out2bodytr2_(i)] += rate;
    }
    if (out2bodytr3_(i) >= 0) {
      ydotg[out2bodytr3_(i)] += rate;
    }
  }

  //grain assisted reactions
  for (int i=0; i<n_gr_; i++) {
    rate = kgr_(i) * y[ingr1_(i)];
    ydotg[ingr1_(i)] -= rate;
    ydotg[ingr2_(i)] -= rate;
    ydotg[outgr_(i)] += rate;
  }

  //special reactions
  for (int i=0; i<n_sr_; i++) {
    rate = ksr_(i);
    for (int jin=0; jin<n_insr_; jin++) {
      if (insr_(i, jin) >= 0) {
        ydotg[insr_(i, jin)] -= rate;
      }
    }
    for (int jout=0; jout<n_outsr_; jout++) {
      if (outsr_(i, jout) >= 0) {
        ydotg[outsr_(i, jout)] += rate;
      }
    }
  }

	//set ydot to return
	for (int i=0; i<NSCALARS; i++) {
    //return in code units
		ydot[i] = ydotg[i] * unit_time_in_s_;
	}

  //throw error if nan, or inf, or large value occurs
  for (int i=0; i<NSCALARS; i++) {
    if ( isnan(ydot[i]) || isinf(ydot[i]) ) {
      printf("ydot: ");
      for (int j=0; j<NSCALARS; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), ydot[j]);
      }
      printf("abundances: ");
      for (int j=0; j<NSCALARS; j++) {
        printf("%s: %.2e  ", species_names[j].c_str(), y[j]);
      }
      printf("\n");
      OutputRates(stdout);
      printf("rad_ = ");
      for (int ifreq=0; ifreq < n_freq_; ++ifreq) {
        printf("%.2e  ", rad_(ifreq));
      }
      printf("\n");
      printf("nH_ = %.2e\n", nH_);
      printf("ED = %.2e\n", ED);
      printf("E_ergs = %.2e\n", E_ergs);
      printf("unit_E_in_cgs_ = %.2e\n", unit_E_in_cgs_);
      Real y_e;
      if (ispec_map_.find("e-") != ispec_map_.end()) {
        y_e = y0[ispec_map_["e-"]];
      } else {
        y_e = 0.;
      }
      printf("T = %.2e\n", E_ergs/Thermo::CvCold(y0[ispec_map_["H2"]], xHe_, y_e));
      std::stringstream msg;
      msg << "ChemNetwork (kida): RHS(ydot): nan or inf" << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED){
  Real E_ergs = ED * unit_E_in_cgs_ / nH_; //ernergy per hydrogen atom
  //isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }
  Real T = 0.;
  Real dEdt = 0.;
	Real y0[NSCALARS];
  //correct negative abundance to zero
  for (int i=0; i<NSCALARS; i++) {
    if (y[i] < 0) {
      y0[i] = 0;
    } else {
      y0[i] = y[i];
    }
  }
  //abundances
  const Real y_H2 = y0[ispec_map_["H2"]];
  const Real y_H = y0[ispec_map_["H"]];
  Real y_e, y_He, y_Hplus, kcr_H, kcr_H2, kcr_He, kgr_H, kph_H2;
  if (ispec_map_.find("e-") != ispec_map_.end()) {
    y_e = y0[ispec_map_["e-"]];
  } else {
    y_e = 0.;
  }
  if (ispec_map_.find("He") != ispec_map_.end()) {
    y_He = y0[ispec_map_["He"]];
  } else {
    y_He = 0.;
  }
  if (ispec_map_.find("H+") != ispec_map_.end()) {
    y_Hplus = y0[ispec_map_["H+"]];
  } else {
    y_Hplus = 0.;
  }
  if (icr_H_ >= 0) {
    kcr_H = kcr_(icr_H_);
  } else {
    kcr_H = 0.;
  }
  if (icr_H2_ >= 0) {
    kcr_H2 = kcr_(icr_H2_);
  } else {
    kcr_H2 = 0.;
  }
  if (icr_He_ >= 0) {
    kcr_He = kcr_(icr_He_);
  } else {
    kcr_He = 0.;
  }
  if (igr_H_ >= 0) {
    kgr_H = kgr_(igr_H_);
  } else {
    kgr_H = 0.;
  }
  if (iph_H2_ >= 0) {
    kph_H2 = kph_(iph_H2_);
  } else {
    kph_H2 = 0.;
  }
  //temperature
  T = E_ergs / Thermo::CvCold(y_H2, xHe_, y_e);
  //apply temperature floor, incase of very small or negative energy
	if (T < temp_min_rates_) {
		T = temp_min_rates_;
  }


  //--------------------------heating-----------------------------
	Real GCR, GPE, GH2gr, dot_xH2_photo, GH2pump, GH2diss;
  //CR heating
  GCR = Thermo::HeatingCr(y_e,  nH_, y_H,  y_He,  y_H2, kcr_H, kcr_He, kcr_H2);
  //photo electric effect on dust
  GPE = Thermo::HeatingPE(rad_(index_gpe_), zdg_, T, nH_*y_e);
  //H2 formation on dust grains
  GH2gr = Thermo::HeatingH2gr(y_H,  y_H2, nH_, T, kgr_H);
  //H2 UV pumping
  dot_xH2_photo = kph_H2 * y_H2;
  GH2pump = Thermo::HeatingH2pump(y_H,  y_H2, nH_, T, dot_xH2_photo);
  //H2 photo dissiociation.
  GH2diss = Thermo::HeatingH2diss(dot_xH2_photo);

  //--------------------------cooling-----------------------------
	Real LCII, LCI, LOI, LHotGas, LCOR, LH2, LDust, LRec, LH2diss, LHIion;
	Real vth, nCO, grad_small;
  Real NCOeff, gradeff;
  Real k2body_H2_H, k2body_H2_H2, k2body_H_e;
  if (i2body_H2_H_ >= 0) {
    k2body_H2_H = k2body_(i2body_H2_H_);
  } else {
    k2body_H2_H = 0.;
  }
  if (i2body_H2_H2_ >= 0) {
    k2body_H2_H2 = k2body_(i2body_H2_H2_);
  } else {
    k2body_H2_H2 = 0.;
  }
  if (i2body_H_e_ >= 0) {
    k2body_H_e = k2body_(i2body_H_e_);
  } else {
    k2body_H_e = 0.;
  }
	if (T < temp_min_cool_) {
		LCII = 0.;
		LCI = 0;
		LOI = 0.;
		LHotGas = 0;
		LCOR = 0;
		LH2 = 0;
		LDust = 0;
		LRec = 0;
		LH2diss = 0;
		LHIion = 0;
	} else {
		// C+ fine structure line 
    if (ispec_map_.find("C+") != ispec_map_.end()) {
      LCII = Thermo::CoolingCII(y0[ispec_map_["C+"]],
                                nH_*y_H,  nH_*y_H2, nH_*y_e, T);
    } else {
      LCII = 0.;
    }
		// CI fine structure line 
    if (ispec_map_.find("C") != ispec_map_.end()) {
      LCI = Thermo::CoolingCI(y0[ispec_map_["C"]], nH_*y_H, nH_*y_H2, nH_*y_e, T);
    } else {
      LCI = 0.;
    }
		// OI fine structure line 
    if (ispec_map_.find("O") != ispec_map_.end()) {
      LOI = Thermo::CoolingOI(y0[ispec_map_["O"]], nH_*y_H, nH_*y_H2, nH_*y_e, T);
    } else {
      LOI = 0.;
    }
		// cooling of hot gas: radiative cooling, free-free.
		LHotGas = Thermo::CoolingHotGas(nH_, T, zdg_);
		// CO rotational lines 
    if (ispec_map_.find("CO") != ispec_map_.end()) {
      // Calculate effective CO column density
      Real y_CO = y0[ispec_map_["CO"]];
      vth = sqrt(2. * Thermo::kb_ * T / ChemistryUtility::mCO);
      nCO = nH_ * y_CO;
      grad_small = vth/Leff_CO_max_;
      gradeff = std::max(gradv_, grad_small);
      NCOeff = nCO / gradeff;
      LCOR = Thermo::CoolingCOR(y_CO, nH_*y_H,  nH_*y_H2, nH_*y_e, T, NCOeff);
    } else {
      LCOR = 0.;
    }
		// H2 vibration and rotation lines 
    LH2 = Thermo::CoolingH2(y_H2, nH_*y_H, nH_*y_H2, nH_*y_He,
                            nH_*y_Hplus, nH_*y_e, T);
		// dust thermo emission 
		LDust = Thermo::CoolingDustTd(zdg_, nH_, T, temp_dust_thermo_);
		// reconbination of e on PAHs 
		LRec = Thermo::CoolingRec(zdg_, T,  nH_*y_e, rad_(index_gpe_));
		// collisional dissociation of H2 
		LH2diss = Thermo::CoolingH2diss(y_H, y_H2, k2body_H2_H, k2body_H2_H2);
		// collisional ionization of HI 
		LHIion = Thermo::CoolingHIion(y_H,  y_e, k2body_H_e);
  }

  dEdt = (GCR + GPE + GH2gr + GH2pump + GH2diss)
          - (LCII + LCI + LOI + LHotGas + LCOR 
              + LH2 + LDust + LRec + LH2diss + LHIion);
  //return in code units
  Real dEDdt = dEdt * nH_ / unit_E_in_cgs_ * unit_time_in_s_;
	if ( isnan(dEdt) || isinf(dEdt) ) {
    if ( isnan(LCOR) || isinf(LCOR) ) {
      printf("NCOeff=%.2e, gradeff=%.2e, gradv_=%.2e, vth=%.2e, nH_=%.2e, nCO=%.2e\n",
          NCOeff, gradeff, gradv_, vth, nH_, nCO);
    }
		printf("GCR=%.2e, GPE=%.2e, GH2gr=%.2e, GH2pump=%.2e GH2diss=%.2e\n",
				GCR , GPE , GH2gr , GH2pump , GH2diss);
		printf("LCII=%.2e, LCI=%.2e, LOI=%.2e, LHotGas=%.2e, LCOR=%.2e\n",
				LCII , LCI , LOI , LHotGas , LCOR);
		printf("LH2=%.2e, LDust=%.2e, LRec=%.2e, LH2diss=%.2e, LHIion=%.2e\n",
				LH2 , LDust , LRec , LH2diss , LHIion);
		printf("T=%.2e, dEdt=%.2e, E=%.2e, dEergsdt=%.2e, E_ergs=%.2e, Cv=%.2e, nH=%.2e\n",
        T, dEDdt, ED, dEdt, E_ergs, Thermo::CvCold(y_H2, xHe_, y_e), nH_);
		for (int i=0; i<NSCALARS; i++) {
			printf("%s: %.2e  ", species_names[i].c_str(), y0[i]);
		}
		printf("\n");
    std::stringstream msg;
    msg << "ChemNetwork (kida): dEdt: nan or inf number" << std::endl;
    ATHENA_ERROR(msg);
	}
#ifdef DEBUG
  if (output_thermo) {
    printf("NCOeff=%.2e, gradeff=%.2e, gradv_=%.2e, vth=%.2e, nH_=%.2e, nCO=%.2e\n",
        NCOeff, gradeff, gradv_, vth, nH_, nCO);
		printf("GCR=%.2e, GPE=%.2e, GH2gr=%.2e, GH2pump=%.2e GH2diss=%.2e\n",
				GCR , GPE , GH2gr , GH2pump , GH2diss);
		printf("LCII=%.2e, LCI=%.2e, LOI=%.2e, LHotGas=%.2e, LCOR=%.2e\n",
				LCII , LCI , LOI , LHotGas , LCOR);
		printf("LH2=%.2e, LDust=%.2e, LRec=%.2e, LH2diss=%.2e, LHIion=%.2e\n",
				LH2 , LDust , LRec , LH2diss , LHIion);
		printf("T=%.2e, dEdt=%.2e, E=%.2e, dEergsdt=%.2e, E_ergs=%.2e, Cv=%.2e, nH=%.2e\n",
        T, dEDdt, ED, dEdt, E_ergs, Thermo::CvCold(y_H2, xHe_, y_e), nH_);
		for (int i=0; i<NSCALARS; i++) {
			printf("%s: %.2e  ", species_names[i].c_str(), y0[i]);
		}
		printf("\n");
    output_thermo = false;
  }
#endif
  return dEDdt;
}

void ChemNetwork::SetGrad_v(const int k, const int j, const int i) {
  AthenaArray<Real> &w = pmy_mb_->phydro->w;
  Real dvdx, dvdy, dvdz, dvdr_avg, di1, di2;
  Real dx1, dx2, dy1, dy2, dz1, dz2;
  Real dndx, dndy, dndz, gradn;
  //velocity gradient, same as LVG approximation in RADMC-3D when calculating
  //CO line emission.
  //vx
  di1 = w(IVX, k, j, i+1) - w(IVX, k, j, i);
  dx1 = ( pmy_mb_->pcoord->dx1f(i+1)+pmy_mb_->pcoord->dx1f(i) )/2.;
  di2 = w(IVX, k, j, i) - w(IVX, k, j, i-1);
  dx2 = ( pmy_mb_->pcoord->dx1f(i)+pmy_mb_->pcoord->dx1f(i-1) )/2.;
  dvdx = (di1/dx1 + di2/dx2)/2.;
  //vy
  di1 = w(IVY, k, j+1, i) - w(IVY, k, j, i);
  dy1 = ( pmy_mb_->pcoord->dx2f(j+1)+pmy_mb_->pcoord->dx2f(j) )/2.;
  di2 = w(IVY, k, j, i) - w(IVY, k, j-1, i);
  dy2 = ( pmy_mb_->pcoord->dx2f(j)+pmy_mb_->pcoord->dx2f(j-1) )/2.;
  dvdy = (di1/dy1 + di2/dy2)/2.;
  //vz
  di1 = w(IVZ, k+1, j, i) - w(IVZ, k, j, i);
  dz1 = ( pmy_mb_->pcoord->dx3f(k+1)+pmy_mb_->pcoord->dx3f(k) )/2.;
  di2 = w(IVZ, k, j, i) - w(IVZ, k-1, j, i);
  dz2 = ( pmy_mb_->pcoord->dx3f(k)+pmy_mb_->pcoord->dx3f(k-1) )/2.;
  dvdz = (di1/dz1 + di2/dz2)/2.;
  dvdr_avg = ( fabs(dvdx) + fabs(dvdy) + fabs(dvdz) ) / 3.;
  //asign gradv_, in cgs.
  gradv_ = dvdr_avg * unit_vel_in_cms_ / unit_length_in_cm_;
  return;
}


//default: no special rates
void __attribute__((weak)) ChemNetwork::UpdateRatesSpecial(const Real y[NSCALARS],
                                                           const Real E) {
  // do nothing
  return;
}