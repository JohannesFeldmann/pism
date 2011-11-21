// Copyright (C) 2007-2011 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __icePSTexModel_hh
#define __icePSTexModel_hh

#include <petsc.h>

#include "basal_resistance.hh"
#include "iceEISModel.hh"
#include "PISMMohrCoulombYieldStress.hh"

class DiagnosticTimeseries;

//! Derived class for Plastic till ice Stream with Thermocoupling (PST) experiments.
/*!
Runs numerical experiments reported in \ref BBssasliding .

There are four experiments, P1,P2,P3,P4.  Each has four (three for P2) ice
streams generated by lowered till friction angle.

(This derived class supercedes an older class which produced results presented
by Bueler at AGU 2007 and at NYU in Feb 2008.  IcePSTexModel has no "lake" or
"fjord", which were present in the old version.)
 */
class IcePSTexModel : public IceEISModel {

public:
  IcePSTexModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides);
  virtual ~IcePSTexModel();
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode allocate_basal_yield_stress();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode initFromFile(const char *fname);
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode additionalAtEndTimestep();

protected:
  char exper_chosen_name[10];
  int exper_chosen;
  
  PetscErrorCode setBedElev();

  // PST scalar diagnostic series
  PetscErrorCode prepare_series();
  char seriesname[PETSC_MAX_PATH_LEN];   // file name: "ser_pst_foo.nc" if "-o foo.nc"
  DiagnosticTimeseries
     *ivol, *iarea,                      // ice vol and area; mildly redundant 
     *dt_ser,                            // time step; mildly redundant 
     *maxcbar,                           // max speed anywhere
     *avup0, *avup1, *avup2, *avup3,     // upstream speeds; avup3 not written for P3
     *avdwn0, *avdwn1, *avdwn2, *avdwn3; // downstream speeds; avdwn3 DITTO
};

class PSTYieldStress : public PISMMohrCoulombYieldStress
{
public:
  PSTYieldStress(IceGrid &g, const NCConfigVariable &conf, int e, string name)
    : PISMMohrCoulombYieldStress(g, conf),
      experiment(e), experiment_name(name) {}
  virtual ~PSTYieldStress() {}

  virtual PetscErrorCode init(PISMVars &vars);
protected:
  int experiment;
  string experiment_name;
  PetscErrorCode init_till_phi();
  PetscScalar phiLocal(const PetscScalar width, 
         const PetscScalar x, const PetscScalar y,
         const PetscScalar STRONG, const PetscScalar UP, const PetscScalar DOWN);
  int sectorNumberP2(const PetscScalar x, const PetscScalar y);
};

#endif /* __icePSTexModel_hh */

