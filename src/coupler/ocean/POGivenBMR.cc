// Copyright (C) 2011, 2012 Constantine Khroulev
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

#include "POGivenBMR.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "pism_options.hh"

PetscErrorCode POGivenBMR::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool regrid = true;
  int start = -1;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean model 'BMR' (which reads 'shelfbmassflux' from file) ...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters(""); CHKERRQ(ierr);

  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  ierr = mass_flux.set_attrs("climate_forcing",
                       "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                       "m s-1", ""); CHKERRQ(ierr);

  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  // read time-independent data right away:
  if (mass_flux.get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  // adjustment of basal melt rates
  ierr = PISMOptionsIsSet("-adjust_bmr", adjust_bmr_set); CHKERRQ(ierr);

  if (adjust_bmr_set) {
    ierr = verbPrintf(2, grid.com,
                      "* Sub-shelf mass flux will be adjusted according to reference ice shelf base elevation"); CHKERRQ(ierr);

    ierr = find_pism_input(filename, regrid, start); CHKERRQ(ierr);

    ierr = melt_ref_thk.create(grid, "melt_ref_thk", false); CHKERRQ(ierr);
    ierr = melt_ref_thk.set_attrs("model_state",
            "reference ice geometry",
            "m",
            ""); CHKERRQ(ierr); // no CF standard_name ??
    // ierr = variables.add(melt_ref_thk); CHKERRQ(ierr); Would have to be done in "iceModel.cc"

    ref_openocean_shelfthk  = config.get("ref_openocean_shelfthk");
    meltrate_increase_per_K = config.get("meltrate_increase_per_K"); // m/(year*K)

    ierr = PISMOptionsReal("-adjust_bmr", "ref_openocean_shelfthk",
         ref_openocean_shelfthk, adjust_bmr_set); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-bmr_per_K", "meltrate_increase_per_K",
         meltrate_increase_per_K, bmr_per_K_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
                      "\n  - Reference ice thickness for initially open ocean area: %f\n  - Meltrate increase per K: %f\n",
          ref_openocean_shelfthk, meltrate_increase_per_K); CHKERRQ(ierr);

    // read reference ice geometry from file
    ierr = verbPrintf(2, grid.com,
          "  - Reading reference ice geometry ('melt_ref_thk') from '%s' ... \n",
          filename.c_str()); CHKERRQ(ierr);
    if (regrid) {
      ierr = melt_ref_thk.regrid(filename.c_str(), true); CHKERRQ(ierr); // fails if not found!
    } else {
      ierr = melt_ref_thk.read(filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
    }
    string ref_shelfbaseelev_history = "read from " + filename + "\n";

    ierr = melt_ref_thk.set_attr("history", ref_shelfbaseelev_history); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode POGivenBMR::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = mass_flux.at_time(t); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGivenBMR::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar T0 = config.get("water_melting_point_temperature"), // K
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * (*ice_thickness)(i,j); // FIXME issue #15
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access();   CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGivenBMR::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  PetscReal dbmrdz;
  PetscReal shelfbaseelev, ref_shelfbaseelev, reference_thickness;
  // convert input reference thk to shelfbaseelev
  //PetscReal ref_openocean_shelfbaseelev = - ( ice_rho / sea_water_rho ) * ref_openocean_shelfthk;

  ierr = mass_flux.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  if (adjust_bmr_set) {

    PetscScalar **H;
    ierr = melt_ref_thk.begin_access(); CHKERRQ(ierr);
    ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        reference_thickness = melt_ref_thk(i,j);

        //ierr = verbPrintf(3, grid.com,
          //"reference thickness orig = %f\n", reference_thickness); CHKERRQ(ierr);

        if ( reference_thickness == 0)
          reference_thickness = ref_openocean_shelfthk;

        //ierr = verbPrintf(3, grid.com,
        //  "reference thickness adapted = %f\n", reference_thickness); CHKERRQ(ierr);

        shelfbaseelev     = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME issue #15
        ref_shelfbaseelev = - ( ice_rho / sea_water_rho ) * reference_thickness;

        // first order correction to melt rate with
        // bmr(z) = bmr0(z) + dbmr/dz * (z-z0)
        // db/dz is a function of bmr predominantely, we use an exponential fit here
        // parameters for yearly melt rates
        //dbmrdz = 0.0321799 *( 1- 0.82526363*exp(0.02508303*mass_flux(i,j)*secpera));
        dbmrdz = -0.03337955 + 0.02736375*exp(-0.02269549*mass_flux(i,j)*secpera);
        //array([-0.03337955,  0.02736375,  0.02269549])
        //dT_pmp = beta_CC_grad * ( ref_shelfbaseelev - shelfbaseelev );
        if (mass_flux(i,j) != 0)
          ierr = verbPrintf(3, grid.com, "b0, dbdz = %e, %e\n", mass_flux(i,j)*secpera, dbmrdz ); CHKERRQ(ierr);
        //result(i,j) = mass_flux(i,j) + ( meltrate_increase_per_K * dT_pmp ) / secpera;
        result(i,j) = mass_flux(i,j) + dbmrdz/secpera * (shelfbaseelev - ref_shelfbaseelev) ;

        //ierr = verbPrintf(3, grid.com, "dbmrdz = %e\n", dbmrdz); CHKERRQ(ierr);
        //ierr = verbPrintf(3, grid.com, "adjusted bmelt, delta bmelt = %e, %e\n",
        //                  result(i,j), dbmrdz * ( ref_shelfbaseelev - shelfbaseelev ) / secpera); CHKERRQ(ierr);
        //ierr = verbPrintf(3, grid.com, "z0, z, dz = %e, %e, %e\n",
        //                  ref_shelfbaseelev, shelfbaseelev, ref_shelfbaseelev - shelfbaseelev); CHKERRQ(ierr);

       }
    }

    ierr = ice_thickness->end_access(); CHKERRQ(ierr);
    ierr = melt_ref_thk.end_access(); CHKERRQ(ierr);


  } else {

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        result(i,j) = mass_flux(i,j);

      }
    }
  }

  ierr = mass_flux.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

void POGivenBMR::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (keyword == "medium" || keyword == "big") {
  result["shelfbtemp"] = shelfbtemp;
  result["shelfbmassflux"] = shelfbmassflux;
  }
}

PetscErrorCode POGivenBMR::define_variables(set<string> vars, const PIO &nc,
                                            PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.define(nc, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode POGivenBMR::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (adjust_bmr_set) {
    ierr = melt_ref_thk.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}
