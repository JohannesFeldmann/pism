// Copyright (C) 2011, 2012 PISM Authors
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

#include "POConstant.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

POConstant::POConstant(IceGrid &g, const NCConfigVariable &conf)
  : PISMOceanModel(g, conf) {

  shelfbmassflux.init_2d("shelfbmassflux", g);
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("m s-1");
  shelfbmassflux.set_glaciological_units("m year-1");

  shelfbtemp.init_2d("shelfbtemp", g);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");

  mymeltrate = 0.0;
  meltrate_set = false;
}

PetscErrorCode POConstant::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (!config.get_flag("is_dry_simulation")) {
    ierr = verbPrintf(2, grid.com, "* Initializing the constant ocean model...\n"); CHKERRQ(ierr);
  }

  ierr = PetscOptionsBegin(grid.com, "", "Ocean model", ""); CHKERRQ(ierr);

  ierr = PISMOptionsReal("-shelf_base_melt_rate",
                          "Specifies a sub shelf ice-equivalent melt rate in meters/year",
			  mymeltrate, meltrate_set); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (meltrate_set) {
    ierr = verbPrintf(2, grid.com,
                      "    - option '-shelf_base_melt_rate' seen, "
                      "setting basal sub shelf basal melt rate to %.2f m/year ... \n",
                      mymeltrate); CHKERRQ(ierr);
  }

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }
  //inserted
  bed_topography = dynamic_cast<IceModelVec2S*>(vars.get("topg"));
  if (!bed_topography) { SETERRQ(grid.com, 1, "ERROR: bed topography is not available"); }
  mask_array = dynamic_cast<IceModelVec2S*>(vars.get("mask"));
  if (!mask_array) { SETERRQ(grid.com, 1, "ERROR: mask not available"); }
  gl_mask_array = dynamic_cast<IceModelVec2S*>(vars.get("gl_mask"));
  if (!gl_mask_array) { SETERRQ(grid.com, 1, "ERROR: gl_mask not available"); }
  //inserted
  return 0;
}

PetscErrorCode POConstant::sea_level_elevation(PetscReal &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POConstant::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar T0 = config.get("water_melting_point_temperature"), // K
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME task #7297
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  
  return 0;                                 
}

//! Computes mass flux in ice-equivalent m s-1, from assumption that basal heat flux rate converts to mass flux.
PetscErrorCode POConstant::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal L = config.get("water_latent_heat_fusion"),
    rho = config.get("ice_density"),
    meltrate;
  
  //inserted
  scalearray.resize(2);
  ierr = PISMOptionsRealArray("-scale_bmr_gl", "bmr_gl_fact, topg_thresh",
  				scalearray, scale_bmr_gl_set); CHKERRQ(ierr);

  PetscScalar **vbed;
  PetscScalar **vmask;
  PetscScalar **vgl_mask;
  PetscReal bmr_gl_fact = scalearray[0], topg_thresh = scalearray[1];

  const bool sub_gl = config.get_flag("sub_groundingline");

  ierr = PISMOptionsIsSet("-gl_strip", "gl_strip", gl_strip_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-gl_strip_large", "gl_strip_large", gl_strip_large_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-scale_subgl", "scale_subgl", scale_subgl_set); CHKERRQ(ierr);

  ierr = bed_topography->get_array(vbed); CHKERRQ(ierr);
  ierr = mask_array->get_array(vmask); CHKERRQ(ierr);
  ierr = gl_mask_array->get_array(vgl_mask); CHKERRQ(ierr);
  // inserted

  ierr = result.begin_access(); CHKERRQ(ierr);

  if (meltrate_set) {

    meltrate = convert(mymeltrate,"m year-1","m s-1");
    // meltrate = convert(mymeltrate,"m year-1","m s-1");
    ierr = result.set(meltrate); CHKERRQ(ierr);

  } 
  if (scale_bmr_gl_set) {
    // ierr = verbPrintf(2, grid.com, "scale_subgl!!! \n"); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

	if (gl_strip_set) {
	  if (vmask[i][j] == 3 && (vmask[i+1][j] == 2 || vmask[i-1][j] == 2 || 
				   vmask[i][j-1] == 2 || vmask[i][j-1] == 2 ||
				   vmask[i+1][j+1] == 2 || vmask[i+1][j-1] == 2 || 
				   vmask[i-1][j+1] == 2 || vmask[i-1][j-1] == 2)) {
	    // scale bmr at thin strip of 1 grid cell width along grounding line where topg is beneath threshold
	    if (vbed[i][j] < topg_thresh) {
	      result(i,j) = bmr_gl_fact * result(i,j);	    
	    }
	  }
	}
	else if (gl_strip_large_set) {
	  if (vmask[i][j] == 3 && (vmask[i+1][j] == 2 || vmask[i-1][j] == 2 || 
				   vmask[i][j-1] == 2 || vmask[i][j-1] == 2 ||
				   vmask[i+1][j+1] == 2 || vmask[i+1][j-1] == 2 || 
				   vmask[i-1][j+1] == 2 || vmask[i-1][j-1] == 2 ||
				   vmask[i+2][j] == 2 || vmask[i-2][j] == 2 || 
				   vmask[i][j-2] == 2 || vmask[i][j-2] == 2)) {
	    // scale bmr at thin strip of 2 grid cells width along grounding line where topg is beneath threshold
	    if (vbed[i][j] < topg_thresh) {
	      result(i,j) = bmr_gl_fact * result(i,j);	    
	    }
	  }
	}
	if (scale_subgl_set && vmask[i][j] == 2 && vgl_mask[i][j] < 1 && vgl_mask[i][j] > 0) {
	  // scale bmr at grounded but partly floating ice (i.e. where gl is interpolated) 
	  // where topg is beneath threshold
	  if (vbed[i][j] < topg_thresh) {
	    result(i,j) = bmr_gl_fact + result(i,j);	    
	    // result(i,j) = bmr_gl_fact * result(i,j);	    
	    // ierr = verbPrintf(2, grid.com, "scale_subgl!!! \n"); CHKERRQ(ierr);
	  }
	}
      }
    }
  }

    // commented
    // // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
    // meltrate = config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * rho); // m s-1
    // ierr = result.set(meltrate); CHKERRQ(ierr);
    // commented
  ierr = bed_topography->end_access(); CHKERRQ(ierr);
  ierr = mask_array->end_access(); CHKERRQ(ierr);
  ierr = gl_mask_array->end_access(); CHKERRQ(ierr);

  return 0;
}

void POConstant::add_vars_to_output(string, map<string,NCSpatialVariable> &result) {
  result["shelfbtemp"] = shelfbtemp;
  result["shelfbmassflux"] = shelfbmassflux;
}

PetscErrorCode POConstant::define_variables(set<string> vars, const PIO &nc,
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

PetscErrorCode POConstant::write_variables(set<string> vars, string filename) {
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

  return 0;
}
