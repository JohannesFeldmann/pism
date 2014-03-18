// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir, Andy Aschwanden and Torsten Albrecht
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

#include "POConstantPIK.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

POConstantPIK::POConstantPIK(IceGrid &g, const NCConfigVariable &conf)
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
}

PetscErrorCode POConstantPIK::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (!config.get_flag("is_dry_simulation")) {
    ierr = verbPrintf(2, grid.com, "* Initializing the constant ocean model...\n"); CHKERRQ(ierr);
  }

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  mask_array = dynamic_cast<IceModelVec2S*>(vars.get("mask"));
  if (!mask_array) { SETERRQ(grid.com, 1, "ERROR: mask not available"); }

  // for symmetric (MISMIP-like) experiments: non-zero meltrates only at RHS of domain 
  ierr = PISMOptionsIsSet("-melt1side", "melt1side", melt1side_set); CHKERRQ(ierr);
  // field of melt rates reduced to have zero melt rates where thinnest ice (hopefully at ice shelf front)
  ierr = PISMOptionsIsSet("-offset_meltrate", "offset_meltrate", offset_meltrate_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-offset_meltrate_PIG_TG", "offset_meltrate_PIG_TG", offset_meltrate_PIG_TG_set); CHKERRQ(ierr);

  ierr = PISMOptionsIsSet("-use_waterTemp_PIG_TG", waterTemp_PIG_TG_set); CHKERRQ(ierr); 

  if (waterTemp_PIG_TG_set) {
    waterTemp_array.resize(3);
    // default values:
    waterTemp_array[0] = -1.7; 
    waterTemp_array[1] = -1.7;
    waterTemp_array[2] = -1.7;
    
    ierr = PISMOptionsRealArray("-use_waterTemp_PIG_TG", "waterTemp_PIGn, waterTemp_PIGs, waterTemp_TG",
                              waterTemp_array, waterTemp_PIG_TG_set); CHKERRQ(ierr);
                             
    ierr = verbPrintf(2, grid.com,
                      "* Water temperatures for PIGn, PIGs and TG are are set separately to\n"
                      "     waterTemp_PIG_n=%f K, waterTemp_PIG_s=%f K, waterTemp_TG=%f K \n", waterTemp_array[0], waterTemp_array[1], waterTemp_array[2]); CHKERRQ(ierr);
                              
    if (waterTemp_array.size() != 3) {
      PetscPrintf(grid.com,
                "PISM ERROR: option -use_waterTemp_PIG_TG requires a comma-separated list with 3 numbers; got %d\n",
                waterTemp_array.size());
      PISMEnd();
    }                      
  }

  ierr = PISMOptionsIsSet("-use_meltfactor_PIG_TG", meltfactor_PIG_TG_set); CHKERRQ(ierr); 

  if (meltfactor_PIG_TG_set) {
    meltfactor_array.resize(3);
    // default values:
    meltfactor_array[0] = -1.7; 
    meltfactor_array[1] = -1.7;
    meltfactor_array[2] = -1.7;
    
    ierr = PISMOptionsRealArray("-use_meltfactor_PIG_TG", "meltfactor_PIGn, meltfactor_PIGs, meltfactor_TG",
                              meltfactor_array, meltfactor_PIG_TG_set); CHKERRQ(ierr);
                             
    ierr = verbPrintf(2, grid.com,
                      "* Melt factors for PIGn, PIGs and TG are are set separately to\n"
                      "     meltfactor_PIG_n=%f, meltfactor_PIG_s=%f, meltfactor_TG=%f\n", meltfactor_array[0], meltfactor_array[1], meltfactor_array[2]); CHKERRQ(ierr);
                              
    if (meltfactor_array.size() != 3) {
      PetscPrintf(grid.com,
                "PISM ERROR: option -use_meltfactor_PIG_TG requires a comma-separated list with 3 numbers; got %d\n",
                meltfactor_array.size());
      PISMEnd();
    }                      
  }

  ierr = PISMOptionsIsSet("-use_boundary_PIG_TG", boundary_PIG_TG_set); CHKERRQ(ierr); 

  if (boundary_PIG_TG_set) {
    boundary_array.resize(2);
    // default values:
    boundary_array[0] = 50; 
    boundary_array[1] = 25;
    
    ierr = PISMOptionsRealArray("-use_boundary_PIG_TG", "boundary_innerPIG, boundary_TG",
                              boundary_array, boundary_PIG_TG_set); CHKERRQ(ierr);
                             
    ierr = verbPrintf(2, grid.com,
                      "* Boundaries for inner PIG, and TG ice shelves to\n"
                      "     boundary_innerPIG=%i, boundary_TG=%i \n", boundary_array[0], boundary_array[1]); CHKERRQ(ierr);
                              
    if (boundary_array.size() != 2) {
      PetscPrintf(grid.com,
                "PISM ERROR: option -use_boundary_PIG_TG requires a comma-separated list with 2 numbers; got %d\n",
                boundary_array.size());
      PISMEnd();
    }                      
  }

  return 0;
}

PetscErrorCode POConstantPIK::sea_level_elevation(PetscReal &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POConstantPIK::shelf_base_temperature(IceModelVec2S &result) {
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
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME issue #15
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes mass flux in ice-equivalent m s-1.
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
PetscErrorCode POConstantPIK::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  PetscReal L = config.get("water_latent_heat_fusion"),
    rho_ocean = config.get("sea_water_density"),
    rho_ice = config.get("ice_density");

  const PetscScalar c_p_ocean	  = 3974.0,   // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_T	  = 1e-4;     // m/s, thermal exchange velocity
  //FIXME: gamma_T should be a function of the friction velocity, not a const

  PetscScalar ocean_salinity = 35.0; 

  PetscScalar T_water = -1.7, //Default in PISM-PIK
    T_ocean = 273.15 + T_water;

  // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
  // PetscReal meltrate = config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * rho_ice); // m s-1

  PetscReal meltfactor = 5e-3;
  bool meltfactorSet;
  double meltfactor_pik;
  ierr = PISMOptionsReal("-meltfactor_pik",
                         "Uses as a meltfactor as in sub-shelf-melting parameterization of martin_winkelmann11",
                         meltfactor_pik, meltfactorSet); CHKERRQ(ierr);

  if (meltfactorSet) {
    meltfactor = meltfactor_pik; // default is 5e-3 as in martin_winkelmann11
  }

  if (waterTemp_PIG_TG_set) {
    T_water_PIGn = waterTemp_array[0];
    T_water_PIGs = waterTemp_array[1];
    T_water_TG = waterTemp_array[2];
    // nomelt_PIGn = int((53*5000)/dx + 0.5);
    // nomelt_TGw = int((58*5000)/dx + 0.5);
    // nomelt_TGs = int((86*5000)/dx + 0.5);
    // nomelt_PIGTGn = int((71*5000)/dx + 0.5);
    // nomelt_PIGTGe = int((93*5000)/dx + 0.5);
    // nomelt_PIGTGs = int((75*5000)/dx + 0.5);
    // nomelt_PIGTGw = int((78*5000)/dx + 0.5);
  }

  if (meltfactor_PIG_TG_set) {
    meltfactor_PIGn = meltfactor_array[0];
    meltfactor_PIGs = meltfactor_array[1];
    meltfactor_TG = meltfactor_array[2];
  }			

  if (waterTemp_PIG_TG_set || meltfactor_PIG_TG_set) {
    // use boundaries for PIG north and south (64) and TG (82) which were defined
    // for 5km and scale to used resolution
    dx = grid.dx;
    innerPIGbound = int((50*5000)/dx + 0.5); //addition of 0.5 is because C++ automatically rounds down
    // TGbound = int((35*5000)/dx + 0.5);
    TGbound = int((25*5000)/dx + 0.5);
  }

  if (boundary_PIG_TG_set) {
    innerPIGbound = int((boundary_array[0]*5000)/dx + 0.5);
    TGbound = int((boundary_array[1]*5000)/dx + 0.5);
  }			

  PetscReal add_constant_bmr = 0.0;
  ierr = PISMOptionsReal("-add_constant_bmr",
                         "add constant to bmr field",
                         add_constant_bmr, add_constant_bmr_set); CHKERRQ(ierr);

  PetscScalar **H, **vmask;
  ierr = ice_thickness->get_array(H); CHKERRQ(ierr);
  ierr = mask_array->get_array(vmask); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  PetscScalar result_min = 10000.0, result_min_PIGs = 10000.0, result_min_PIGn = 10000.0, result_min_TG = 10000.0; 
  // set to very high value to make sure that it will be replaced by smaller value later

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // compute T_f[i][j] according to beckmann_goosse03, which has the
      // meaning of the freezing temperature of the ocean water directly
      // under the shelf, (of salinity 35psu) [this is related to the
      // Pressure Melting Temperature, see beckmann_goosse03 eq. 2 for
      // details]
      PetscScalar shelfbaseelev = - (rho_ice / rho_ocean) * H[i][j],
        T_f= 273.15 + (0.0939 -0.057 * ocean_salinity + 7.64e-4 * shelfbaseelev);
      // add 273.15 to get it in Kelvin

      if(waterTemp_PIG_TG_set && j > TGbound){
	if(i <= innerPIGbound){
	  T_ocean = 273.15 + T_water_PIGn;
	  // T_ocean = 273.15 + waterTemp_array[0];
	  meltfactor = meltfactor_PIGn;
	}
	if(i > innerPIGbound){
	  T_ocean = 273.15 + T_water_PIGs;
	  // T_ocean = 273.15 + waterTemp_array[1];
	  meltfactor = meltfactor_PIGs;
	}
      }
      if(waterTemp_PIG_TG_set && j <= TGbound){
	T_ocean = 273.15 + T_water_TG;
	// T_ocean = 273.15 + waterTemp_array[2];
	meltfactor = meltfactor_TG;
      }

      // compute ocean_heat_flux according to beckmann_goosse03
      // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
      PetscScalar oceanheatflux = meltfactor * rho_ocean * c_p_ocean * gamma_T * (T_ocean - T_f); // in W/m^2
      // TODO: T_ocean -> field!

      // shelfbmassflux is positive if ice is freezing on; here it is always negative:
      // same sign as OceanHeatFlux... positive if massflux FROM ice TO ocean
      result(i,j) = oceanheatflux / (L * rho_ice); // m s-1

      if (offset_meltrate_set) {
	if (H[i][j] > 0.0 && vmask[i][j] == 3) {
	  result_min = PetscMin(result_min,result(i,j));
	}
      }

      if (offset_meltrate_PIG_TG_set && vmask[i][j] == 3 && j > TGbound) {
	if(i <= innerPIGbound) {
	  result_min_PIGn = PetscMax(0,PetscMin(result_min_PIGn,result(i,j)));
	}
	if (i > innerPIGbound) {
	  result_min_PIGs = PetscMax(0,PetscMin(result_min_PIGs,result(i,j)));
	}
      }
      if (offset_meltrate_PIG_TG_set && vmask[i][j] == 3 && j <= TGbound) {
	result_min_TG = PetscMax(0,PetscMin(result_min_TG,result(i,j)));
      }

    }
  }
  
  PetscScalar result_min_global, result_min_PIGn_global, result_min_PIGs_global, result_min_TG_global;
  ierr = PISMGlobalMin(&result_min, &result_min_global, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMin(&result_min_PIGn, &result_min_PIGn_global, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMin(&result_min_PIGs, &result_min_PIGs_global, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMin(&result_min_TG, &result_min_TG_global, grid.com); CHKERRQ(ierr);

  // ierr = PetscPrintf(PETSC_COMM_SELF,"!!! PISM_INFO: result_min_PIGn=%12.2f m/a\n",result_min_PIGn_global*secpera); CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_SELF,"!!! PISM_INFO: result_min_PIGs=%12.2f m/a\n",result_min_PIGs_global*secpera); CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_SELF,"!!! PISM_INFO: result_min_TG=%12.2f m/a\n",result_min_TG_global*secpera); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (offset_meltrate_set) {
	if (H[i][j] > 0.0) {
	  result(i,j) = result(i,j) - result_min_global;
	}
      }

      if (offset_meltrate_PIG_TG_set && j > TGbound) {
	if(i <= innerPIGbound) {
	  result(i,j) = result(i,j) - result_min_PIGn_global;
	}
	if(i > innerPIGbound) {
	  result(i,j) = result(i,j) - result_min_PIGs_global;
	}
      }								
      if (offset_meltrate_PIG_TG_set && j <= TGbound){
	result(i,j) = result(i,j) - result_min_TG_global;
      }
      
      if (add_constant_bmr_set) {
	result(i,j) = result(i,j) + add_constant_bmr / secpera;
      }

      if (H[i][j] < 0.01) {
	result(i,j) = 0.0;
      }

      if (melt1side_set) {
      	if (i < ((grid.Mx - 1)/2)) {
      	  result(i,j) = 0.0;
      	}
      }

    }
  }

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = mask_array->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

void POConstantPIK::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (keyword == "medium" || keyword == "big") {
    result["shelfbtemp"] = shelfbtemp;
    result["shelfbmassflux"] = shelfbmassflux;
  }
}

PetscErrorCode POConstantPIK::define_variables(set<string> vars, const PIO &nc,
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

PetscErrorCode POConstantPIK::write_variables(set<string> vars, string filename) {
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
