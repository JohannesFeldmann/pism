// Copyright (C) 2008-2015 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir, Andy Aschwanden and Torsten Albrecht
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include <gsl/gsl_math.h>

#include "POMismip.hh"
#include "base/util/PISMVars.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/IceGrid.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace ocean {

Mismip::Mismip(IceGrid::ConstPtr g)
  : OceanModel(g),
    m_shelfbmassflux(m_sys, "shelfbmassflux"),
    m_shelfbtemp(m_sys, "shelfbtemp")
{
  m_shelfbmassflux.set_string("pism_intent", "climate_state");
  m_shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  m_shelfbmassflux.set_string("units", "kg m-2 s-1");
  m_shelfbmassflux.set_string("glaciological_units", "kg m-2 year-1");

  m_shelfbtemp.set_string("pism_intent", "climate_state");
  m_shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  m_shelfbtemp.set_string("units", "Kelvin");

  // m_meltfactor = m_config->get_double("ocean_pik_melt_factor");
}

Mismip::~Mismip() {
  // empty
}

void Mismip::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the MISMIP+ sub-shelf-melting parameterization...\n");
  
  // set defaults of mismip ocean parameters
  Omega = 0.2; // yr-1
  H_0   = 75.0; // m
  z_0   = -100.0; // m
  ice1  = false;

  options::RealList MOP("-mismip_ocean_parameters",
			"mismip sub-shelf melt parameterization");

  if (MOP.is_set()) {
    if (MOP->size() != 3) {
      throw RuntimeError("option -mismip_ocean_parameters requires an argument"
    			 " (comma-separated list of 3 numbers)");
    }
    Omega = MOP[0];
    H_0   = MOP[1];
    z_0   = MOP[2];
    ice1  = true;

    m_log->message(2,
		   "   Omega = %3.3f yr-1\n"
		   "   H_0 =   %3.3f m\n"
		   "   z_0 =   %3.3f m\n",
		   Omega, H_0, z_0);
  }

  // set defaults of mismip sub-shelf melt experiment Ice2
  meltr = 100.0; // yr-1
  bound = 480000.0; // m
  ice2  = false;

  options::RealList ICE2("-mismip_ocean_ice2",
			"mismip sub-shelf melt experiment Ice2");

  if (ICE2.is_set()) {
    if (ICE2->size() != 2) {
      throw RuntimeError("option -mismip_ocean_ice2 requires an argument"
    			 " (comma-separated list of 2 numbers)");
    }
    meltr = ICE2[0];
    bound = ICE2[1];
    ice2  = true;

    m_log->message(2, 
		   "   \nApplying sub-shelf melting according to MISMIP+ experiment Ice2...\n"
		   "   meltrate = %3.3f m yr-1\n"
		   "   boundary =   %3.3f m\n",
		   meltr, bound);
  }

}

MaxTimestep Mismip::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Mismip::update_impl(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
}

void Mismip::sea_level_elevation_impl(double &result) {
  result = m_sea_level;
}

void Mismip::shelf_base_temperature_impl(IceModelVec2S &result) {
  const double
    T0          = m_config->get_double("water_melting_point_temperature"), // K
    beta_CC     = m_config->get_double("beta_CC"),
    g           = m_config->get_double("standard_gravity"),
    ice_density = m_config->get_double("ice_density");

  const IceModelVec2S &H = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(H);
  list.add(result);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * H(i,j); // FIXME task #7297
    // temp is set to melting point at depth
    result(i,j) = T0 - beta_CC * pressure;
  }
}

//! \brief Computes mass flux in [kg m-2 s-1].
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
void Mismip::shelf_base_mass_flux_impl(IceModelVec2S &result) {
  const double
    sea_water_density = m_config->get_double("sea_water_density"),
    ice_density       = m_config->get_double("ice_density"),
    secpera           = units::convert(m_sys, 1.0, "year", "seconds");

    // L                 = m_config->get_double("water_latent_heat_fusion"),
    // sea_water_density = m_config->get_double("sea_water_density"),
    // ice_density       = m_config->get_double("ice_density"),
    // c_p_ocean         = 3974.0, // J/(K*kg), specific heat capacity of ocean mixed layer
    // gamma_T           = 1e-4,   // m/s, thermal exchange velocity
    // ocean_salinity    = 35.0,   // g/kg
    // T_ocean           = units::convert(m_sys, -1.7, "Celsius", "Kelvin");   //Default in PISM-PIK
  
  //FIXME: gamma_T should be a function of the friction velocity, not a const

  if (ice1 == true) {
    const IceModelVec2S &H = *m_grid->variables().get_2d_scalar("land_ice_thickness");
    const IceModelVec2S &z_base = *m_grid->variables().get_2d_scalar("bedrock_altitude");
    // const IceModelVec2S &z_base = *m_grid->variables().get_2d_scalar("bedrock_surface_elevation");

    IceModelVec::AccessList list;
    list.add(H);
    list.add(z_base);
    list.add(result);
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // compute base elevation of ice shelf (ice draft)
      double
	z_bot = - (ice_density / sea_water_density) * H(i,j);

      // compute thickness of ice-shelf cavity
      double
	H_cav = z_bot - z_base(i,j);

      result(i,j) = Omega * tanh( H_cav / H_0 ) * std::max( ( z_0 - z_bot ), 0.0); // m yr-1

      // convert from [m yr-1] to [kg m-2 s-1]:
      result(i,j) *= ice_density / secpera;

    }
  }

  if (ice2 == true) {
    IceModelVec::AccessList list;
    list.add(result);

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      {
	double x = m_grid->x(i);
	if (x < -bound || x > bound) {
	  result(i,j) = meltr * ice_density / secpera;
	}
	else {
	  result(i,j) = 0.0;
	}
      }
    }
  }

  if (ice1 == true && ice2 == true) {
    // m_log->message(2,
    //            "PISM WARNING: options -mismip_ocean_parameters AND"
    // 		   " -mismip_ocean_ice2 set");
    throw RuntimeError("PISM ERROR: options -mismip_ocean_parameters AND"
		       " -mismip_ocean_ice2 set");
  }

    // add 273.15 to convert from Celsius to Kelvin

    // compute ocean_heat_flux according to beckmann_goosse03
    // positive, if T_oc > T_ice ==> heat flux FROM ocean TO ice
    // double ocean_heat_flux = m_meltfactor * sea_water_density * c_p_ocean * gamma_T * (T_ocean - T_f); // in W/m^2

    // TODO: T_ocean -> field!

    // shelfbmassflux is positive if ice is freezing on; here it is always negative:
    // same sign as ocean_heat_flux (positive if massflux FROM ice TO ocean)
    // result(i,j) = ocean_heat_flux / (L * ice_density); // m s-1

    // convert from [m s-1] to [kg m-2 s-1]:
    // result(i,j) *= ice_density;
}

void Mismip::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert("shelfbtemp");
    result.insert("shelfbmassflux");
  }
}

void Mismip::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype) {
  std::string order = m_grid->ctx()->config()->get_string("output_variable_order");

  if (set_contains(vars, "shelfbtemp")) {
    io::define_spatial_variable(m_shelfbtemp, *m_grid, nc, nctype, order, true);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    io::define_spatial_variable(m_shelfbmassflux, *m_grid, nc, nctype, order, true);
  }
}

void Mismip::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (not tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = m_shelfbtemp;
    shelf_base_temperature(tmp);
    tmp.write(nc);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);
    }

    tmp.metadata() = m_shelfbmassflux;
    tmp.write_in_glaciological_units = true;
    shelf_base_mass_flux(tmp);
    tmp.write(nc);
  }
}

} // end of namespace ocean
} // end of namespace pism
