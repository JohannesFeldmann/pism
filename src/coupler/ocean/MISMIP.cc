// Copyright (C) 2008-2019 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "MISMIP.hh"
#include "pism/util/Vars.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/Geometry.hh"

#include "pism/util/pism_options.hh"

namespace pism {
namespace ocean {

MISMIP::MISMIP(IceGrid::ConstPtr g)
  : CompleteOceanModel(g) {
  // empty
}

MISMIP::~MISMIP() {
  // empty
}

void MISMIP::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
		 "* Initializing the MISMIP+ sub-shelf-melting parameterization...\n");
  {
    {
      // set defaults of mismip ocean parameters
      Omega = 0.2; // yr-1
      H_0   = 75.0; // m
      z_0   = -100.0; // m
      ice1  = false;

      options::RealList MOP("-mismip_ocean_parameters", "mismip sub-shelf melt parameterization", {Omega,H_0,z_0});

      if (MOP.is_set()) {
	if (MOP->size() != 3) {
	  throw RuntimeError(PISM_ERROR_LOCATION, "option -mismip_ocean_parameters requires an argument"
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
    }

    {
      // set defaults of mismip sub-shelf melt experiment Ice2
      meltr = 100.0; // yr-1
      bound = 480000.0; // m
      ice2  = false;

      options::RealList ICE2("-mismip_ocean_ice2", "mismip sub-shelf melt experiment Ice2", {meltr,bound});

      if (ICE2.is_set()) {
	if (ICE2->size() != 2) {
	  throw RuntimeError(PISM_ERROR_LOCATION, "option -mismip_ocean_ice2 requires an argument"
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
  }
}
 

MaxTimestep MISMIP::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean MISMIP");
}

void MISMIP::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  const IceModelVec2S &H = geometry.ice_thickness;

  // Set shelf base temperature to the melting temperature at the base (depth within the
  // ice equal to ice thickness).
  melting_point_temperature(H, *m_shelf_base_temperature);

  mass_flux(H, *m_shelf_base_mass_flux);
  
  m_melange_back_pressure_fraction->set(m_config->get_number("ocean.melange_back_pressure_fraction"));
}

/*!
 * Compute melting temperature at a given depth within the ice.
 */
void MISMIP::melting_point_temperature(const IceModelVec2S &depth,
                                    IceModelVec2S &result) const {
  const double
    T0          = m_config->get_number("constants.fresh_water.melting_point_temperature"), // K
    beta_CC     = m_config->get_number("constants.ice.beta_Clausius_Clapeyron"),
    g           = m_config->get_number("constants.standard_gravity"),
    ice_density = m_config->get_number("constants.ice.density");

  IceModelVec::AccessList list{&depth, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * depth(i,j); // FIXME task #7297
    // result is set to melting point at depth
    result(i,j) = T0 - beta_CC * pressure;
  }
}

//! \brief Computes mass flux in [kg m-2 s-1].
/*!
 * Assumes that mass flux is proportional to the shelf-base heat flux.
 */
void MISMIP::mass_flux(const IceModelVec2S &ice_thickness, IceModelVec2S &result) const {
  const double

    sea_water_density = m_config->get_number("constants.sea_water.density"),
    ice_density       = m_config->get_number("constants.ice.density"),
    secpera           = units::convert(m_sys, 1.0, "year", "seconds");

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
    throw RuntimeError(PISM_ERROR_LOCATION, "PISM ERROR: options -mismip_ocean_parameters AND"
		       " -mismip_ocean_ice2 set");
  }
  
}

} // end of namespace ocean
} // end of namespace pism
