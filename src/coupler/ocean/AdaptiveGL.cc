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

#include "AdaptiveGL.hh"
#include "pism/util/Vars.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace ocean {

AdaptiveGL::AdaptiveGL(IceGrid::ConstPtr g)
  : CompleteOceanModel(g)
  //   melt_mask(m_grid, "melt_mask", WITH_GHOSTS)
{
  melt_mask.create(m_grid, "melt_mask", WITH_GHOSTS);
  melt_mask.set_attrs("model_state",
  		      "mask displaying melting regions",
  		      "", "", "", 0);
  // empty
}

AdaptiveGL::~AdaptiveGL() {
  // empty
}

void AdaptiveGL::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the adaptive GL ocean model...\n");

  melt_rate = options::Real("-meltrate",
                             "Use as a constant melt rate",
                             melt_rate);
  melt_rate = melt_rate / 3.15569259747e7; // change unit to m per s

  experiment = options::Real("-experiment",
                             "Specify the experiment, 1 stands for central stripe, 2 for lateral stripes",
                             experiment);
  length = options::Real("-length",
                             "Specify the length of the applied melt stripe",
                             length);
  width = options::Real("-width",
                             "Specify the width of the applied melt stripe",
                             width);
  lathalf = options::Bool("-lathalf", "Lateral half");

  dist = 0.0;
  dist = options::Real("-dist",
		       "Specify the distance to the grounding line of the applied melt stripe",
		       dist);
  dist = static_cast<int>(dist);
}

MaxTimestep AdaptiveGL::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("ocean adeptive");
}

void AdaptiveGL::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  const IceModelVec2S &H = geometry.ice_thickness;

  // Set shelf base temperature to the melting temperature at the base (depth within the
  // ice equal to ice thickness).
  melting_point_temperature(H, *m_shelf_base_temperature);

  mass_flux(H, *m_shelf_base_mass_flux);

  m_melange_back_pressure_fraction->set(m_config->get_number("ocean.melange_back_pressure_fraction"));
}

const int AdaptiveGL::maskfloating = 3;
const int AdaptiveGL::maskocean    = 4;
const int AdaptiveGL::maskgrounded = 2;

const int AdaptiveGL::expmask_include = 2;  
const int AdaptiveGL::expmask_neighboring = 1;   
const int AdaptiveGL::expmask_unidentified= 0;  

/*!
 * Compute melting temperature at a given depth within the ice.
 */
void AdaptiveGL::melting_point_temperature(const IceModelVec2S &depth,
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
void AdaptiveGL::mass_flux(const IceModelVec2S &ice_thickness, IceModelVec2S &result) {
  const double
    melt_factor       = m_config->get_number("ocean.pik_melt_factor"),
    L                 = m_config->get_number("constants.fresh_water.latent_heat_of_fusion"),
    sea_water_density = m_config->get_number("constants.sea_water.density"),
    ice_density       = m_config->get_number("constants.ice.density"),
    c_p_ocean         = 3974.0, // J/(K*kg), specific heat capacity of ocean mixed layer
    gamma_T           = 1e-4,   // m/s, thermal exchange velocity
    ocean_salinity    = 35.0,   // g/kg
    T_ocean           = units::convert(m_sys, -1.7, "Celsius", "Kelvin"); //Default in PISM-PIK
 
  const IceModelVec2Int &mask = *m_grid->variables().get_2d_mask("mask");
  
  IceModelVec::AccessList list{&ice_thickness, &result};

  list.add(mask);
  list.add(melt_mask);

  melt_mask.set(0);
  melt_mask.update_ghosts();

  
  // identify the melt_mask, depending on the experimemt chosen
  if (experiment==1){ // melting at central part of GL
    m_log->message(2, "Experiment: Melting at the central part of GL, "); 
    m_log->message(2, "Melt rate is %f m per yr \n",melt_rate*3.15569259747e7); // in units m/yr
    m_log->message(2, "Dist is %f boxes \n",dist);
    int mid = (m_grid->My() - 1)/2;
    int midx = (m_grid->Mx() - 1)/2;

    // FIND FIRST ROW, STILL GROUNDED
    for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        if (j >= mid -(length-1)/2 && j <= mid + (length-1)/2){
	    // SHIFT BY "dist" AWAY FROM GROUNDING LINE INTO THE ICE SHELF
	    if (i < midx && mask(i+dist+1,j)==maskgrounded && mask(i+dist,j)==maskfloating) {
	      melt_mask(i,j) = expmask_include;
	    }
	    if (i > midx && mask(i-dist-1,j)==maskgrounded && mask(i-dist,j)==maskfloating) {
	      melt_mask(i,j) = expmask_include;
	    }

	  // if (mask(i,j)==maskgrounded){
	  //   // SHIFT BY "dist" AWAY FROM GROUNDING LINE INTO THE ICE SHELF
	  //   if (i < midx && mask(i-1,j)==maskfloating) {
	  //     melt_mask(i-dist-1,j) = expmask_include;
	  //   }
	  //   if (i > midx && mask(i+1,j)==maskfloating) {
	  //     melt_mask(i+dist+1,j) = expmask_include;
	  //   }
          // }
        } //if
    } //p

    melt_mask.update_ghosts();
    
    // ITERATE INTO SHELF
    for (int k=0; k<width-1; k++){
      m_log->message(2, "%d s iteration over width\n",k);
      //find possible cells
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();  

        if (mask(i,j)!=maskgrounded && melt_mask(i,j)==expmask_unidentified) {
	  if (i < midx && melt_mask(i+1,j) == expmask_include) {
	    melt_mask(i,j) = expmask_neighboring;
	  }
	  if (i > midx && melt_mask(i-1,j) == expmask_include) {
	    melt_mask(i,j) = expmask_neighboring;
	  }
	}

      } //p
      // relabel possible cells

      melt_mask.update_ghosts();

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();  
        if (melt_mask(i,j)==expmask_neighboring ){
          melt_mask(i,j) = expmask_include;
          m_log->message(2, "add label \n");
        }
      } //p

      melt_mask.update_ghosts();

    } // width

  } else if (experiment==2){ // melting along the sides
    m_log->message(2, "Experiment: Melting at lateral parts of GL, ");
    m_log->message(2, "Melt rate is %f m per yr \n",melt_rate*3.15569259747e7); // change units to m per a

    int mid_y = (m_grid->My() - 1)/2;
    int max_y = m_grid->My();
    int mid_left_x = 0;
    int mid_right_x = 0;
    int limit_left_x = 0;
    int limit_right_x = 0;

    // FIND x-coordinates:
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if(j == mid_y) {
          //m_log->message(2, "If j=mid_y.... ");
          if ( mask(i,j)==maskgrounded && mask(i-1,j)!=maskgrounded ){ // x-doordinate of mid_y in 'left' basin
            //m_log->message(2, "If at the left grounding line is true.... ");
            mid_left_x = i;
          }
          if ( mask(i,j)==maskgrounded && mask(i+1,j)!=maskgrounded){ // x-doordinate of mid_y in 'left' basin
            //m_log->message(2, "If at the right grounding line is true.... ");
            mid_right_x = i;
        }
      }
      if(j==0){ 
          if ( mask(i,j)==maskgrounded && mask(i-1,j)!=maskgrounded){ // x-doordinate of mid_y in 'left' basin
            //m_log->message(2, "If at the left grounding line is true.... ");
            limit_left_x = i;
          }
          if (mask(i,j)==maskgrounded && mask(i+1,j)!=maskgrounded){ // x-doordinate of mid_y in 'left' basin
            //m_log->message(2, "If at the right grounding line is true.... ");
            limit_right_x = i;
          }        
      }    
    }

    mid_left_x    = GlobalMax(m_grid->com, mid_left_x);
    mid_right_x   = GlobalMax(m_grid->com, mid_right_x);
    limit_left_x  = GlobalMax(m_grid->com, limit_left_x);
    limit_right_x = GlobalMax(m_grid->com, limit_right_x);


    m_log->message(2, "mid_y=%d, mid_left_x =%d, mid_right_x =%d, limit_left_x=%d, limit_right_x=%d\n ", mid_y, mid_left_x, mid_right_x, limit_left_x, limit_right_x);
    //FIXME still need to find appropriate algo to find the middle!
    int left_x = mid_left_x - static_cast<int>(round(2.0/5.0 * (mid_left_x - limit_left_x)));
    int right_x = mid_right_x + static_cast<int>(round(2.0/5.0 * (limit_right_x - mid_right_x)));
    
    m_log->message(2, "left_x=%d, right_x =%d\n ", left_x, right_x);
    
    melt_mask.update_ghosts();

    // find melt_mask, first row
    for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        if (j < (max_y-dist) && ((i >= left_x -length/2 && i <= left_x + length/2) || (i >= right_x - length/2 && i <=right_x + length/2 ))) {
	  m_log->message(2, "j = %i ", j);
	  // SHIFT BY "dist" AWAY FROM GROUNDING LINE INTO THE ICE SHELF
	  if (j < mid_y && mask(i,j-dist-1)==maskgrounded && mask(i,j-dist)==maskfloating && lathalf==false){
	    melt_mask(i,j) = expmask_include;
	  }
	  if (j > mid_y && mask(i,j+dist+1)==maskgrounded && mask(i,j+dist)==maskfloating){
	    melt_mask(i,j) = expmask_include;
	  }

          // if (mask(i,j)==maskgrounded && ( mask(i,j-1)==maskfloating || mask(i,j+1)==maskfloating) ){
	  //   // SHIFT BY "dist" AWAY FROM GROUNDING LINE INTO THE ICE SHELF
	  //   if (j < mid_y && lathalf==false){
	  //     melt_mask(i,j+dist+1) = expmask_include;
	  //   }
	  //   if (j > mid_y){
	  //     melt_mask(i,j-dist-1) = expmask_include;
	  //   }
          // }

	  
        } //if
    } //p

    melt_mask.update_ghosts();

    //ITERATE INTO SHELF
    for (int k=0; k<width-1; k++){
      //find possible cells
  	    m_log->message(2, "%d s iteration over width\n",k);
      
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        if (mask(i,j)!=maskgrounded && melt_mask(i,j)==expmask_unidentified) {
	  if (j < mid_y && melt_mask(i,j-1)==expmask_include) {
	    melt_mask(i,j) = expmask_neighboring;
	  }
	  if (j > mid_y && melt_mask(i,j+1)==expmask_include) {
	    melt_mask(i,j) = expmask_neighboring;
	  }
        }
      } //p
      
      melt_mask.update_ghosts();
    
      // relabel possible cells
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        if (melt_mask(i,j)==expmask_neighboring ){
          melt_mask(i,j) = expmask_include;
  	  m_log->message(2, "add label \n");
        }
      } //i
      
      melt_mask.update_ghosts();

    } // width

  } // experiment sides


  // set the melt rates 
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (melt_mask(i,j) == expmask_include && mask(i,j) == maskfloating ){
      result(i,j) = melt_rate;
    } else{
      result(i,j) = 0.0;
    }
    
    // convert from [m s-1] to [kg m-2 s-1]:
    result(i,j) *= ice_density;
  
  }
}


} // end of namespace ocean
} // end of namespace pism
