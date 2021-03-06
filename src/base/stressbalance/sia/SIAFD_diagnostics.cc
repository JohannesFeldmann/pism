// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev
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

#include "SIAFD_diagnostics.hh"
#include "PISMBedSmoother.hh"
#include "base/util/PISMVars.hh"

namespace pism {
namespace stressbalance {

void SIAFD::get_diagnostics_impl(std::map<std::string, Diagnostic*> &dict,
                            std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  dict["diffusivity"] = new SIAFD_diffusivity(this);
  dict["diffusivity_staggered"] = new SIAFD_diffusivity_staggered(this);
  dict["schoofs_theta"] = new SIAFD_schoofs_theta(this);
  dict["thksmooth"] = new SIAFD_thksmooth(this);
  dict["topgsmooth"] = new SIAFD_topgsmooth(this);
  dict["h_x"] = new SIAFD_h_x(this);
  dict["h_y"] = new SIAFD_h_y(this);
}

SIAFD_schoofs_theta::SIAFD_schoofs_theta(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "schoofs_theta"));

  set_attrs("multiplier 'theta' in Schoof's (2003) theory of bed roughness in SIA", "",
            "1", "", 0);
  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("valid_max", 1);
}

IceModelVec::Ptr SIAFD_schoofs_theta::compute() {
  const IceModelVec2S *surface = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "schoofs_theta", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->m_bed_smoother->get_theta(*surface, *result);

  return result;
}


SIAFD_topgsmooth::SIAFD_topgsmooth(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "topgsmooth"));
  set_attrs("smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

IceModelVec::Ptr SIAFD_topgsmooth::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "topgsmooth", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  result->copy_from(model->m_bed_smoother->get_smoothed_bed());

  return result;
}

SIAFD_thksmooth::SIAFD_thksmooth(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "thksmooth"));
  set_attrs("thickness relative to smoothed bed elevation in Schoof's (2003) theory of bed roughness in SIA",
            "", "m", "m", 0);
}

IceModelVec::Ptr SIAFD_thksmooth::compute() {
  const IceModelVec2S *surface, *thickness;
  const IceModelVec2Int *mask;

  surface   = m_grid->variables().get_2d_scalar("surface_altitude");
  thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  mask      = m_grid->variables().get_2d_mask("mask");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "thksmooth", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  model->m_bed_smoother->get_smoothed_thk(*surface, *thickness, *mask,
                                          *result);
  return result;
}



SIAFD_diffusivity::SIAFD_diffusivity(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "diffusivity"));

  set_attrs("diffusivity of SIA mass continuity equation", "",
            "m2 s-1", "m2 s-1", 0);
}

IceModelVec::Ptr SIAFD_diffusivity::compute() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "diffusivity", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->compute_diffusivity(*result);

  return result;
}

SIAFD_diffusivity_staggered::SIAFD_diffusivity_staggered(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(SpatialVariableMetadata(m_sys, "diffusivity_i"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "diffusivity_j"));

  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (i-offset)", "",
            "m2 s-1", "m2 s-1", 0);
  set_attrs("diffusivity of SIA mass continuity equation on the staggered grid (j-offset)", "",
            "m2 s-1", "m2 s-1", 1);
}

IceModelVec::Ptr SIAFD_diffusivity_staggered::compute() {
  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "diffusivity", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->write_in_glaciological_units = true;

  model->compute_diffusivity_staggered(*result);

  return result;
}

SIAFD_h_x::SIAFD_h_x(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(SpatialVariableMetadata(m_sys, "h_x_i"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "h_x_j"));

  set_attrs("the x-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the x-component of the surface gradient, j-offset", "",
            "", "", 1);
}

IceModelVec::Ptr SIAFD_h_x::compute() {

  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "h_x", WITH_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->write_in_glaciological_units = true;

  model->compute_surface_gradient(model->m_work_2d_stag[0],
                                  model->m_work_2d_stag[1]);

  result->copy_from(model->m_work_2d_stag[0]);

  return result;
}

SIAFD_h_y::SIAFD_h_y(SIAFD *m)
  : Diag<SIAFD>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(SpatialVariableMetadata(m_sys, "h_y_i"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "h_y_j"));

  set_attrs("the y-component of the surface gradient, i-offset", "",
            "", "", 0);
  set_attrs("the y-component of the surface gradient, j-offset", "",
            "", "", 1);
}

IceModelVec::Ptr SIAFD_h_y::compute() {

  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "h_y", WITH_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->write_in_glaciological_units = true;

  model->compute_surface_gradient(model->m_work_2d_stag[0],
                                  model->m_work_2d_stag[1]);

  result->copy_from(model->m_work_2d_stag[1]);

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
