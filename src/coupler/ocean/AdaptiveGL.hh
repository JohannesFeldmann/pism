// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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

#ifndef _POADAPTIVEGL_H_
#define _POADAPTIVEGL_H_

#include "CompleteOceanModel.hh"

namespace pism {
namespace ocean {
//! \brief Implements constant melting in a stripe close to the GL or along the boundaries.
class AdaptiveGL : public CompleteOceanModel {
public:
  AdaptiveGL(IceGrid::ConstPtr g);
  virtual ~AdaptiveGL();

private:
  MaxTimestep max_timestep_impl(double t) const;
  void update_impl(const Geometry &geometry, double my_t, double my_dt);
  void init_impl(const Geometry &geometry);

  void melting_point_temperature(const IceModelVec2S &depth, IceModelVec2S &result) const;
  void mass_flux(const IceModelVec2S &ice_thickness, IceModelVec2S &result);

  IceModelVec2Int melt_mask;

  double melt_rate,
         experiment,
         width,
         length,
         dist;
  static const int  maskfloating,
                    maskocean,
                    maskgrounded,
                    expmask_include,
                    expmask_neighboring,
                    expmask_unidentified;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POADAPTIVEGL_H_ */
