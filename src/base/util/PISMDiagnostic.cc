/* Copyright (C) 2015, 2016, 2017 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMDiagnostic.hh"
#include "base/util/PISMTime.hh"
#include "error_handling.hh"
#include "io/io_helpers.hh"
#include "base/util/Logger.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

Diagnostic::Ptr Diagnostic::wrap(const IceModelVec2S &input) {
  return Ptr(new DiagWithDedicatedStorage<IceModelVec2S>(input));
}

Diagnostic::Ptr Diagnostic::wrap(const IceModelVec2V &input) {
  return Ptr(new DiagWithDedicatedStorage<IceModelVec2V>(input));
}

Diagnostic::Diagnostic(IceGrid::ConstPtr g)
  : m_grid(g), m_sys(g->ctx()->unit_system()), m_config(g->ctx()->config()) {
  m_dof = 1;
  m_fill_value = m_config->get_double("output.fill_value");
}

Diagnostic::~Diagnostic() {
  // empty
}

//! \brief Update a cumulative quantity needed to compute a rate of change.
//! So far we there is only one such quantity: the rate of change of the ice
//! thickness.
void Diagnostic::update(double dt) {
  this->update_impl(dt);
}

void Diagnostic::update_impl(double dt) {
  (void) dt;
  // empty
}

void Diagnostic::reset() {
  this->reset_impl();
}

void Diagnostic::reset_impl() {
  // empty
}

//! Get the number of NetCDF variables corresponding to a diagnostic quantity.
unsigned int Diagnostic::n_variables() const {
  return m_dof;
}

void Diagnostic::init(const PIO &input, unsigned int time) {
  this->init_impl(input, time);
}

void Diagnostic::define_state(const PIO &output) const {
  this->define_state_impl(output);
}

void Diagnostic::write_state(const PIO &output) const {
  this->write_state_impl(output);
}

void Diagnostic::init_impl(const PIO &input, unsigned int time) {
  (void) input;
  (void) time;
  // empty
}

void Diagnostic::define_state_impl(const PIO &output) const {
  (void) output;
  // empty
}

void Diagnostic::write_state_impl(const PIO &output) const {
  (void) output;
  // empty
}

//! Get a metadata object corresponding to variable number N.
SpatialVariableMetadata& Diagnostic::metadata(unsigned int N) {
  if (N >= m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "variable metadata index %d is out of bounds",
                                  N);
  }

  return m_vars[N];
}

void Diagnostic::define(const PIO &file, IO_Type default_type) const {
  this->define_impl(file, default_type);
}

//! Define NetCDF variables corresponding to a diagnostic quantity.
void Diagnostic::define_impl(const PIO &file, IO_Type default_type) const {
  for (unsigned int j = 0; j < m_dof; ++j) {
    io::define_spatial_variable(m_vars[j], *m_grid, file,
                                default_type,
                                m_grid->ctx()->config()->get_string("output.variable_order"),
                                true);
  }
}

//! \brief A method for setting common variable attributes.
void Diagnostic::set_attrs(const std::string &long_name,
                           const std::string &standard_name,
                           const std::string &units,
                           const std::string &glaciological_units,
                           unsigned int N) {
  if (N >= m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "N (%d) >= m_dof (%d)", N, m_dof);
  }

  m_vars[N].set_string("pism_intent", "diagnostic");

  m_vars[N].set_string("long_name", long_name);

  m_vars[N].set_string("standard_name", standard_name);

  m_vars[N].set_string("units", units);

  if (not glaciological_units.empty()) {
    m_vars[N].set_string("glaciological_units", glaciological_units);
  }
}

IceModelVec::Ptr Diagnostic::compute() const {
  // use the name of the first variable
  std::vector<std::string> names;
  for (auto &v : m_vars) {
    names.push_back(v.get_name());
  }
  std::string all_names = join(names, ",");

  m_grid->ctx()->log()->message(3, "-  Computing %s...\n", all_names.c_str());
  IceModelVec::Ptr result = this->compute_impl();
  m_grid->ctx()->log()->message(3, "-  Done computing %s.\n", all_names.c_str());
  result->write_in_glaciological_units = true;
  return result;
}

TSDiagnostic::TSDiagnostic(IceGrid::ConstPtr g, const std::string &name, Kind kind)
  : m_kind(kind),
    m_grid(g),
    m_config(g->ctx()->config()),
    m_sys(g->ctx()->unit_system()),
    m_ts(*g, name, g->ctx()->config()->get_string("time.dimension_name")) {

  m_current_time = 0;
  m_accumulator  = 0.0;
  m_start        = 0;

  m_buffer_size = (size_t)m_config->get_double("output.timeseries.buffer_size");

  m_ts.dimension_metadata().set_string("calendar", m_grid->ctx()->time()->calendar());
  m_ts.dimension_metadata().set_string("long_name", m_config->get_string("time.dimension_name"));
  m_ts.dimension_metadata().set_string("axis", "T");
  m_ts.dimension_metadata().set_string("units", m_grid->ctx()->time()->CF_units_string());
}

TSDiagnostic::~TSDiagnostic() {
  flush();
}

void TSDiagnostic::evaluate_regular(double t0, double t1, double v0, double v1) {

  // skip times before the beginning of this time step
  while (m_current_time < m_times->size() and (*m_times)[m_current_time] < t0) {
    m_current_time += 1;
  }

  while (m_current_time < m_times->size() and (*m_times)[m_current_time] <= t1) {
    const unsigned int k = m_current_time;
    m_current_time += 1;

    // skip the first time: it defines the beginning of a reporting interval
    if (k == 0) {
      continue;
    }

    const double
      t_s = (*m_times)[k - 1],
      t_e = (*m_times)[k];

    const double v = v0 + (v1 - v0) * (t_e - t0) / (t1 - t0);
    m_ts.append(v, t_s, t_e);
  }
}

void TSDiagnostic::evaluate_rate(double t0, double t1, double change) {
  assert(t1 > t0);

  // skip times before the beginning of this time step
  while (m_current_time < m_times->size() and (*m_times)[m_current_time] < t0) {
    m_current_time += 1;
  }

  // number of requested times in this time step
  unsigned int N = 0;

  // loop through requested times that are within this time step
  while (m_current_time < m_times->size() and (*m_times)[m_current_time] <= t1) {
    const unsigned int k = m_current_time;
    m_current_time += 1;

    N += 1;

    // skip the first time: it defines the beginning of a reporting interval
    if (k == 0) {
      continue;
    }

    const double
      t_s = (*m_times)[k - 1],
      t_e = (*m_times)[k];

    double rate = 0.0;
    if (N == 1) {
      // it is the right end-point of the first reporting interval in this time step: count the
      // contribution from the last time step plus the one from the beginning of this time step
      const double
        total_change = m_accumulator + change * (t_e - t0) / (t1 - t0),
        dt           = t_e - t_s;

      rate = total_change / dt;

    } else {
      // this reporting interval is completely contained within the time step, so the rate of change
      // does not depend on its length
      rate = change / (t1 - t0);
    }

    m_ts.append(rate, t_s, t_e);
    m_accumulator = 0.0;
  }

  if (N == 0) {
    // if this time step contained no requested times we need to add the whole change to the
    // accumulator
    m_accumulator += change;
  } else {
    // if this time step contained some requested times we need to add the change since the last one
    // to the accumulator
    const double dt = t1 - (*m_times)[m_current_time - 1];
    m_accumulator += change * (dt / (t1 - t0));
  }
}

void TSDiagnostic::update(double t0, double t1) {
  double value = this->compute(t0, t1);

  assert(t1 > t0);

  m_v.push_back(value);

  if (m_v.size() == 2) {

    if (m_kind == SNAPSHOT) {
      evaluate_regular(t0, t1, m_v[0], m_v[1]);
    } else {
      evaluate_rate(t0, t1, value);
    }

    m_v.pop_front();
  }
}

void TSDiagnostic::define(const PIO &file) const {
  io::define_timeseries(m_ts.metadata(), file, PISM_DOUBLE);
  io::define_time_bounds(m_ts.bounds_metadata(), file, PISM_DOUBLE);
}


void TSDiagnostic::flush() {

  const std::vector<double> &times = m_ts.times();

  if (times.empty()) {
    return;
  }

  std::string dimension_name = m_ts.dimension_metadata().get_name();

  PIO file(m_grid->com, "netcdf3", m_output_filename, PISM_READWRITE); // OK to use netcdf3
  unsigned int len = file.inq_dimlen(dimension_name);

  if (len > 0) {
    double last_time;
    file.inq_dim_limits(dimension_name, NULL, &last_time);
    if (last_time < times.front()) {
      m_start = len;
    }
  }

  if (len == m_start) {
    io::write_timeseries(file, m_ts.dimension_metadata(), m_start, times);
    io::write_time_bounds(file, m_ts.bounds_metadata(), m_start, m_ts.bounds());
  }
  io::write_timeseries(file, m_ts.metadata(), m_start, m_ts.values());

  m_start += m_ts.times().size();

  m_ts.reset();
}

void TSDiagnostic::init(const std::string &output_filename,
                        std::shared_ptr<std::vector<double>> requested_times) {
  {
    // Get the number of records in the file (for appending):
    if (io::file_exists(m_grid->com, output_filename)) {
      PIO nc(m_grid->com, "netcdf3", output_filename, PISM_READONLY); // OK to use netcdf3
      m_start = nc.inq_dimlen(m_ts.dimension_metadata().get_name());
      // if (len > 0) {
      //   std::vector<double> tmp;

      //   if (nc.inq_var(name())) {
      //     nc.get_1d_var(name(), len - 1, 1, tmp);
      //   }
      // }
    } else {
      m_start = 0;
    }
  }

  m_output_filename = output_filename;
  m_times = requested_times;
}

const VariableMetadata &TSDiagnostic::metadata() const {
  return m_ts.metadata();
}

void TSDiagnostic::define_state(const PIO &file) const {
  // FIXME
}

void TSDiagnostic::write_state(const PIO &file) const {
  // FIXME
}

} // end of namespace pism
