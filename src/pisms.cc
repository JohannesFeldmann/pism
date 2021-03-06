// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "Ice sheet driver for EISMINT II, and other constant climate, simplified geometry\n"
  "intercomparison simulations.\n";

#include <petscsys.h>

#include "base/iceModel.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfig.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_options.hh"
#include "eismint/iceEISModel.hh"
#include "base/util/Context.hh"
#include "base/util/Logger.hh"
#include "base/util/PISMTime.hh"
#include "base/enthalpyConverter.hh"
#include "base/util/io/PIO.hh"

using namespace pism;

Context::Ptr pisms_context(MPI_Comm com) {
  // unit system
  units::System::Ptr sys(new units::System);

  // logger
  Logger::Ptr logger = logger_from_options(com);

  // configuration parameters
  Config::Ptr config = config_from_options(com, *logger, sys);

  config->set_string("calendar", "none");
  config->set_double("grid_Lx", 750e3);
  config->set_double("grid_Ly", 750e3);
  config->set_string("grid_periodicity", "none");
  config->set_string("sia_flow_law", "pb");

  set_config_from_options(*config);

  print_config(*logger, 3, *config);

  Time::Ptr time = time_from_options(com, config, sys);

  EnthalpyConverter::Ptr EC = enthalpy_converter_from_options(*config);

  return Context::Ptr(new Context(com, sys, config, EC, time, logger, "pisms"));
}

IceGrid::Ptr pisms_grid(Context::Ptr ctx) {
  options::String input_file("-i", "Specifies a PISM input file");
  options::forbidden("-bootstrap");

  if (input_file.is_set()) {
    Periodicity p = string_to_periodicity(ctx->config()->get_string("grid_periodicity"));

    // get grid from a PISM input file
    std::vector<std::string> names;
    names.push_back("enthalpy");
    names.push_back("temp");

    return IceGrid::FromFile(ctx, input_file, names, p);
  } else {
    // use defaults from the configuration database
    GridParameters P(ctx->config());
    P.horizontal_size_from_options();
    P.horizontal_extent_from_options();
    P.vertical_grid_from_options(ctx->config());
    P.ownership_ranges_from_options(ctx->size());

    return IceGrid::Ptr(new IceGrid(ctx, P));
  }
}


int main(int argc, char *argv[]) {

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  try {
    verbosityLevelFromOptions();

    verbPrintf(2,com, "PISMS %s (simplified geometry mode)\n",
               PISM_Revision);

    if (options::Bool("-version", "stop after printing print PISM version")) {
      return 0;
    }

    std::string usage =
      "  pisms [-eisII x] [OTHER PISM & PETSc OPTIONS]\n"
      "where major option chooses type of simplified experiment:\n"
      "  -eisII x    choose EISMINT II experiment (x = A|B|C|D|E|F|G|H|I|J|K|L)\n";

    std::vector<std::string> required;
    required.clear(); // no actually required options; "-eisII A" is default

    bool done = show_usage_check_req_opts(com, "pisms", required, usage);
    if (done) {
      return 0;
    }

    std::string experiment = options::Keyword("-eisII", "EISMINT II experiment name",
                                              "A,B,C,D,E,F,G,H,I,J,K,L", "A");

    Context::Ptr ctx = pisms_context(com);
    Config::Ptr config = ctx->config();

    IceGrid::Ptr g = pisms_grid(ctx);
    IceEISModel m(g, ctx, experiment[0]);

    m.init();

    m.run();

    verbPrintf(2,com, "... done with run \n");

    // provide a default output file name if no -o option is given.
    m.writeFiles("unnamed.nc");

    print_unused_parameters(*ctx->log(), 3, *config);
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
