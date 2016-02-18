// Copyright (C) 2009--2016 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#include "base/util/IceGrid.hh"
#include "SSAFEM.hh"
#include "FETools.hh"
#include "base/util/Mask.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/rheology/FlowLaw.hh"
#include "base/util/pism_options.hh"
#include "base/util/error_handling.hh"
#include "base/util/PISMVars.hh"

#include "node_types.hh"

namespace pism {
namespace stressbalance {

/** The Q1 finite element SSA solver.
 *
 *
 *
 */
SSAFEM::SSAFEM(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e)
  : SSA(g, e),
    m_element_index(*g),
    m_element(*g),
    m_quadrature(g->dx(), g->dy(), 1.0),
    m_quadrature_vector(g->dx(), g->dy(), 1.0) {

  m_driving_stress_x = NULL;
  m_driving_stress_y = NULL;

  PetscErrorCode ierr;

  m_dirichletScale = 1.0;
  m_beta_ice_free_bedrock = m_config->get_double("beta_ice_free_bedrock");

  ierr = SNESCreate(m_grid->com, m_snes.rawptr());
  PISM_CHK(ierr, "SNESCreate");

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian.
  m_callback_data.da = *m_da;
  m_callback_data.ssa = this;

#if PETSC_VERSION_LT(3,5,0)
  typedef PetscErrorCode (*JacobianCallback)(DMDALocalInfo*, void*, Mat, Mat, MatStructure*, void*);
  typedef PetscErrorCode (*FunctionCallback)(DMDALocalInfo*, void*, void*, void*);

  ierr = DMDASNESSetFunctionLocal(*m_da, INSERT_VALUES,
                                  (FunctionCallback)function_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  ierr = DMDASNESSetJacobianLocal(*m_da,
                                  (JacobianCallback)jacobian_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetJacobianLocal");
#else
  ierr = DMDASNESSetFunctionLocal(*m_da, INSERT_VALUES,
                                  (DMDASNESFunction)function_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  ierr = DMDASNESSetJacobianLocal(*m_da,
                                  (DMDASNESJacobian)jacobian_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetJacobianLocal");
#endif

  ierr = DMSetMatType(*m_da, "baij");
  PISM_CHK(ierr, "DMSetMatType");

  ierr = DMSetApplicationContext(*m_da, &m_callback_data);
  PISM_CHK(ierr, "DMSetApplicationContext");

  ierr = SNESSetDM(m_snes, *m_da);
  PISM_CHK(ierr, "SNESSetDM");

  // Default of maximum 200 iterations; possibly overridden by command line options
  int snes_max_it = 200;
  ierr = SNESSetTolerances(m_snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                           snes_max_it, PETSC_DEFAULT);
  PISM_CHK(ierr, "SNESSetTolerances");

  ierr = SNESSetFromOptions(m_snes);
  PISM_CHK(ierr, "SNESSetFromOptions");

  // Allocate m_coefficients, which contains coefficient data at the
  // quadrature points of all the elements. There are nElement
  // elements, and Quadrature2x2::Nq quadrature points.
  int nElements = m_element_index.element_count();
  m_coefficients.resize(fem::Quadrature2x2::Nq * nElements);

  m_coeffs.create(m_grid, "ssa_coefficients", WITH_GHOSTS, 1);

  m_node_type.create(m_grid, "node_type", WITH_GHOSTS, 1);
  m_node_type.set_attrs("internal", // intent
                        "node types: interior, boundary, exterior", // long name
                        "", ""); // no units or standard name

  // ElementMap::nodal_values() expects a ghosted IceModelVec2S. Ghosts if this field are never
  // assigned to and not communocated, though.
  m_boundary_integral.create(m_grid, "boundary_integral", WITH_GHOSTS, 1);
  m_boundary_integral.set_attrs("internal", // intent
                                "residual contribution from lateral boundaries", // long name
                                "", ""); // no units or standard name
}

SSA* SSAFEMFactory(IceGrid::ConstPtr g, EnthalpyConverter::Ptr ec) {
  return new SSAFEM(g, ec);
}

SSAFEM::~SSAFEM() {
  // empty
}

// Initialize the solver, called once by the client before use.
void SSAFEM::init_impl() {

  SSA::init_impl();

  // Use explicit driving stress if the surface elevation is not available.
  if (m_surface == NULL) {
    m_driving_stress_x = m_grid->variables().get_2d_scalar("ssa_driving_stress_x");
    m_driving_stress_y = m_grid->variables().get_2d_scalar("ssa_driving_stress_y");
  }

  m_log->message(2,
                 "  [using the SNES-based finite element method implementation]\n");

  // process command-line options
  {
    m_dirichletScale = 1.0e9;
    m_dirichletScale = options::Real("-ssa_fe_dirichlet_scale",
                                     "Enforce Dirichlet conditions with this additional scaling",
                                     m_dirichletScale);

  }

  // On restart, SSA::init() reads the SSA velocity from a PISM output file
  // into IceModelVec2V "velocity". We use that field as an initial guess.
  // If we are not restarting from a PISM file, "velocity" is identically zero,
  // and the call below clears SSAX.

  m_velocity_global.copy_from(m_velocity);

  // Store coefficient data at the quadrature points.
  cache_inputs();
  cache_inputs_new();
}

//! Solve the SSA.  The FEM solver exchanges time for memory by computing
//! the value of the various coefficients at each of the quadrature points
//! only once.  When running in an ice-model context, at each time step,
//! SSA::update is called, which calls SSAFEM::solve.  Since coefficients
//! have generally changed between timesteps, we need to recompute coefficeints
//! at the quad points. On the other hand, in the context of inversion,
//! coefficients will not change between iteration and there is no need to
//! recompute the values at the quad points.  So there are two different solve
//! methods, SSAFEM::solve() and SSAFEM::solve_nocache().  The only difference
//! is that SSAFEM::solve() recomputes the cached values of the coefficients at
//! quadrature points before calling SSAFEM::solve_nocache().
void SSAFEM::solve() {

  TerminationReason::Ptr reason = solve_with_reason();
  if (reason->failed()) {
    throw RuntimeError::formatted("SSAFEM solve failed to converge (SNES reason %s)",
                                  reason->description().c_str());
  } else if (getVerbosityLevel() > 2) {
    m_stdout_ssa += "SSAFEM converged (SNES reason " + reason->description() + ")";
  }
}

TerminationReason::Ptr SSAFEM::solve_with_reason() {

  // Set up the system to solve (store coefficient data at the quadrature points):
  cache_inputs();
  cache_inputs_new();

  return solve_nocache();
}

//! Solve the SSA without first recomputing the values of coefficients at quad
//! points.  See the disccusion of SSAFEM::solve for more discussion.
TerminationReason::Ptr SSAFEM::solve_nocache() {
  PetscErrorCode ierr;

  m_epsilon_ssa = m_config->get_double("epsilon_ssa");

  options::String filename("-ssa_view", "");
  if (filename.is_set()) {
    petsc::Viewer viewer;
    ierr = PetscViewerASCIIOpen(m_grid->com, filename->c_str(), viewer.rawptr());
    PISM_CHK(ierr, "PetscViewerASCIIOpen");

    ierr = PetscViewerASCIIPrintf(viewer, "SNES before SSASolve_FE\n");
    PISM_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = SNESView(m_snes, viewer);
    PISM_CHK(ierr, "SNESView");

    ierr = PetscViewerASCIIPrintf(viewer, "solution vector before SSASolve_FE\n");
    PISM_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = VecView(m_velocity_global.get_vec(), viewer);
    PISM_CHK(ierr, "VecView");
  }

  m_stdout_ssa.clear();
  if (getVerbosityLevel() >= 2) {
    m_stdout_ssa = "  SSA: ";
  }

  // Solve:
  ierr = SNESSolve(m_snes, NULL, m_velocity_global.get_vec());
  PISM_CHK(ierr, "SNESSolve");

  // See if it worked.
  SNESConvergedReason snes_reason;
  ierr = SNESGetConvergedReason(m_snes, &snes_reason); PISM_CHK(ierr, "SNESGetConvergedReason");

  TerminationReason::Ptr reason(new SNESTerminationReason(snes_reason));
  if (not reason->failed()) {

    // Extract the solution back from SSAX to velocity and communicate.
    m_velocity.copy_from(m_velocity_global);
    m_velocity.update_ghosts();

    bool view_solution = options::Bool("-ssa_view_solution", "view solution of the SSA system");
    if (view_solution) {
      petsc::Viewer viewer;
      ierr = PetscViewerASCIIOpen(m_grid->com, filename->c_str(), viewer.rawptr());
      PISM_CHK(ierr, "PetscViewerASCIIOpen");

      ierr = PetscViewerASCIIPrintf(viewer, "solution vector after SSASolve\n");
      PISM_CHK(ierr, "PetscViewerASCIIPrintf");

      ierr = VecView(m_velocity_global.get_vec(), viewer);
      PISM_CHK(ierr, "VecView");
    }

  }

  return reason;
}

//! Initialize stored data from the coefficients in the SSA.  Called by SSAFEM::solve.
/**
   This method is should be called after SSAFEM::init and whenever any geometry or temperature
   related coefficients have changed. The method stores the values of the coefficients at quadrature
   points of each element so that these interpolated values do not need to be computed during each
   outer iteration of the nonlinear solve.

   In addition to coefficients at quadrature points we store "node types" used to identify interior
   elements, exterior elements, and boundary faces.
*/
void SSAFEM::cache_inputs() {

  using fem::Quadrature2x2;
  using fem::Germ;

  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  std::vector<double> Enth_q[Nq];
  const double *Enth_e[4];

  const double
    ice_density      = m_config->get_double("ice_density"),
    standard_gravity = m_config->get_double("standard_gravity"),
    rho_g = ice_density * standard_gravity;

  for (unsigned int q=0; q<Nq; q++) {
    Enth_q[q].resize(m_grid->Mz());
  }

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list;
  list.add(*m_enthalpy);
  bool explicit_driving_stress = (m_driving_stress_x != NULL) && (m_driving_stress_y != NULL);
  if (explicit_driving_stress) {
    list.add(*m_driving_stress_x);
    list.add(*m_driving_stress_y);
  } else {
    // The class SSA ensures in this case that 'surface' is available
    list.add(*m_surface);
  }

  list.add(*m_thickness);
  list.add(*m_bed);
  list.add(*m_tauc);

  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j=ys; j<ys+ym; j++) {
      for (int i=xs; i<xs+xm; i++) {
        double hq[Nq], hxq[Nq], hyq[Nq];
        double ds_xq[Nq], ds_yq[Nq];
        double tmp[Nk];

        m_element.reset(i, j);

        if (explicit_driving_stress) {
          m_element.nodal_values(*m_driving_stress_x, tmp);
          m_quadrature.quadrature_point_values(tmp, ds_xq);

          m_element.nodal_values(*m_driving_stress_y, tmp);
          m_quadrature.quadrature_point_values(tmp, ds_yq);
        } else {
          m_element.nodal_values(*m_surface, tmp);
          m_quadrature.quadrature_point_values(tmp, hq, hxq, hyq);
        }

        double Hq[Nq], bq[Nq], taucq[Nq];
        m_element.nodal_values(*m_thickness, tmp);
        m_quadrature.quadrature_point_values(tmp, Hq);

        m_element.nodal_values(*m_bed, tmp);
        m_quadrature.quadrature_point_values(tmp, bq);

        m_element.nodal_values(*m_tauc, tmp);
        m_quadrature.quadrature_point_values(tmp, taucq);

        const int ij = m_element_index.flatten(i, j);
        Coefficients *coefficients = &m_coefficients[ij*Nq];
        for (unsigned int q = 0; q < Nq; q++) {
          coefficients[q].thickness  = Hq[q];
          coefficients[q].tauc = taucq[q];
          if (explicit_driving_stress) {
            coefficients[q].driving_stress.u = ds_xq[q];
            coefficients[q].driving_stress.v = ds_yq[q];
          } else {
            coefficients[q].driving_stress.u = -rho_g * Hq[q]*hxq[q];
            coefficients[q].driving_stress.v = -rho_g * Hq[q]*hyq[q];
          }

          coefficients[q].mask = gc.mask(m_sea_level, bq[q], coefficients[q].thickness);
        }

        // In the following, we obtain the averaged hardness value from enthalpy by
        // interpolating enthalpy in each column over a quadrature point and then
        // taking the average over the column.  A faster approach would be to take
        // the column average over each element nodes and then interpolate to the
        // quadrature points. Does this make a difference?

        // Obtain the values of enthalpy at each vertical level at each of the vertices
        // of the current element.
        Enth_e[0] = m_enthalpy->get_column(i, j);
        Enth_e[1] = m_enthalpy->get_column(i+1, j);
        Enth_e[2] = m_enthalpy->get_column(i+1, j+1);
        Enth_e[3] = m_enthalpy->get_column(i, j+1);

        // We now want to interpolate to the quadrature points at each of the
        // vertical levels.  It would be nice to use quadrature::computeTestFunctionValues,
        // but the way we have just obtained the values at the element vertices
        // using getInternalColumn doesn't make this straightforward.  So we compute the values
        // by hand.
        const Germ (*test)[Nk] = m_quadrature.test_function_values();
        for (unsigned int k = 0; k < m_grid->Mz(); k++) {
          Enth_q[0][k] = Enth_q[1][k] = Enth_q[2][k] = Enth_q[3][k] = 0;
          for (unsigned int q = 0; q < Nq; q++) {
            for (unsigned int p = 0; p < Nk; p++) {
              Enth_q[q][k] += test[q][p].val * Enth_e[p][k];
            }
          }
        }

        // Now, for each column over a quadrature point, find the averaged_hardness.
        for (unsigned int q = 0; q < Nq; q++) {
          // Evaluate column integrals in flow law at every quadrature point's column
          coefficients[q].hardness = rheology::averaged_hardness(*m_flow_law,
                                                                 coefficients[q].thickness,
                                                                 m_grid->kBelowHeight(coefficients[q].thickness),
                                                                 &(m_grid->z()[0]),
                                                                 &(Enth_q[q])[0]);
        }

      } // j-loop
    } // i-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();

  const bool use_cfbc = m_config->get_boolean("calving_front_stress_boundary_condition");
  if (use_cfbc) {
    // Note: the call below uses ghosts of m_thickness.
    compute_node_types(*m_thickness,
                       m_config->get_double("mask_icefree_thickness_standard"),
                       m_node_type);
  } else {
    m_node_type.set(NODE_INTERIOR);
  }

  cache_residual_cfbc();
}

void SSAFEM::cache_inputs_new() {
  const std::vector<double> &z = m_grid->z();

  IceModelVec::AccessList list;
  list.add(m_coeffs);
  list.add(*m_enthalpy);
  list.add(*m_thickness);
  list.add(*m_bed);
  list.add(*m_tauc);

  bool explicit_driving_stress = (m_driving_stress_x != NULL) && (m_driving_stress_y != NULL);
  if (explicit_driving_stress) {
    list.add(*m_driving_stress_x);
    list.add(*m_driving_stress_y);
  }

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double thickness = (*m_thickness)(i, j);

      Vector2 driving_stress;
      if (explicit_driving_stress) {
        driving_stress.u = (*m_driving_stress_x)(i, j);
        driving_stress.v = (*m_driving_stress_y)(i, j);
      }

      const double *enthalpy = m_enthalpy->get_column(i, j);
      double hardness = rheology::averaged_hardness(*m_flow_law, thickness,
                                                    m_grid->kBelowHeight(thickness),
                                                    &z[0], enthalpy);

      Coeffs c;
      c.thickness      = thickness;
      c.bed            = (*m_bed)(i, j);
      c.sea_level      = m_sea_level; // FIXME: use a 2D field
      c.tauc           = (*m_tauc)(i, j);
      c.hardness       = hardness;
      c.driving_stress = driving_stress;

      m_coeffs(i, j) = c;
    } // loop over owned grid points
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_coeffs.update_ghosts();

  const bool use_cfbc = m_config->get_boolean("calving_front_stress_boundary_condition");
  if (use_cfbc) {
    // Note: the call below uses ghosts of m_thickness.
    compute_node_types(*m_thickness,
                       m_config->get_double("mask_icefree_thickness_standard"),
                       m_node_type);
  } else {
    m_node_type.set(NODE_INTERIOR);
  }

  cache_residual_cfbc();

}

/** @brief Compute the "(regularized effective viscosity) x (ice thickness)" and effective viscous
 *  bed strength from the current solution, at a single quadrature point.
 *
 * @param[in] coefficients SSA coefficients at the current quadrature point
 * @param[in] U the value of the solution
 * @param[in] U_x x-derivatives of velocity components
 * @param[in] U_y y-derivatives of velocity components
 * @param[out] nuH product of the ice viscosity and thickness @f$ \nu H @f$
 * @param[out] dnuH derivative of @f$ \nu H @f$ with respect to the
 *                  second invariant @f$ \gamma @f$. Set to NULL if
 *                  not desired.
 * @param[out] beta basal drag coefficient @f$ \beta @f$
 * @param[out] dbeta derivative of @f$ \beta @f$ with respect to the
 *                   second invariant @f$ \gamma @f$. Set to NULL if
 *                   not desired.
 */
void SSAFEM::PointwiseNuHAndBeta(const Coefficients &c,
                                 const Vector2 &U,
                                 const Vector2 &U_x,
                                 const Vector2 &U_y,
                                 double *nuH, double *dnuH,
                                 double *beta, double *dbeta) {

  if (c.thickness < strength_extension->get_min_thickness()) {
    *nuH = strength_extension->get_notional_strength();
    if (dnuH) {
      *dnuH = 0;
    }
  } else {
    m_flow_law->effective_viscosity(c.hardness, secondInvariant_2D(U_x, U_y),
                                    nuH, dnuH);

    *nuH  = m_epsilon_ssa + *nuH * c.thickness;

    if (dnuH) {
      *dnuH *= c.thickness;
    }
  }

  if (mask::grounded_ice(c.mask)) {
    m_basal_sliding_law->drag_with_derivative(c.tauc, U.u, U.v, beta, dbeta);
  } else {
    *beta = 0;

    if (mask::ice_free_land(c.mask)) {
      *beta = m_beta_ice_free_bedrock;
    }

    if (dbeta) {
      *dbeta = 0;
    }
  }
}


//! Compute and cache residual contributions from the integral over the lateral boundary.
/**

   This method computes FIXME.

 */
void SSAFEM::cache_residual_cfbc() {

  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = fem::BoundaryQuadrature2::Nq;
  const unsigned int n_sides = fem::BoundaryQuadrature2::n_sides;
  using mask::ocean;

  const bool
    use_cfbc          = m_config->get_boolean("calving_front_stress_boundary_condition"),
    is_dry_simulation = m_config->get_boolean("is_dry_simulation");

  const double
    ice_density      = m_config->get_double("ice_density"),
    ocean_density    = m_config->get_double("sea_water_density"),
    standard_gravity = m_config->get_double("standard_gravity");

  // Reset the boundary integral so that all values are overwritten.
  m_boundary_integral.set(0.0);

  if (not use_cfbc) {
    // If CFBC is not used then we're done.
    return;
  }

  fem::BoundaryQuadrature2 bq(m_grid->dx(), m_grid->dy());

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list(m_node_type);
  list.add(*m_thickness);
  list.add(*m_bed);
  list.add(m_boundary_integral);

  // Iterate over the elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {
        // Initialize the map from global to local degrees of freedom for this element.
        m_element.reset(i, j);

        int node_type[Nk];
        m_element.nodal_values(m_node_type, node_type);

        // an element is "interior" if all its nodes are interior or boundary
        const bool interior_element = (node_type[0] < NODE_EXTERIOR and
                                       node_type[1] < NODE_EXTERIOR and
                                       node_type[2] < NODE_EXTERIOR and
                                       node_type[3] < NODE_EXTERIOR);

        // residual contributions at element nodes
        std::vector<Vector2> I(Nk);

        if (not interior_element) {
          // not an interior element: the contribution is zero
          continue;
        }

        double H_nodal[Nk];
        m_element.nodal_values(*m_thickness, H_nodal);

        double b_nodal[Nk];
        m_element.nodal_values(*m_bed, b_nodal);

        // nodes corresponding to a given element "side"
        const int nodes[n_sides][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        // storage for test function values psi[0] for the first "end" of a side, psi[1] for the
        // second
        double psi[2] = {0.0, 0.0};

        // loop over element sides
        for (unsigned int side = 0; side < n_sides; ++side) {

          // nodes incident to the current side
          const int
            n0 = nodes[side][0],
            n1 = nodes[side][1];

          if (not (node_type[n0] == NODE_BOUNDARY and node_type[n1] == NODE_BOUNDARY)) {
            // not a boundary side; skip it
            continue;
          }

          // in our case (i.e. uniform spacing in x and y directions) JxW is the same at all
          // quadrature points along a side.
          const double JxW = bq.weighted_jacobian(side);

          for (unsigned int q = 0; q < Nq; ++q) {

            // test functions at nodes incident to the current side, evaluated at the quadrature point
            // q
            psi[0] = bq.germ(side, q, n0).val;
            psi[1] = bq.germ(side, q, n1).val;

            // Compute ice thickness and bed elevation at a quadrature point. This uses a 1D basis
            // expansion on the side.
            const double
              H   = H_nodal[n0] * psi[0] + H_nodal[n1] * psi[1],
              bed = b_nodal[n0] * psi[0] + b_nodal[n1] * psi[1];

            const bool floating = ocean(gc.mask(m_sea_level, bed, H));

            // ocean pressure difference at a quadrature point
            const double dP = ocean_pressure_difference(floating, is_dry_simulation,
                                                        H, bed, m_sea_level,
                                                        ice_density, ocean_density,
                                                        standard_gravity);

            // This integral contributes to the residual at 2 nodes (the ones incident to the
            // current side). This is is written in a way that allows *adding* (... += ...) the
            // boundary contribution in the residual computation. Note also that we need to worry
            // about the direction of the normal to the lateral boundary.
            switch (side) {
            case 0:             // south
              I[n0].v += - JxW * (- psi[0] * dP);
              I[n1].v += - JxW * (- psi[1] * dP);
              break;
            case 1:             // east
              I[n0].u += JxW * (- psi[0] * dP);
              I[n1].u += JxW * (- psi[1] * dP);
              break;
            case 2:             // north
              I[n0].v += JxW * (- psi[0] * dP);
              I[n1].v += JxW * (- psi[1] * dP);
              break;
            case 3:             // west
              I[n0].u += - JxW * (- psi[0] * dP);
              I[n1].u += - JxW * (- psi[1] * dP);
              break;
            default:
              throw RuntimeError::formatted("invalid side number: %d", side);
            }
          } // q-loop

        } // loop over element sides

        m_element.add_residual_contribution(&I[0], m_boundary_integral);

      } // i-loop
    } // j-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


//! Implements the callback for computing the SNES local function.
/*!
 * Compute the residual \f[r_{ij}= G(x, \psi_{ij}) \f] where \f$G\f$
 * is the weak form of the SSA, \f$x\f$ is the current approximate
 * solution, and the \f$\psi_{ij}\f$ are test functions.
 *
 * The weak form of the SSA system is
 */
void SSAFEM::compute_local_function(Vector2 const *const *const velocity_global,
                                    Vector2 **residual_global) {
  using namespace fem;

  const bool use_cfbc = m_config->get_boolean("calving_front_stress_boundary_condition");

  const unsigned int Nk = ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  IceModelVec::AccessList list(m_node_type);

  // Set the boundary contribution of the residual. This is computed at the nodes, so we don't want
  // to set it using ElementMap::add_residual_contribution() because that would lead to
  // double-counting. Also note that without CFBC m_boundary_integral is exactly zero.
  list.add(m_boundary_integral);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    residual_global[j][i] = m_boundary_integral(i, j);
  }

  // Start access to Dirichlet data if present.
  DirichletData_Vector dirichlet_data(m_bc_mask, m_bc_values, m_dirichletScale);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.weighted_jacobian();

  // Storage for the current solution and its derivatives at quadrature points.
  Vector2 U[Nq], U_x[Nq], U_y[Nq];

  // An Nq by Nk array of test function values.
  const Germ (*test)[Nk] = m_quadrature.test_function_values();

  // Iterate over the elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {

        // Initialize the map from global to local degrees of freedom for this element.
        m_element.reset(i, j);

        int node_type[Nk];
        m_element.nodal_values(m_node_type, node_type);

        // an element is "interior" if all its nodes are interior or boundary
        const bool interior_element = (node_type[0] < NODE_EXTERIOR and
                                       node_type[1] < NODE_EXTERIOR and
                                       node_type[2] < NODE_EXTERIOR and
                                       node_type[3] < NODE_EXTERIOR);

        if (use_cfbc and (not interior_element)) {
          // an exterior element in the CFBC case
          continue;
        }

        // Note: without CFBC all elements are "interior".

        // Storage for the solution and residuals at element nodes.
        Vector2 velocity_nodal[Nk];
        Vector2 residual[Nk];

        // Index into coefficient storage in m_coefficients
        const int ij = m_element_index.flatten(i, j);

        // Coefficients and weights for this quadrature point.
        const Coefficients *coefficients = &m_coefficients[ij*Nq];

        // Obtain the value of the solution at the nodes adjacent to the element.
        m_element.nodal_values(velocity_global, velocity_nodal);

        // These values now need to be adjusted if some nodes in the element have
        // Dirichlet data.
        if (dirichlet_data) {
          // Set elements of velocity_nodal that correspond to Dirichlet nodes to prescribed values.
          dirichlet_data.enforce(m_element, velocity_nodal);
          // mark Dirichlet nodes in m_element so that they are not touched by add_residual_contribution()
          // below
          dirichlet_data.constrain(m_element);
        }

        // Zero out the element-local residual in preparation for updating it.
        for (unsigned int k = 0; k < Nk; k++) {
          residual[k].u = 0;
          residual[k].v = 0;
        }

        // Compute the solution values and symmetric gradient at the quadrature points.
        m_quadrature_vector.quadrature_point_values(velocity_nodal, // input
                                                    U, U_x, U_y);         // outputs

        // loop over quadrature points on this element:
        for (unsigned int q = 0; q < Nq; q++) {

          double eta = 0.0, beta = 0.0;
          PointwiseNuHAndBeta(coefficients[q], U[q], U_x[q], U_y[q], // inputs
                              &eta, NULL, &beta, NULL);              // outputs

          // The next few lines compute the actual residual for the element.
          const Vector2
            tau_b = U[q] * (- beta), // basal shear stress
            tau_d = coefficients[q].driving_stress; // gravitational driving stress

          const double
            u_x          = U_x[q].u,
            v_y          = U_y[q].v,
            u_y_plus_v_x = U_y[q].u + U_x[q].v;

          // Loop over test functions.
          for (unsigned int k = 0; k < Nk; k++) {
            const Germ &psi = test[q][k];

            residual[k].u += JxW[q] * (eta * (psi.dx * (4.0 * u_x + 2.0 * v_y) + psi.dy * u_y_plus_v_x)
                                       - psi.val * (tau_b.u + tau_d.u));
            residual[k].v += JxW[q] * (eta * (psi.dx * u_y_plus_v_x + psi.dy * (2.0 * u_x + 4.0 * v_y))
                                       - psi.val * (tau_b.v + tau_d.v));
          } // k
        } // q

        m_element.add_residual_contribution(residual, residual_global);
      } // j-loop
    } // i-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if (dirichlet_data) {
    dirichlet_data.fix_residual(velocity_global, residual_global);
  }

  if (use_cfbc) {
    // Prescribe homogeneous Dirichlet B.C. at ice-free nodes. This uses the fact that m_node_type
    // can be used as a "Dirichlet B.C. mask", i.e. ice-free nodes (and only ice-free nodes) are
    // marked with ones.
    DirichletData_Vector dirichlet_ice_free(&m_node_type, NULL, m_dirichletScale);
    dirichlet_ice_free.fix_residual_homogeneous(residual_global);
  }

  monitor_function(velocity_global, residual_global);
}

void SSAFEM::monitor_function(Vector2 const *const *const velocity_global,
                              Vector2 const *const *const residual_global) {
  PetscErrorCode ierr;
  bool monitorFunction = options::Bool("-ssa_monitor_function", "monitor the SSA residual");
  if (not monitorFunction) {
    return;
  }

  ierr = PetscPrintf(m_grid->com,
                     "SSA Solution and Function values (pointwise residuals)\n");
  PISM_CHK(ierr, "PetscPrintf");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      ierr = PetscSynchronizedPrintf(m_grid->com,
                                     "[%2d, %2d] u=(%12.10e, %12.10e)  f=(%12.4e, %12.4e)\n",
                                     i, j,
                                     velocity_global[j][i].u, velocity_global[j][i].v,
                                     residual_global[j][i].u, residual_global[j][i].v);
      PISM_CHK(ierr, "PetscSynchronizedPrintf");
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

#if PETSC_VERSION_LT(3,5,0)
  ierr = PetscSynchronizedFlush(m_grid->com);
  PISM_CHK(ierr, "PetscSynchronizedFlush");
#else
  ierr = PetscSynchronizedFlush(m_grid->com, NULL);
  PISM_CHK(ierr, "PetscSynchronizedFlush");
#endif
}


//! Implements the callback for computing the SNES local Jacobian.
/*!
  Compute the Jacobian

  @f[ J_{ij}{kl} \frac{d r_{ij}}{d x_{kl}}= G(x, \psi_{ij}) @f]

  where \f$G\f$ is the weak form of the SSA, \f$x\f$ is the current
  approximate solution, and the \f$\psi_{ij}\f$ are test functions.

*/
void SSAFEM::compute_local_jacobian(Vector2 const *const *const velocity_global, Mat Jac) {

  using fem::Quadrature2x2;
  const unsigned int Nk = fem::ShapeQ1::Nk;
  const unsigned int Nq = Quadrature2x2::Nq;

  const bool use_cfbc = m_config->get_boolean("calving_front_stress_boundary_condition");

  PetscErrorCode ierr;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(Jac);
  PISM_CHK(ierr, "MatZeroEntries");

  IceModelVec::AccessList list(m_node_type);

  // Start access to Dirichlet data if present.
  fem::DirichletData_Vector dirichlet_data(m_bc_mask, m_bc_values, m_dirichletScale);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.weighted_jacobian();

  // Storage for the current solution at quadrature points.
  Vector2 U[Nq], U_x[Nq], U_y[Nq];

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const fem::Germ (*test)[Nk] = m_quadrature.test_function_values();

  // Loop through all the elements.
  int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {
        // Initialize the map from global to local degrees of freedom for this element.
        m_element.reset(i, j);

        int node_type[Nk];
        m_element.nodal_values(m_node_type, node_type);
        // an element is "interior" if all its nodes are interior or boundary
        const bool interior_element = (node_type[0] < NODE_EXTERIOR and
                                       node_type[1] < NODE_EXTERIOR and
                                       node_type[2] < NODE_EXTERIOR and
                                       node_type[3] < NODE_EXTERIOR);

        if (use_cfbc and (not interior_element)) {
          // an exterior element in the CFBC case
          continue;
        }

        // Index into the coefficient storage array.
        const int ij = m_element_index.flatten(i, j);

        // Coefficients at quadrature points in the current element:
        const Coefficients *coefficients = &m_coefficients[ij*Nq];

        // Values of the solution at the nodes of the current element.
        Vector2 velocity_local[Nk];
        // Obtain the value of the solution at the adjacent nodes to the element.
        m_element.nodal_values(velocity_global, velocity_local);

        // These values now need to be adjusted if some nodes in the element have
        // Dirichlet data.
        if (dirichlet_data) {
          dirichlet_data.enforce(m_element, velocity_local);
          dirichlet_data.constrain(m_element);
        }
        // Compute the values of the solution at the quadrature points.
        m_quadrature_vector.quadrature_point_values(velocity_local, U, U_x, U_y);

        // Element-local Jacobian matrix (there are Nk vector valued degrees
        // of freedom per element, for a total of (2*Nk)*(2*Nk) = 16
        // entries in the local Jacobian.
        double K[2*Nk][2*Nk];
        // Build the element-local Jacobian.
        ierr = PetscMemzero(K, sizeof(K));
        PISM_CHK(ierr, "PetscMemzero");

        for (unsigned int q = 0; q < Nq; q++) {
          const double
            jw           = JxW[q],
            u            = U[q].u,
            v            = U[q].v,
            u_x          = U_x[q].u,
            v_y          = U_y[q].v,
            u_y_plus_v_x = U_y[q].u + U_x[q].v;

          double eta = 0.0, deta = 0.0, beta = 0.0, dbeta = 0.0;
          PointwiseNuHAndBeta(coefficients[q], U[q], U_x[q], U_y[q],
                              &eta, &deta, &beta, &dbeta);

          for (unsigned int l = 0; l < Nk; l++) { // Trial functions

            // Current trial function and its derivatives:
            const fem::Germ &phi = test[q][l];

            // Derivatives of \gamma with respect to u_l and v_l:
            const double
              gamma_u = (2.0 * u_x + v_y) * phi.dx + 0.5 * u_y_plus_v_x * phi.dy,
              gamma_v = 0.5 * u_y_plus_v_x * phi.dx + (u_x + 2.0 * v_y) * phi.dy;

            // Derivatives of \eta = \nu*H with respect to u_l and v_l:
            const double
              eta_u = deta * gamma_u,
              eta_v = deta * gamma_v;

            // Derivatives of the basal shear stress term (\tau_b):
            const double
              taub_xu = -dbeta * u * u * phi.val - beta * phi.val,  // x-component, derivative with respect to u_l
              taub_xv = -dbeta * u * v * phi.val,                   // x-component, derivative with respect to u_l
              taub_yu = -dbeta * v * u * phi.val,                   // y-component, derivative with respect to v_l
              taub_yv = -dbeta * v * v * phi.val - beta * phi.val;  // y-component, derivative with respect to v_l

            for (unsigned int k = 0; k < Nk; k++) {   // Test functions

              // Current test function and its derivatives:

              const fem::Germ &psi = test[q][k];

              if (eta == 0) {
                ierr = PetscPrintf(PETSC_COMM_SELF, "eta=0 i %d j %d q %d k %d\n", i, j, q, k);
                PISM_CHK(ierr, "PetscPrintf");
              }

              // u-u coupling
              K[k*2 + 0][l*2 + 0] += jw * (eta_u * (psi.dx * (4 * u_x + 2 * v_y) + psi.dy * u_y_plus_v_x)
                                           + eta * (4 * psi.dx * phi.dx + psi.dy * phi.dy) - psi.val * taub_xu);
              // u-v coupling
              K[k*2 + 0][l*2 + 1] += jw * (eta_v * (psi.dx * (4 * u_x + 2 * v_y) + psi.dy * u_y_plus_v_x)
                                           + eta * (2 * psi.dx * phi.dy + psi.dy * phi.dx) - psi.val * taub_xv);
              // v-u coupling
              K[k*2 + 1][l*2 + 0] += jw * (eta_u * (psi.dx * u_y_plus_v_x + psi.dy * (2 * u_x + 4 * v_y))
                                           + eta * (psi.dx * phi.dy + 2 * psi.dy * phi.dx) - psi.val * taub_yu);
              // v-v coupling
              K[k*2 + 1][l*2 + 1] += jw * (eta_v * (psi.dx * u_y_plus_v_x + psi.dy * (2 * u_x + 4 * v_y))
                                           + eta * (psi.dx * phi.dx + 4 * psi.dy * phi.dy) - psi.val * taub_yv);

            } // l
          } // k
        } // q
        m_element.add_jacobian_contribution(&K[0][0], Jac);
      } // j
    } // i
  } catch (...) {
    loop.failed();
  }
  loop.check();


  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.
  if (dirichlet_data) {
    dirichlet_data.fix_jacobian(Jac);
  }

  if (use_cfbc) {
    // Prescribe homogeneous Dirichlet B.C. at ice-free nodes. This uses the fact that m_node_type
    // can be used as a "Dirichlet B.C. mask", i.e. ice-free nodes (and only ice-free nodes) are
    // marked with ones.
    fem::DirichletData_Vector dirichlet_ice_free(&m_node_type, NULL, m_dirichletScale);
    dirichlet_ice_free.fix_jacobian(Jac);
  }

  ierr = MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyBegin");

  ierr = MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyEnd");

  ierr = MatSetOption(Jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");

  ierr = MatSetOption(Jac, MAT_SYMMETRIC, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");

  monitor_jacobian(Jac);
}

void SSAFEM::monitor_jacobian(Mat Jac) {
  PetscErrorCode ierr;
  bool mon_jac = options::Bool("-ssa_monitor_jacobian", "monitor the SSA Jacobian");

  if (not mon_jac) {
    return;
  }

  // iter has to be a PetscInt because it is used in the
  // SNESGetIterationNumber() call below.
  PetscInt iter = 0;
  ierr = SNESGetIterationNumber(m_snes, &iter);
  PISM_CHK(ierr, "SNESGetIterationNumber");

  char file_name[PETSC_MAX_PATH_LEN];
  snprintf(file_name, PETSC_MAX_PATH_LEN, "PISM_SSAFEM_J%d.m", (int)iter);

  m_log->message(2,
                 "writing Matlab-readable file for SSAFEM system A xsoln = rhs to file `%s' ...\n",
                 file_name);

  petsc::Viewer viewer(m_grid->com);

  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PISM_CHK(ierr, "PetscViewerSetType");

  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PISM_CHK(ierr, "PetscViewerPushFormat");

  ierr = PetscViewerFileSetName(viewer, file_name);
  PISM_CHK(ierr, "PetscViewerFileSetName");

  ierr = PetscObjectSetName((PetscObject) Jac, "A");
  PISM_CHK(ierr, "PetscObjectSetName");

  ierr = MatView(Jac, viewer);
  PISM_CHK(ierr, "MatView");

  ierr = PetscViewerPopFormat(viewer);
  PISM_CHK(ierr, "PetscViewerPopFormat");
}

//!
PetscErrorCode SSAFEM::function_callback(DMDALocalInfo *info,
                                         Vector2 const *const *const velocity,
                                         Vector2 **residual,
                                         CallbackData *fe) {
  try {
    (void) info;
    fe->ssa->compute_local_function(velocity, residual);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)fe->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

#if PETSC_VERSION_LT(3,5,0)
PetscErrorCode SSAFEM::jacobian_callback(DMDALocalInfo *info, const Vector2 **velocity,
                                         Mat A, Mat J, MatStructure *str, CallbackData *fe) {
  try {
    (void) A;
    (void) info;
    fe->ssa->compute_local_jacobian(velocity, J);
    *str = SAME_NONZERO_PATTERN;
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)fe->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}
#else
PetscErrorCode SSAFEM::jacobian_callback(DMDALocalInfo *info,
                                         Vector2 const *const *const velocity,
                                         Mat A, Mat J, CallbackData *fe) {
  try {
    (void) A;
    (void) info;
    fe->ssa->compute_local_jacobian(velocity, J);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)fe->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}
#endif
} // end of namespace stressbalance
} // end of namespace pism
