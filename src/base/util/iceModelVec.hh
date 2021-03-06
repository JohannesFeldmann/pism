// Copyright (C) 2008--2015 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#ifndef __IceModelVec_hh
#define __IceModelVec_hh

#include <petscvec.h>
#include <gsl/gsl_interp.h>

#include "pism_memory.hh"
#include "VariableMetadata.hh"
#include "base/util/petscwrappers/Viewer.hh"
#include "Vector2.hh"
#include "StarStencil.hh"
#include "base/util/petscwrappers/DM.hh"
#include "base/util/petscwrappers/Vec.hh"
#include "base/util/IceGrid.hh"

namespace pism {

class IceGrid;
class PIO;

//! What "kind" of a vector to create: with or without ghosts.
enum IceModelVecKind {WITHOUT_GHOSTS=0, WITH_GHOSTS=1};

struct Range {
  double min, max;
};

//! \brief Abstract class for reading, writing, allocating, and accessing a
//! DA-based PETSc Vec (2D and 3D fields) from within IceModel.
/*!
  @anchor icemodelvec_use

  This class represents 2D and 3D fields in PISM. Its methods common to all
  the derived classes can be split (roughly) into six kinds:

  - memory allocation (create)
  - point-wise access (begin_access(), end_access())
  - arithmetic (range(), norm(), add(), shift(), scale(), set(), ...)
  - setting or reading metadata (set_attrs(), metadata())
  - file input/output (read, write, regrid)
  - tracking whether a field was updated (get_state_counter(), inc_state_counter())

  ## Memory allocation

  Creating an IceModelVec... object does not allocate memory for storing it
  (some IceModelVecs serve as "references" and don't have their own storage).
  To complete IceModelVec... creation, use the "create()" method:

  \code
  IceModelVec2S var;
  ierr = var.create(grid, "var_name", WITH_GHOSTS); CHKERRQ(ierr);
  // var is ready to use
  \endcode

  ("WITH_GHOSTS" means "can be used in computations using map-plane neighbors
  of grid points.)

  It is usually a good idea to set variable metadata right after creating it.
  The method set_attrs() is used throughout PISM to set commonly used
  attributes.

  ## Point-wise access

  PETSc performs some pointer arithmetic magic to allow convenient indexing of
  grid point values. Because of this one needs to surround the code using row,
  column or level indexes with begin_access() and end_access() calls:

  \code
  double foo;
  int i = 0, j = 0;
  IceModelVec2S var;
  // assume that var was allocated
  ierr = var.begin_access(); CHKERRQ(ierr);
  foo = var(i,j) * 2;
  ierr = var.end_access(); CHKERRQ(ierr);
  \endcode

  Please see [this page](@ref computational_grid) for a discussion of the
  organization of PISM's computational grid and examples of for-loops you will
  probably put between begin_access() and end_access().

  To ensure that ghost values are up to date add the following call
  before the code using ghosts:

  \code
  ierr = var.update_ghosts(); CHKERRQ(ierr);
  \endcode

  ## Reading and writing variables

  PISM can read variables either from files with data on a grid matching the
  current grid (read()) or, using bilinear interpolation, from files
  containing data on a different (but compatible) grid (regrid()).

  To write a field to a "prepared" NetCDF file, use write(). (A file is prepared
  if it contains all the necessary dimensions, coordinate variables and global
  metadata.)

  If you need to "prepare" a file, do:
  \code
  PIO nc(grid.com, grid.config.get_string("output_format"));

  std::string time_name = config.get_string("time_dimension_name");
  ierr = nc.open(filename, PISM_READWRITE); CHKERRQ(ierr); // append == false
  ierr = nc.def_time(time_name, grid.time->calendar(),
  grid.time->CF_units_string()); CHKERRQ(ierr);
  ierr = nc.append_time(time_name, grid.time->current()); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  \endcode

  A note about NetCDF write performance: due to limitations of the NetCDF
  (classic, version 3) format, it is significantly faster to
  \code
  for (all variables)
  var.define(...);

  for (all variables)
  var.write(...);
  \endcode

  as opposed to

  \code
  for (all variables) {
  var.define(...);
  var.write(...);
  }
  \endcode

  IceModelVec::define() is here so that we can use the first approach.

  ## Tracking if a field changed

  It is possible to track if a certain field changed with the help of
  get_state_counter() and inc_state_counter() methods.

  For example, PISM's SIA code re-computes the smoothed bed only if the bed
  deformation code updated it:

  \code
  if (bed->get_state_counter() > bed_state_counter) {
  ierr = bed_smoother->preprocess_bed(...); CHKERRQ(ierr);
  bed_state_counter = bed->get_state_counter();
  }
  \endcode

  The state counter is **not** updated automatically. For the code snippet above
  to work, a bed deformation model has to call inc_state_counter() after an
  update.
*/
class IceModelVec {
public:
  IceModelVec();
  virtual ~IceModelVec();

  typedef PISM_SHARED_PTR(IceModelVec) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec) ConstPtr;


  virtual bool was_created() const;
  IceGrid::ConstPtr get_grid() const;
  unsigned int get_ndims() const;
  //! \brief Returns the number of degrees of freedom per grid point.
  unsigned int get_ndof() const;
  unsigned int get_stencil_width() const;
  std::vector<double> get_levels() const;

  virtual Range range() const;
  double norm(int n) const;
  std::vector<double> norm_all(int n) const;
  virtual void  add(double alpha, const IceModelVec &x);
  virtual void  squareroot();
  virtual void  shift(double alpha);
  virtual void  scale(double alpha);
  // This is used in Python code (as a local-to-global replacement),
  // but we should be able to get rid of it.
  void copy_to_vec(petsc::DM::Ptr destination_da, Vec destination) const;
  void copy_from_vec(Vec source);
  virtual void copy_from(const IceModelVec &source);
  Vec get_vec();
  petsc::DM::Ptr get_dm() const;
  virtual void  set_name(const std::string &name);
  virtual const std::string& get_name() const;
  virtual void  set_attrs(const std::string &my_pism_intent, const std::string &my_long_name,
                          const std::string &my_units, const std::string &my_standard_name,
                          int component = 0);
  virtual void  read_attributes(const std::string &filename, int component = 0);
  virtual void  define(const PIO &nc, IO_Type output_datatype = PISM_DOUBLE) const;

  void read(const std::string &filename, unsigned int time);
  void read(const PIO &nc, unsigned int time);

  void  write(const std::string &filename) const;
  void  write(const PIO &nc) const;

  void  regrid(const std::string &filename, RegriddingFlag flag,
               double default_value = 0.0);
  void  regrid(const PIO &nc, RegriddingFlag flag,
               double default_value = 0.0);

  virtual void  begin_access() const;
  virtual void  end_access() const;
  virtual void  update_ghosts();
  virtual void  update_ghosts(IceModelVec &destination) const;

  void  set(double c);

  SpatialVariableMetadata& metadata(unsigned int N = 0);

  const SpatialVariableMetadata& metadata(unsigned int N = 0) const;

  int get_state_counter() const;
  void inc_state_counter();
  void set_time_independent(bool flag);

  bool   m_report_range;                 //!< If true, report range when regridding.
  bool   write_in_glaciological_units;
  //!< \brief If true, data is written to a file in "human-friendly" units.

protected:

  void global_to_local(petsc::DM::Ptr dm, Vec source, Vec destination) const;
  virtual void read_impl(const PIO &nc, unsigned int time);
  virtual void regrid_impl(const PIO &nc, RegriddingFlag flag,
                                     double default_value = 0.0);
  virtual void write_impl(const PIO &nc) const;
  std::vector<double> zlevels;

  petsc::Vec  m_v;                       //!< Internal storage
  std::string m_name;

  //! stores metadata (NetCDF variable attributes)
  std::vector<SpatialVariableMetadata> m_metadata;

  IceGrid::ConstPtr m_grid;

  unsigned int m_dof;                     //!< number of "degrees of freedom" per grid point
  unsigned int m_da_stencil_width;      //!< stencil width supported by the DA
  bool m_has_ghosts;            //!< m_has_ghosts == true means "has ghosts"
  petsc::DM::Ptr m_da;          //!< distributed mesh manager (DM)

  bool begin_end_access_use_dof;

  //! It is a map, because a temporary IceModelVec can be used to view
  //! different quantities
  mutable std::map<std::string,petsc::Viewer::Ptr> map_viewers;

  mutable void *array;  // will be cast to double** or double*** in derived classes

  mutable int m_access_counter;           // used in begin_access() and end_access()
  int m_state_counter;            //!< Internal IceModelVec "revision number"

  virtual void checkCompatibility(const char *function, const IceModelVec &other) const;

  //! \brief Check the array indices and warn if they are out of range.
  void check_array_indices(int i, int j, unsigned int k) const;
  void reset_attrs(unsigned int N);
  NormType int_to_normtype(int input) const;

  void get_dof(petsc::DM::Ptr da_result, Vec result, unsigned int n,
               unsigned int count=1) const;
  void set_dof(petsc::DM::Ptr da_source, Vec source, unsigned int n,
               unsigned int count=1);
private:
  // disable copy constructor and the assignment operator:
  IceModelVec(const IceModelVec &other);
  IceModelVec& operator=(const IceModelVec&);
public:
  //! Dump an IceModelVec to a file. *This is for debugging only.*
  //! Uses const char[] to make it easier to call it from gdb.
  void dump(const char filename[]) const;

public:

  //! Makes sure that we call begin_access() and end_access() for all accessed IceModelVecs.
  class AccessList {
  public:
    AccessList();
    AccessList(const IceModelVec &v);
    ~AccessList();
    void add(const IceModelVec &v);
  private:
    std::vector<const IceModelVec*> m_vecs;
  };
};

class IceModelVec2S;

/** Class for a 2d DA-based Vec.

    As for the difference between IceModelVec2 and IceModelVec2S, the
    former can store fields with more than 1 "degree of freedom" per grid
    point (such as 2D fields on the "staggered" grid, with the first
    degree of freedom corresponding to the i-offset and second to
    j-offset). */
class IceModelVec2 : public IceModelVec {
public:
  IceModelVec2();

  typedef PISM_SHARED_PTR(IceModelVec2) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec2) ConstPtr;

  virtual void view(int viewer_size) const;
  virtual void view(petsc::Viewer::Ptr v1, petsc::Viewer::Ptr v2) const;
  // component-wise access:
  virtual void get_component(unsigned int n, IceModelVec2S &result) const;
  virtual void set_component(unsigned int n, const IceModelVec2S &source);
  inline double& operator() (int i, int j, int k);
  inline const double& operator() (int i, int j, int k) const;
  void create(IceGrid::ConstPtr my_grid, const std::string &my_short_name,
              IceModelVecKind ghostedp, unsigned int stencil_width, int dof);
protected:
  virtual void read_impl(const PIO &nc, const unsigned int time);
  virtual void regrid_impl(const PIO &nc, RegriddingFlag flag,
                                     double default_value = 0.0);
  virtual void write_impl(const PIO &nc) const;
};

class IceModelVec2V;

/** A class for storing and accessing scalar 2D fields.
    IceModelVec2S is just IceModelVec2 with "dof == 1" */
class IceModelVec2S : public IceModelVec2 {
  friend class IceModelVec2V;
  friend class IceModelVec2Stag;
public:
  IceModelVec2S();

  typedef PISM_SHARED_PTR(IceModelVec2S) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec2S) ConstPtr;

  static Ptr To2DScalar(IceModelVec::Ptr input);

  // does not need a copy constructor, because it does not add any new data members
  using IceModelVec2::create;
  void create(IceGrid::ConstPtr my_grid, const std::string &my_name,
              IceModelVecKind ghostedp, int width = 1);
  petsc::Vec::Ptr allocate_proc0_copy() const;
  void put_on_proc0(Vec onp0) const;
  void get_from_proc0(Vec onp0);
  virtual void copy_from(const IceModelVec &source);
  double** get_array();
  virtual void set_to_magnitude(const IceModelVec2S &v_x, const IceModelVec2S &v_y);
  virtual void set_to_magnitude(const IceModelVec2V &input);
  virtual void mask_by(const IceModelVec2S &M, double fill = 0.0);
  virtual void add(double alpha, const IceModelVec &x);
  virtual void add(double alpha, const IceModelVec &x, IceModelVec &result) const;
  virtual double sum() const;
  virtual double min() const;
  virtual double max() const;
  virtual double absmax() const;
  virtual double diff_x(int i, int j) const;
  virtual double diff_y(int i, int j) const;
  virtual double diff_x_stagE(int i, int j) const;
  virtual double diff_y_stagE(int i, int j) const;
  virtual double diff_x_stagN(int i, int j) const;
  virtual double diff_y_stagN(int i, int j) const;
  virtual double diff_x_p(int i, int j) const;
  virtual double diff_y_p(int i, int j) const;

  //! Provides access (both read and write) to the internal double array.
  /*!
    Note that i corresponds to the x direction and j to the y.
  */
  inline double& operator() (int i, int j);
  inline const double& operator()(int i, int j) const;
  inline StarStencil<double> star(int i, int j) const;
protected:
  void put_on_proc0(Vec parallel, Vec onp0) const;
  void get_from_proc0(Vec onp0, Vec parallel);
};


//! \brief A simple class "hiding" the fact that the mask is stored as
//! floating-point scalars (instead of integers).
class IceModelVec2Int : public IceModelVec2S {
public:

  typedef PISM_SHARED_PTR(IceModelVec2Int) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec2Int) ConstPtr;

  inline int as_int(int i, int j) const;
  inline StarStencil<int> int_star(int i, int j) const;
};

/** Class for storing and accessing 2D vector fields used in IceModel.
    IceModelVec2V is IceModelVec2 with "dof == 2". (Plus some extra methods, of course.)
*/
class IceModelVec2V : public IceModelVec2 {
public:
  IceModelVec2V();
  ~IceModelVec2V();

  typedef PISM_SHARED_PTR(IceModelVec2V) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec2V) ConstPtr;

  static Ptr ToVector(IceModelVec::Ptr input);

  void create(IceGrid::ConstPtr my_grid, const std::string &my_short_name,
              IceModelVecKind ghostedp, unsigned int stencil_width = 1);
  virtual void copy_from(const IceModelVec &source);
  virtual void add(double alpha, const IceModelVec &x);
  virtual void add(double alpha, const IceModelVec &x, IceModelVec &result) const;

  // I/O:
  Vector2** get_array();
  inline Vector2& operator()(int i, int j);
  inline const Vector2& operator()(int i, int j) const;
  inline StarStencil<Vector2> star(int i, int j) const;
};

//! \brief A class for storing and accessing internal staggered-grid 2D fields.
//! Uses dof=2 storage. This class is identical to IceModelVec2V, except that
//! components are not called `u` and `v` (to avoid confusion).
class IceModelVec2Stag : public IceModelVec2 {
public:
  IceModelVec2Stag();

  typedef PISM_SHARED_PTR(IceModelVec2Stag) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec2Stag) ConstPtr;

  static Ptr ToStaggered(IceModelVec::Ptr input);

  void create(IceGrid::ConstPtr my_grid, const std::string &my_short_name,
              IceModelVecKind ghostedp, unsigned int stencil_width = 1);
  virtual void staggered_to_regular(IceModelVec2S &result) const;
  virtual void staggered_to_regular(IceModelVec2V &result) const;
  virtual std::vector<double> absmaxcomponents() const;

  //! Returns the values at interfaces of the cell i,j using the staggered grid.
  /*! The ij member of the return value is set to 0, since it has no meaning in
    this context.
  */
  inline StarStencil<double> star(int i, int j) const;
};

//! \brief A virtual class collecting methods common to ice and bedrock 3D
//! fields.
class IceModelVec3D : public IceModelVec {
public:
  IceModelVec3D();
  virtual ~IceModelVec3D();

  void set_column(int i, int j, double c);
  void set_column(int i, int j, double *valsIN);
  double* get_column(int i, int j);
  const double* get_column(int i, int j) const;

  virtual double getValZ(int i, int j, double z) const;
  virtual bool isLegalLevel(double z) const;

  inline double& operator() (int i, int j, int k);
  inline const double& operator() (int i, int j, int k) const;
protected:
  void allocate(IceGrid::ConstPtr mygrid, const std::string &my_short_name,
                IceModelVecKind ghostedp, const std::vector<double> &levels,
                unsigned int stencil_width = 1);
private:
  gsl_interp_accel *m_bsearch_accel;
};


//! Class for a 3d DA-based Vec for ice scalar quantities.
class IceModelVec3 : public IceModelVec3D {
public:
  IceModelVec3();
  virtual ~IceModelVec3();

  typedef PISM_SHARED_PTR(IceModelVec3) Ptr;
  typedef PISM_SHARED_PTR(const IceModelVec3) ConstPtr;

  static Ptr To3DScalar(IceModelVec::Ptr input);

  void create(IceGrid::ConstPtr mygrid, const std::string &my_short_name,
              IceModelVecKind ghostedp,
              unsigned int stencil_width = 1);

  void  getHorSlice(Vec &gslice, double z) const; // used in iMmatlab.cc
  void  getHorSlice(IceModelVec2S &gslice, double z) const;
  void  getSurfaceValues(IceModelVec2S &gsurf, const IceModelVec2S &myH) const;
};

/** 
 * Convert a PETSc Vec from the units in `from` into units in `to` (in place).
 *
 * @param v data to convert
 * @param system unit system
 * @param spec1 source unit specification string
 * @param spec2 destination unit specification string 
 */
void convert_vec(Vec v, units::System::Ptr system,
                 const std::string &spec1, const std::string &spec2);

} // end of namespace pism

// include inline methods; contents are wrapped in namespace pism {...}
#include "IceModelVec_inline.hh"

#endif /* __IceModelVec_hh */

