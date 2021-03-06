project (Pism)
cmake_minimum_required (VERSION 2.8.0) # 2.8.0 is almost 3 years old, so this should be OK
if (NOT ${CMAKE_VERSION} VERSION_LESS 3.0.0)
  message(STATUS "CMAKE_VERSION >= 3.0.0. Setting policy 0026 behavior to 'old'...")
  # This has to do with the way we compute dependencies of PISM's
  # Python bindings. Enabling the "old" behavior allows us to support
  # both old and new (>= 3.0) CMake versions.
  cmake_policy(SET CMP0026 OLD)
endif()

include ("CMake/PISM_CMake_macros.cmake")
list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")

set (Pism_BRANCH "stable")

# Set Pism_REVISION_TAG
pism_set_revision_tag()

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Set the install prefix
pism_set_install_prefix()

# Set Pism_DEFAULT_CONFIG_FILE (*after* we set the CMAKE_INSTALL_PREFIX above).
pism_check_build_dir_location()
set (Pism_DEFAULT_CONFIG_FILE "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}/pism_config.nc" CACHE STRING "" FORCE)
mark_as_advanced (Pism_DEFAULT_CONFIG_FILE)
file (WRITE ${PROJECT_BINARY_DIR}/.petscrc "-config ${PROJECT_BINARY_DIR}/pism_config.nc")

# The default options cache
option (Pism_BUILD_EXTRA_EXECS "Build extra executables (mostly testing/verification)" OFF)
option (BUILD_SHARED_LIBS      "Build shared Pism libraries"                           ON)
option (Pism_BUILD_PYTHON_BINDINGS "Build python bindings" OFF)
option (Pism_USE_PROJ4 "Use Proj.4 to compute cell areas, longitudes, and latitudes." OFF)
option (Pism_USE_PARALLEL_NETCDF4 "Enables parallel NetCDF-4 I/O." OFF)
option (Pism_USE_PNETCDF "Enables parallel NetCDF-3 I/O using PnetCDF." OFF)
option (Pism_USE_PARALLEL_HDF5 "Enables parallel HDF5 I/O." OFF)
option (Pism_USE_TAO "Use TAO in inverse solvers." OFF)

option (Pism_TEST_USING_VALGRIND "Add extra regression tests using valgrind" OFF)
mark_as_advanced (Pism_TEST_USING_VALGRIND)

option (Pism_ADD_FPIC "Add -fPIC to C++ compiler flags (CMAKE_CXX_FLAGS). Try turning it off if it does not work." ON)
option (Pism_LINK_STATICALLY "Set CMake flags to try to ensure that everything is linked statically")
option (Pism_LOOK_FOR_LIBRARIES "Specifies whether PISM should look for libraries. (Disable this on Crays.)" ON)
option (Pism_USE_TR1 "Use the std::tr1 namespace to access shared pointer definitions. Disable to get shared pointers from the std namespace (might be needed with some compilers)." ON)

# Use rpath by default; this has to go first, because rpath settings may be overridden later.
pism_use_rpath()

if (Pism_LINK_STATICALLY)
  pism_strictly_static()
endif ()

# Deal with build types
set (Pism_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
# Force the default cache entry to obey Pism_BUILD_TYPE
set (CMAKE_BUILD_TYPE ${Pism_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
mark_as_advanced (CMAKE_BUILD_TYPE) # hide the default entry
if (Pism_BUILD_TYPE MATCHES "Debug")
  set (Pism_BUILD_EXTRA_EXECS ON CACHE BOOL "Build extra executables (mostly testing/verification)" FORCE)
  option (Pism_DEBUG              "Enables extra checks in the code."          ON)
  option (Pism_PEDANTIC_WARNINGS  "Compile with pedantic warnings."            ON)
  option (Pism_GPROF_FLAGS        "Add flags necessary to profile with gprof." OFF)
endif (Pism_BUILD_TYPE MATCHES "Debug")

# Enable extra checks in the code:
if (Pism_DEBUG)
  message (STATUS "Adding -DPISM_DEBUG=1 to compiler flags.")
  add_definitions (-DPISM_DEBUG=1)
else()
  message (STATUS "Adding -DPISM_DEBUG=0 to compiler flags.")
  add_definitions (-DPISM_DEBUG=0)
endif ()

# Add -fPIC to CXX flags; note that this will not show up in CMakeCache.txt
if (Pism_ADD_FPIC)
  if (NOT CMAKE_CXX_FLAGS MATCHES "-fPIC")
    message (STATUS "Adding -fPIC to C++ compiler flags.")
    set (CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compiler flags" FORCE)
  endif()
endif ()

if (Pism_PEDANTIC_WARNINGS)
  pism_set_pedantic_flags()
endif (Pism_PEDANTIC_WARNINGS)

if (Pism_GPROF_FLAGS)
  set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
endif ()

# Look for libraries using find_package(...), etc. Run CMake with -DPism_LOOK_FOR_LIBRARIES=OFF
# to build on systems that rely on the module system to set all compiler and linker flags.
if (Pism_LOOK_FOR_LIBRARIES)
  pism_find_prerequisites()
endif()

# Set Pism_EXTERNAL_LIBS and include directories.
pism_set_dependencies()

if (Pism_BUILD_PYTHON_BINDINGS)
  find_package(Python REQUIRED)
  find_package(PETSc4Py REQUIRED)
  find_package(SWIG REQUIRED)

  mark_as_advanced (SWIG_DIR SWIG_EXECUTABLE SWIG_VERSION)

  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import distutils.sysconfig; print distutils.sysconfig.get_python_lib(plat_specific=True, prefix='${CMAKE_INSTALL_PREFIX}')"
    OUTPUT_VARIABLE PISM_INSTALL_PYTHON_MODULE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(PISM_INSTALL_PYTHON_MODULE_DIR ${PISM_INSTALL_PYTHON_MODULE_DIR}
    CACHE PATH "Python extension module installation directory." )
endif ()

# Do cell area computations the right way if we have proj.4.
if (Pism_USE_PROJ4)
  add_definitions (-DPISM_USE_PROJ4=1)
else()
  add_definitions (-DPISM_USE_PROJ4=0)
endif()

# Use TAO in inverse solvers.
if (Pism_USE_TAO)
  add_definitions (-DPISM_USE_TAO=1)
endif()

if (Pism_USE_TR1)
  message (STATUS "Adding -DPISM_USE_TR1=1 to compiler flags.")
  add_definitions (-DPISM_USE_TR1=1)
endif ()

add_custom_target (etags
  COMMAND find . -regextype posix-extended -regex ".*\\.(c|cc|h|hh)" | xargs etags
  WORKING_DIRECTORY ${Pism_SOURCE_DIR}
  VERBATIM
)

# re-run tests that failed
add_custom_target (retest
  COMMAND ${CMAKE_CTEST_COMMAND} --rerun-failed --output-on-failure
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  VERBATIM
)

# run Python tests
add_custom_target (test-python
  COMMAND ${CMAKE_CTEST_COMMAND} -R "Python:"
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  VERBATIM
)

# Install helper scripts residing in util/
install (PROGRAMS
  util/check_stationarity.py
  util/fill_missing.py
  util/flowline.py
  util/nc2cdo.py
  util/nc2mat.py
  util/nccmp.py
  util/pism_config_editor.py
  DESTINATION ${Pism_BIN_DIR})

add_subdirectory (src)
add_subdirectory (site-packages)

if (Pism_BUILD_PDFS)
  add_subdirectory (doc)
else()
  add_subdirectory (doc EXCLUDE_FROM_ALL)
endif()

# PISM regression testing
ENABLE_TESTING()
include(CTest)
add_subdirectory (test)
add_subdirectory (test/regression)
