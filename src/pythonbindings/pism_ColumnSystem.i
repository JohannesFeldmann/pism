/* Classes related to column-wise systems and interpolation in the
 * column (from the storage grid to the fine computational grid).
 *
 * We use these Python wrappers to test C++ code.
 */

%{
#include "base/energy/enthSystem.hh"
#include "base/util/ColumnInterpolation.hh"
%}

/* wrap the enthalpy solver to make testing easier */
%include "base/columnSystem.hh"
%rename(get_lambda) pism::energy::enthSystemCtx::lambda;
%include "base/energy/enthSystem.hh"

%include "base/util/ColumnInterpolation.hh"
