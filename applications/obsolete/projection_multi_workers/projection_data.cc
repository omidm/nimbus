/*
 * The internal data structures of PhysBAM projection.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include "projection_data.h"

namespace PhysBAM {
template class ProjectionData<VECTOR<float,2> >;
template class ProjectionData<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ProjectionData<VECTOR<double,2> >;
template class ProjectionData<VECTOR<double,3> >;
#endif
template class ProjectionInternalData<VECTOR<float,2> >;
template class ProjectionInternalData<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ProjectionInternalData<VECTOR<double,2> >;
template class ProjectionInternalData<VECTOR<double,3> >;
#endif
}  // namespace PhysBAM
