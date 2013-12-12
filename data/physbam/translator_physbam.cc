/*
  * Copyright 2013 Stanford University.
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions
  * are met:
  *
  * - Redistributions of source code must retain the above copyright
  *   notice, this list of conditions and the following disclaimer.
  *
  * - Redistributions in binary form must reproduce the above copyright
  *   notice, this list of conditions and the following disclaimer in the
  *   documentation and/or other materials provided with the
  *   distribution.
  *
  * - Neither the name of the copyright holders nor the names of
  *   its contributors may be used to endorse or promote products derived
  *   from this software without specific prior written permission.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
  * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
  * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
  * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
  * OF THE POSSIBILITY OF SUCH DAMAGE.
  */

/***********************************************************************
 * AUTHOR: Philip Levis <pal>
 *   FILE: .//translator_physbam.cc
 *   DATE: Wed Dec 11 11:29:18 2013
 *  DESCR:
 ***********************************************************************/
#include <algorithm>
#include "data/physbam/translator_physbam.h"

namespace nimbus {
/**
 * \fn nimbus::TranslatorPhysBAM::TranslatorPhysBAM(Worker *worker)
 * \brief Brief description.
 * \param worker
 * \return
*/
template <class VECTOR_TYPE>
TranslatorPhysBAM<VECTOR_TYPE>::TranslatorPhysBAM(Worker *worker)
: worker_(worker) {}


/**
 * \fn FACE_ARRAY_TYPE * nimbus::TranslatorPhysBAM::MakeFaceArray(GeometricRegion *region,
                                         PLdoVector *objects)
 * \brief Brief description.
 * \param region
 * \param objects
 * \return
*/



template <class VECTOR_TYPE>
PhysBAM::ARRAY<typename VECTOR_TYPE::SCALAR, PhysBAM::FACE_INDEX<VECTOR_TYPE::dimension> >*
TranslatorPhysBAM<VECTOR_TYPE>::MakeFaceArray(GeometricRegion *region,
                                              CPdiVector *objects) {
  PhysBAM::VECTOR<int_dimension_t, 3> vector;
  vector(X_COORD) = region->dx();
  vector(Y_COORD) = region->dy();
  vector(Z_COORD) = region->dz();

  // Create a FACE_ARRAY of the right size.
  PhysBAM::RANGE<PhysBAM::VECTOR<int, 3> > range(0, region->dx(), 0, region->dy(), 0, region->dz());

  PhysBAM::ARRAY<typename VECTOR_TYPE::SCALAR, PhysBAM::FACE_INDEX<VECTOR_TYPE::dimension> > fa();
  // FACE_ARRAY_TYPE fa();
  fa.Resize(range);
  CPdiVector::iterator iter = objects->begin();


  /*
  for(; iter != objects->end(); ++iter) {
    PhysicalDataObject* obj = *iter;
    VECTOR_TYPE<int_dimension_t, 3> overlap = GetOverlapSize(obj->region(), region);

  }
  */

  return fa;
}


/**
 * \fn VECTOR<int_dimension_t, 3> nimbus::TranslatorPhysBAM::GetSourceOffset(GeometricRegion *src,
                                           GeometricRegion *dest)
 * \brief Brief description.
 * \param src
 * \param dest
 * \return
*/
template <class VECTOR_TYPE>
PhysBAM::VECTOR<int_dimension_t, 3>
TranslatorPhysBAM<VECTOR_TYPE>::GetSourceOffset(GeometricRegion *src,
                                                GeometricRegion *dest) {
  PhysBAM::VECTOR<int_dimension_t, 3> result;

  // If source is > than dest, its offset is zero (it's contained),
  // otherwise the offset is the difference between the values.
  int_dimension_t x = dest->x() - src->x();
  int_dimension_t y = dest->y() - src->y();
  int_dimension_t z = dest->z() - src->z();
  result(X_COORD) = (x >= 0)? x:0;
  result(Y_COORD) = (y >= 0)? y:0;
  result(Z_COORD) = (z >= 0)? z:0;

  return result;
}


/**
 * \fn Dimension3Vector nimbus::TranslatorPhysBAM::GetDestOffset(GeometricRegion *src,
                                         GeometricRegion *dest)
 * \brief Brief description.
 * \param src
 * \param dest
 * \return
*/
template <class VECTOR_TYPE>
PhysBAM::VECTOR<int_dimension_t, 3>
TranslatorPhysBAM<VECTOR_TYPE>::GetDestOffset(GeometricRegion *src,
                                              GeometricRegion *dest) {
  PhysBAM::VECTOR<int_dimension_t, 3> result;

  // If dest is > than source, its offset is zero (it's contained),
  // otherwise the offset is the difference between the values.
  int_dimension_t x = src->x() - dest->x();
  int_dimension_t y = src->y() - dest->y();
  int_dimension_t z = src->z() - dest->z();
  result(X_COORD) = (x >= 0)? x:0;
  result(Y_COORD) = (y >= 0)? y:0;
  result(Z_COORD) = (z >= 0)? z:0;

  return result;
}


/**
 * \fn Dimension3Vector nimbus::TranslatorPhysBAM::GetOverlapSize(GeometricRegion *src,
                                          GeometricRegion *dest)
 * \brief Brief description.
 * \param src
 * \param dest
 * \return
*/
template <class VECTOR_TYPE>
PhysBAM::VECTOR<int_dimension_t, 3>
TranslatorPhysBAM<VECTOR_TYPE>::GetOverlapSize(GeometricRegion *src,
                                               GeometricRegion *dest) {
  PhysBAM::VECTOR<int_dimension_t, 3> result;

  int_dimension_t x_start = std::max(src->x(), dest->x());
  int_dimension_t x_end   = std::min(src->x() + src->dx(), dest->x() + dest->dx());
  int_dimension_t x_size = x_end - x_start;

  int_dimension_t y_start = std::max(src->y(), dest->y());
  int_dimension_t y_end   = std::min(src->y() + src->dy(), dest->y() + dest->dy());
  int_dimension_t y_size = y_end - y_start;

  int_dimension_t z_start = std::max(src->y(), dest->z());
  int_dimension_t z_end   = std::min(src->z() + src->dz(), dest->z() + dest->dz());
  int_dimension_t z_size = z_end - z_start;

  result(X_COORD) = (x_size >= 0)? x_size:0;
  result(Y_COORD) = (y_size >= 0)? y_size:0;
  result(Z_COORD) = (z_size >= 0)? z_size:0;

  return result;
}

/**
 * \fn bool nimbus::TranslatorPhysBAM::HasOverlap(VECTOR<int_dimension_t, 3>)
 * \brief Brief description.
 * \param src
 * \param dest
 * \return
*/
template <class VECTOR_TYPE>
bool TranslatorPhysBAM<VECTOR_TYPE>::HasOverlap(Dimension3Vector vector) {
  return (vector(X_COORD) > 0 && vector(Y_COORD) > 0 && vector(Z_COORD) > 0);
}

}  // namespace nimbus
