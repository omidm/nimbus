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

 /*
  * TranslatorPhysBAM is a class for translating Nimbus physical objects into
  * PhysBAM data objects. It is templated to make it a little easier
  * to handle PhysBAM's generality, taking a simple template parameter
  * of a PhysBAM VECTOR. This VECTOR is typically a 2D or 3D float: it
  * represents a point in space. The class derives the scalar type
  * (typically float) from this VECTOR, as well as the dimensionality.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
#define NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_

#include <algorithm>

#include "data/physbam/physbam_include.h"

#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "worker/physical_data_instance.h"
#include "worker/worker.h"

namespace nimbus {

template <class VECTOR_TYPE> class TranslatorPhysBAM {
  public:
    typedef VECTOR_TYPE TV;
    typedef typename TV::SCALAR scalar_t;
    typedef PhysBAM::VECTOR<int_dimension_t, 3> Dimension3Vector;
    typedef typename PhysBAM::FACE_INDEX<TV::dimension> FaceIndex;
    typedef typename PhysBAM::ARRAY<SCALAR_TYPE, FaceIndex > FACE_ARRAY_TYPE;

    enum {
      X_COORD = 0,
      Y_COORD = 1,
      Z_COORD = 2
    };

    explicit TranslatorPhysBAM() {}
    virtual ~TranslatorPhysBAM() {}

    /* Produce an array of scalars fitting the geometric region,
       based on the data in the vector of objects. Returns NULL on
       an error.*/
    virtual FACE_ARRAY_TYPE* MakeFaceArray(GeometricRegion* region,
                                           CPdiVector* objects);

 private:
    /* Returns the X,Y,Z offset in the source of the overlap of
       the src and dest regions.*/
    virtual Dimension3Vector GetSourceOffset(GeometricRegion* src,
                                             GeometricRegion* dest);
    /* Returns the X,Y,Z offset in the dest of the overlap of
       the src and dest regions.*/
    virtual Dimension3Vector GetDestOffset(GeometricRegion* src,
                                           GeometricRegion* dest);
    virtual bool HasOverlap(Dimension3Vector vector);

 public:
    virtual FACE_ARRAY_TYPE* MakeFaceArray(GeometricRegion* region,
                                           CPdiVector* objects) {
      Dimension3Vector vector;
      vector(X_COORD) = region->dx();
      vector(Y_COORD) = region->dy();
      vector(Z_COORD) = region->dz();

      // Create a FACE_ARRAY of the right size.
      PhysBAM::RANGE<PhysBAM::VECTOR<int, 3> > range(0, region->dx(),
                                                     0, region->dy(),
                                                     0, region->dz());

      FACE_ARRAY_TYPE fa();
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

  private:
    virtual Dimension3Vector GetSourceOffset(GeometricRegion* src,
                                             GeometricRegion* dest);
    virtual Dimension3Vector GetDestOffset(GeometricRegion* src,
                                             GeometricRegion* dest);
    virtual bool HasOverlap(Dimension3Vector vector);


    virtual Dimension3Vector GetSourceOffset(GeometricRegion* src,
                                             GeometricRegion* dest) {
      Dimension3Vector result;

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

    virtual Dimension3Vector GetDestOffset(GeometricRegion* src,
                                           GeometricRegion* dest) {
      Dimension3Vector result;

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

    virtual Dimension3Vector GetOverlapSize(GeometricRegion* src,
                                            GeometricRegion* dest) {
      Dimension3Vector result;

      int_dimension_t x_start = std::max(src->x(), dest->x());
      int_dimension_t x_end   = std::min(src->x() + src->dx(),
                                         dest->x() + dest->dx());
      int_dimension_t x_size = x_end - x_start;

      int_dimension_t y_start = std::max(src->y(), dest->y());
      int_dimension_t y_end   = std::min(src->y() + src->dy(),
                                         dest->y() + dest->dy());
      int_dimension_t y_size = y_end - y_start;

      int_dimension_t z_start = std::max(src->y(), dest->z());
      int_dimension_t z_end   = std::min(src->z() + src->dz(),
                                         dest->z() + dest->dz());
      int_dimension_t z_size = z_end - z_start;

      result(X_COORD) = (x_size >= 0)? x_size:0;
      result(Y_COORD) = (y_size >= 0)? y_size:0;
      result(Z_COORD) = (z_size >= 0)? z_size:0;

      return result;
    }

    virtual bool HasOverlap(Dimension3Vector overlapSize) {
      return (vector(X_COORD) > 0 &&
              vector(Y_COORD) > 0 &&
              vector(Z_COORD) > 0);
    }
};

}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
