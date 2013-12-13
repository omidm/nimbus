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
#include "data/physbam/physbam_data.h"

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
    typedef typename PhysBAM::ARRAY<scalar_t, FaceIndex > FaceArray;

    enum {
      X_COORD = 1,
      Y_COORD = 2,
      Z_COORD = 3
    };

    explicit TranslatorPhysBAM() {}
    virtual ~TranslatorPhysBAM() {}

    virtual FaceArray* ReadFaceArray(GeometricRegion* region,
                                     CPdiVector* objects) {
      Dimension3Vector vec;
      vec(X_COORD) = region->dx();
      vec(Y_COORD) = region->dy();
      vec(Z_COORD) = region->dz();

      // Create a FACE_ARRAY of the right size.
      PhysBAM::RANGE<PhysBAM::VECTOR<int, 3> > range(0, region->dx(),
                                                     0, region->dy(),
                                                     0, region->dz());

      FaceArray* fa = new FaceArray();
      fa->Resize(range);

      if (objects != NULL) {
        CPdiVector::iterator iter = objects->begin();
        for (; iter != objects->end(); ++iter) {
          const PhysicalDataInstance* obj = *iter;
          Dimension3Vector overlap = GetOverlapSize(obj->region(), region);
          if (HasOverlap(overlap)) {
            dbg(DBG_TRANSLATE, "Incorporating physical object %lu into FaceArray.\n", obj->id());
            PhysBAMData* data = static_cast<PhysBAMData*>obj->data();
            scalar_t* buffer = static_cast<scalar_t*>(data->Buffer());

            Dimension3Vector dest = GetOffset(region, obj->region());
            Dimension3Vector src = GetOffset(obj->region(), region);

            for (int x = 0; x < overlap(X_COORD); x++) {
              for (int y = 0; y < overlap(Y_COORD); y++) {
                for (int z = 0; z < overlap(Z_COORD); z++) {
                  int source_x = x + src(X_COORD);
                  int source_y = y + src(Y_COORD);
                  int source_z = z + src(Z_COORD);

                  // The underlying Nimbus data objects are stored as
                  // arrays of structs: the x, y, and z face values for
                  // a given cell are stored contiguously. Using arrays
                  // of structs in case PhysBAM uses struct of arrays;
                  // that way only 4 cache lines will be used.
                  int source_index =
                    (source_z * (obj->region()->y() * obj->region()->x())) +
                    (source_y * (obj->region()->x())) +
                    source_x;
                  source_index *= 3;  // We are in three dimensions

                  int dest_x = x + dest(X_COORD);
                  int dest_y = y + dest(Y_COORD);
                  int dest_z = z + dest(Z_COORD);
                  dest_index =
                    (dest_z * (region->y() * region->x())) +
                    (dest_y * (region->x())) +
                     dest_x;

                  FaceArray::TV_INT destinationIndex(dest_x, dest_y, dest_z);

                  // The PhysBAM FACE_ARRAY object abstracts away whether
                  // the data is stored in struct of array or array of struct
                  // form (in practice, usually struct of arrays.
                  fa(X_COORD, destinationIndex) = buffer[source_index];
                  fa(Y_COORD, destinationIndex) = buffer[source_index + 1];
                  fa(Z_COORD, destinationIndex) = buffer[source_index + 2];
                }
              }
            }
          }
        }
      }

      return fa;
    }

    /** Take a FaceArray described by region and write it out to the
     *  PhysicalDataInstance objects in the objects array. */
    virtual bool WriteFaceArray(GeometricRegion* region,
                                CPdiVector* objects,
                                FaceArray* fa) {
      int_dimension_t region_size = region->dx() * region->dy() *
                                    region->dz() * 3;
      if (region_size == fa.buffer_size) {
        dbg(DBG_ERROR, "ERROR: writing a face array of size %i for a region of size %i and the two sizes should be equal.\n", fa.buffer_size, region_size);  // NOLINT
        return FALSE;
      }

      if (objects != NULL) {
        CPdiVector::iterator iter = objects->begin();

        // Loop over the Nimbus objects, copying the relevant PhysBAM data
        // into each one
        for (; iter != objects->end(); ++iter) {
          const PhysicalDataInstance* obj = *iter;
          Dimension3Vector overlap = GetOverlapSize(obj->region(), region);
          if (!HasOverlap(overlap)) {continue;}

          dbg(DBG_TRANSLATE, "Saving FaceArray into physical object %lu.\n", obj->id());
          PhysBAMData* data = static_cast<PhysBAMData*>obj->data();
          scalar_t* buffer = data->buffer();
          Dimension3Vector src = GetOffset(region, obj->region());
          Dimension3Vector dest = GetOffset(obj->region(), region);

          for (int x = 0; x < overlap(X_COORD); x++) {
            for (int y = 0; y < overlap(Y_COORD); y++) {
              for (int z = 0; z < overlap(Z_COORD); z++) {
                int dest_x = x + dest(X_COORD);
                int dest_y = y + dest(Y_COORD);
                int dest_z = z + dest(Z_COORD);

                // The underlying Nimbus data objects are stored as
                // arrays of structs: the x, y, and z face values for
                // a given cell are stored contiguously. Using arrays
                // of structs in case PhysBAM uses struct of arrays;
                // that way only 4 cache lines will be used.
                int destination_index =
                  (dest_z * (obj->region()->y() * obj->region()->x())) +
                  (dest_y * (obj->region()->x())) +
                  dest_x;
                destination_index *= 3;  // We are in three dimensions

                int source_x = x + src(X_COORD);
                int source_y = y + src(Y_COORD);
                int source_z = z + src(Z_COORD);
                source_index =
                  (source_z * (source->y() * source->x())) +
                  (source_y * (source->x())) +
                  source_x;

                FaceArray::TV_INT sourceIndex(source_x, source_y, source_z);

                // The PhysBAM FACE_ARRAY object abstracts away whether
                // the data is stored in struct of array or array of struct
                // form (in practice, usually struct of arrays
                buffer[destination_index]     = fa(X_COORD, sourceIndex);
                buffer[destination_index + 1] = fa(Y_COORD, sourceIndex);
                buffer[destination_index + 2] = fa(Z_COORD, sourceIndex);
              }
            }
          }
        }
      }
      delete fa;
    }

  private:
    /* Return a vector describing what the offset of b
       within a, such that a.x + offset = b.x. If
       offset is negative, return 0. */
    virtual Dimension3Vector GetOffset(GeometricRegion* a,
                                       GeometricRegion* b) {
      Dimension3Vector result;

      // If source is > than dest, its offset is zero (it's contained),
      // otherwise the offset is the difference between the values.
      int_dimension_t x = b->x() - a->x();
      int_dimension_t y = b->y() - a->y();
      int_dimension_t z = b->z() - a->z();
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
      return (overlapSize(X_COORD) > 0 &&
              overlapSize(Y_COORD) > 0 &&
              overlapSize(Z_COORD) > 0);
    }
};

}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
