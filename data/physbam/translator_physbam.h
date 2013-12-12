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
  * This class requires a pointer to the Application class because it
  * needs to be able to translate LogicalDataObjects into physical
  * data objects.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
#define NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_

#include "data/physbam/physbam_include.h"

#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "worker/physical_data_instance.h"
#include "worker/worker.h"

namespace nimbus {

  template <class VECTOR_TYPE> class TranslatorPhysBAM {
  public:
    typedef VECTOR_TYPE TV;
    typedef typename TV::SCALAR SCALAR_TYPE;
    typedef typename PhysBAM::FACE_INDEX<TV::dimension> FaceIndex;
    typedef typename PhysBAM::ARRAY<SCALAR_TYPE, FaceIndex > FACE_ARRAY_TYPE; // NOLINT
    typedef PhysBAM::VECTOR<int_dimension_t, 3> Dimension3Vector;

    enum {
      X_COORD = 0,
      Y_COORD = 1,
      Z_COORD = 2
    };

    explicit TranslatorPhysBAM(Worker* worker);
    virtual ~TranslatorPhysBAM() {}

    /* Produce an array of scalars fitting the geometric region,
       based on the data in the vector of objects. Returns NULL on
       an error.*/
    virtual FACE_ARRAY_TYPE* MakeFaceArray(GeometricRegion* region,
                                           CPdiVector* objects);

  private:
    Worker* worker_;

    virtual Dimension3Vector GetSourceOffset(GeometricRegion* src,
                                            GeometricRegion* dest);

    virtual Dimension3Vector GetDestOffset(GeometricRegion* src,
                                           GeometricRegion* dest);

    virtual Dimension3Vector GetOverlapSize(GeometricRegion* src,
                                            GeometricRegion* dest);

    virtual bool HasOverlap(Dimension3Vector overlapSize);
  };

}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
