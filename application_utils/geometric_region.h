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
  * A geometric region within a Nimbus simulation. Defined by integers in
  * x,y,z and dx, dy, dz.
  */

#ifndef NIMBUS_APPLICATION_UTILS_GEOMETRIC_REGION_H_
#define NIMBUS_APPLICATION_UTILS_GEOMETRIC_REGION_H_

#include "shared/nimbus_types.h"

namespace nimbus {

  class GeometricRegion {
  public:
    GeometricRegion(int_dimension_t x,
                    int_dimension_t y,
                    int_dimension_t z,
                    int_dimension_t dx,
                    int_dimension_t dy,
                    int_dimension_t dz);

    virtual ~GeometricRegion() {}


    virtual bool intersects(GeometricRegion* region);

    /* Returns whether this region covers (encompasses) the argument. */
    virtual bool covers(GeometricRegion* region);
    virtual bool adjacent(GeometricRegion* region);
    virtual bool adjacentOrIntersects(GeometricRegion* region);

    int_dimension_t x();
    int_dimension_t y();
    int_dimension_t z();

    int_dimension_t dx();
    int_dimension_t dy();
    int_dimension_t dz();

  private:
    int_dimension_t x_;
    int_dimension_t y_;
    int_dimension_t z_;
    int_dimension_t dx_;
    int_dimension_t dy_;
    int_dimension_t dz_;
  };
}  // namespace nimbus

#endif  // NIMBUS_APPLICATION_UTILS_GEOMETRIC_REGION_H_
