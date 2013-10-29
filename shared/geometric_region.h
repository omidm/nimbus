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

#ifndef NIMBUS_SHARED_GEOMETRIC_REGION_H_
#define NIMBUS_SHARED_GEOMETRIC_REGION_H_

#include <boost/shared_ptr.hpp>
#include <string>
#include "shared/nimbus_types.h"
#include "shared/protobufs/ldomessage.pb.h"

namespace nimbus {

  class GeometricRegion {
  public:
    GeometricRegion(int_dimension_t x,
                    int_dimension_t y,
                    int_dimension_t z,
                    int_dimension_t dx,
                    int_dimension_t dy,
                    int_dimension_t dz);

    // Argument is a pointer to an array of 6 int_dimension_t
    explicit GeometricRegion(const int_dimension_t* values);
    // Transform a protobuf GeometricRegionMessage into a GeometricRegion
    explicit GeometricRegion(const GeometricRegionMessage* msg);
    // Read in as a serialized GeometricRegionMessage from an istream
    explicit GeometricRegion(std::istream* is);
    // Read in as a serialized GeometricRegionMessage from a string
    explicit GeometricRegion(const std::string& data);

    virtual ~GeometricRegion() {}

    int_dimension_t x();
    int_dimension_t y();
    int_dimension_t z();
    int_dimension_t dx();
    int_dimension_t dy();
    int_dimension_t dz();

    /* Covers returns whether this region covers (encompasses) the argument. */
    virtual bool Covers(GeometricRegion* region);
    virtual bool Intersects(GeometricRegion* region);
    virtual bool Adjacent(GeometricRegion* region);
    virtual bool AdjacentOrIntersects(GeometricRegion* region);

    virtual void FillInMessage(GeometricRegionMessage* msg);

  private:
    int_dimension_t x_;
    int_dimension_t y_;
    int_dimension_t z_;
    int_dimension_t dx_;
    int_dimension_t dy_;
    int_dimension_t dz_;

    void fillInValues(const int_dimension_t* values);
    void fillInValues(const GeometricRegionMessage* msg);
  };
}  // namespace nimbus

#endif  // NIMBUS_SHARED_GEOMETRIC_REGION_H_
