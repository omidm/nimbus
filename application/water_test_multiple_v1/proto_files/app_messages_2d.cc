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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "app_face_array_2d.pb.h"
#include "app_messages_2d.h"
#include "physbam_data_include.h"
#include "physbam_serialize_data_arrays_2d.h"
#include "physbam_deserialize_data_arrays_2d.h"
#include "shared/geometric_region.h"

namespace water_app_data {

    // serialize
    void make_pb_object(
            ::physbam_pb::FaceArray2 *fa,
            ::nimbus::GeometricRegion *region,
            ::communication::AppFaceArray2d *app_fa) {
        assert(fa);
        assert(region);
        assert(app_fa);
        ::communication::PhysbamFaceArray2d *pb_fa = app_fa->mutable_data();
        ::communication::GeometricRegionMessage *pb_region
            = app_fa->mutable_region();
        ::physbam_pb::make_pb_object(fa, pb_fa);
        make_pb_object(region, pb_region);
    }

    void make_pb_object(
            ::nimbus::GeometricRegion *region,
            ::communication::GeometricRegionMessage *rm) {
        rm->set_x(region->x());
        rm->set_y(region->y());
        rm->set_z(region->z());
        rm->set_dx(region->dx());
        rm->set_dy(region->dy());
        rm->set_dz(region->dz());
    }

    // deserialize
    void make_app_object(
            ::physbam_pb::FaceArray2 *fa,
            ::nimbus::GeometricRegion *region,
            const ::communication::AppFaceArray2d &app_fa) {
        ::communication::PhysbamFaceArray2d pb_fa = app_fa.data();
        ::communication::GeometricRegionMessage r
            = app_fa.region();
        ::physbam_pb::make_physbam_object(fa, pb_fa);
        make_app_object(region, &r);
    }

    void make_app_object(
            ::nimbus::GeometricRegion *region,
            ::communication::GeometricRegionMessage *rm) {
        region->Rebuild(rm->x(), rm->y(), rm->z(), rm->dx(), rm->dy(), rm->dz());
    }

} // namespace physbam_pb
