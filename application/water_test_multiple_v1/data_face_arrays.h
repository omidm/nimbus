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

#ifndef NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_DATA_FACE_ARRAYS_H_
#define NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_DATA_FACE_ARRAYS_H_

#include "data_utils.h"
#include "physbam_include.h"
#include "shared/nimbus.h"

#define face_array_id 20

namespace water_app_data {

    typedef float TF;

    /* Face array for storing quantities like face velocities.
    */
    template <class TV>
        class FaceArray : public SimData {

            private:
                typedef typename TV::SCALAR T;
                typedef typename TV::template REBIND<int>::TYPE TV_INT;
                typedef ::PhysBAM::GRID<TV> T_GRID;
                typedef ::PhysBAM::RANGE<TV> T_RANGE;
                typedef
                    ::PhysBAM::ARRAY<T, ::PhysBAM::FACE_INDEX< TV::dimension> >
                    T_FACE_ARRAY;
                typedef ::PhysBAM::BOUNDARY_UNIFORM
                    < ::PhysBAM::GRID<TV>, T> T_BOUNDARY;
                typedef ::PhysBAM::RANGE<TV_INT> T_BOX;

                TV_INT size_;
                T_GRID *grid_;
                T_FACE_ARRAY *data_;

                void Glue_Face_Array(T_FACE_ARRAY *from, T_BOX *box);
                void Update_Face_Array(T_FACE_ARRAY* from, T_BOX *box);

            public:

                FaceArray(TV_INT size);
                virtual void Create();
                virtual void Destroy();
                virtual ::nimbus::Data* Clone();
                virtual int get_debug_info();
                virtual bool Serialize(SerializedData *ser_data);
                virtual bool DeSerialize(const SerializedData& ser_data, Data **result);

                static void Extend_Array(
                        T_FACE_ARRAY *original_vel,
                        T_FACE_ARRAY *extended_vel,
                        T_BOUNDARY *boundary,
                        int bandwidth,
                        TF time,
                        bool extrapolate);

                /* This needs the center region right now due to the way it is
                 * implemented. */
                static void Glue_Regions(
                        T_FACE_ARRAY* result,
                        std::vector<FaceArray* > *parts,
                        int bandwidth,
                        int offset);

                static void Update_Regions(
                        T_FACE_ARRAY *updated,
                        std::vector<FaceArray* > *parts,
                        int bandwidth,
                        int offset);

                T_GRID *grid() {
                    return grid_;
                }

                T_FACE_ARRAY *data() {
                    return data_;
                }

                void set_grid(T_GRID *grid) {
                    grid_ = grid;
                }

                void set_data(T_FACE_ARRAY *data) {
                    data_ = data;
                }
        };

} // namespace water_app_data

#endif // NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_DATA_FACE_ARRAYS_H_
