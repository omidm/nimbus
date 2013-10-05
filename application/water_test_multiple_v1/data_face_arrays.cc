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

#include "data_face_arrays.h"
#include "physbam_include.h"
#include "proto_files/app_face_array_2d.pb.h"
#include "proto_files/physbam_serialize_data_arrays_2d.h"
#include "proto_files/physbam_serialize_data_common_2d.h"
#include "proto_files/physbam_deserialize_data_arrays_2d.h"
#include "proto_files/physbam_deserialize_data_common_2d.h"
#include "shared/nimbus.h"
#include "string.h"


namespace water_app_data {

    template <class TV> FaceArray<TV>::
        FaceArray(int size) : size_(size), grid(0), data(0) {
        }

    template <class TV> void FaceArray<TV>::
        Create() {
            std::cout << "Creating FaceArray\n";
            if (!grid)
                delete grid;
            if (!data)
                delete data;
            if (size_ == 0) {
                grid = new T_GRID();
                data = new T_FACE_ARRAY();
            }
            else {
                grid = new T_GRID(TV_INT::All_Ones_Vector()*size_,
                        T_RANGE::Unit_Box(), true);
                assert(grid);
                data = new  T_FACE_ARRAY(*grid);
                assert(data);
            }
        }

    template <class TV> void FaceArray<TV>::
        Destroy() {
            printf("Destroying face array\n");
            delete grid;
            delete data;
        }

    template <class TV> ::nimbus::Data* FaceArray<TV>::
        Clone() {
            std::cout << "Cloning facearray\n";
            return new FaceArray<TV>(size_);
        }

    template <class TV> int FaceArray<TV>::
        get_debug_info() {
            return face_array_id;
        }

    template <class TV> bool FaceArray<TV>::
        Serialize(SerializedData *ser_data) {
            assert(ser_data);
            ::communication::AppFaceArray2d pb_fa;
            pb_fa.set_size(size_);
            if (grid != NULL) {
                ::communication::PhysbamGrid2d *pb_fa_grid =
                    pb_fa.mutable_grid();
                ::physbam_pb::make_pb_object(grid, pb_fa_grid);
            }
            if (data != NULL) {
                ::communication::PhysbamFaceArray2d *pb_fa_data =
                    pb_fa.mutable_data();
                ::physbam_pb::make_pb_object(data, pb_fa_data);
            }
            std::string ser;
            if (!pb_fa.SerializeToString(&ser)) {
                ser_data->set_size(0);
                return false;
            }
            ser_data->set_size(ser.length());
            char *buffer = new char(ser.length() + 1);
            strcpy(buffer, ser.c_str());
            ser_data->set_data_ptr(buffer);
            return true;
        }

    // DeSerialize should be called only after Create has been called. Create
    // ensures that required memory has been allocated.
    template <class TV> bool FaceArray<TV>::
        DeSerialize(const SerializedData &ser_data, Data **result) {
            assert(result);
            const char *buffer = ser_data.data_ptr();
            const int buff_size = ser_data.size();
            if (buff_size <= 0)
                return false;
            assert(buffer);
            ::communication::AppFaceArray2d pb_fa;
            *result = new FaceArray<TV>(0);
            FaceArray<TV> *des = (FaceArray<TV> *)(*result);
            if (*result == NULL)
                return false;
            pb_fa.ParseFromString((std::string)buffer);
            if (pb_fa.has_grid()) {
                ::physbam_pb::make_physbam_object(des->grid, pb_fa.grid());
            }
            if (pb_fa.has_data())
                ::physbam_pb::make_physbam_object(des->data, pb_fa.data());
            return true;
        }

    template <class TV> void FaceArray<TV>::
        Initialize_Ghost_Regions(
                T_FACE_ARRAY *extended_vel,
                int bandwidth,
                bool extrpolate) {
            (*extended_vel) = (*data);
            extended_vel->Resize(*grid,  bandwidth, false);
        }

    template <class TV> void FaceArray<TV>::
        Put_Face_Array(FaceArray<TV>* to,
                FaceArray<TV>* from, ::PhysBAM::RANGE<TV_INT2>& box) {
            TV_INT2 i, j;
            for (int axis = 1; axis <= 2; axis++) {
                int dx = 0;
                int dy = 0;
                if (axis == 1)
                    dx = 1;
                else
                    dy = 1;

                j.x = 1;
                for(i.x = box.min_corner.x; i.x <= (box.max_corner.x + dx); i.x++)
                {
                    j.y = 1;
                    for(i.y = box.min_corner.y; i.y <= (box.max_corner.y + dy); i.y++) {
                        (*(to->data))(axis, i) = (*(from->data))(axis, j);
                        j.y++;
                    }
                    j.x++;
                }
            }
        }

    template <class TV> void FaceArray<TV>::
        Fill_Ghost_Cells(FaceArray<TV>* result,
                std::vector<FaceArray<TV>* > parts, int bandwidth) {

            TV_INT2 min_corner, max_corner;
            ::PhysBAM::RANGE<TV_INT2> box;

            int len_x = result->grid->numbers_of_cells(1);
            int len_y = result->grid->numbers_of_cells(2);

            min_corner = TV_INT2(1 - bandwidth, 1 - bandwidth);
            max_corner = TV_INT2(0, 0);
            box = ::PhysBAM::RANGE<TV_INT2>(min_corner, max_corner);
            Put_Face_Array(result, parts[1], box);

            min_corner = TV_INT2(len_x + 1, 1 - bandwidth);
            max_corner = TV_INT2(len_x + bandwidth, 0);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[3], box);

            min_corner = TV_INT2(1 - bandwidth, len_y + 1);
            max_corner = TV_INT2(0, len_y + bandwidth);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner,max_corner);
            Put_Face_Array(result, parts[6], box);

            min_corner = TV_INT2(len_x + 1, len_y + 1);
            max_corner = TV_INT2(len_x + bandwidth, len_y + bandwidth);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[8], box);

            min_corner = TV_INT2(1, 1 - bandwidth);
            max_corner = TV_INT2(len_x, 0);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[2], box);

            min_corner = TV_INT2(1 - bandwidth, 1);
            max_corner = TV_INT2(0, len_y);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[4], box);

            min_corner = TV_INT2(len_x + 1, 1);
            max_corner = TV_INT2(len_x + bandwidth, len_y);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[5], box);

            min_corner = TV_INT2(1, len_y + 1);
            max_corner = TV_INT2(len_x, len_y + bandwidth);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[7], box);

            min_corner = TV_INT2(1, 1);
            max_corner = TV_INT2(len_x, len_y);
            box = ::PhysBAM::RANGE<TV_INT2 >(min_corner, max_corner);
            Put_Face_Array(result, parts[0], box);
        }

    namespace {
        typedef ::PhysBAM::VECTOR<float, 2> TVF2;
    } // namespace
    template class ::water_app_data::FaceArray<TVF2>;

} // namespace water_app_data
