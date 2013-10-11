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

    namespace {
        typedef ::PhysBAM::VECTOR<float, 2> TVF2;
    } // namespace

    template <class TV> FaceArray<TV>::
        FaceArray(TV_INT size) :
            size_(size),
            grid_(0),
            data_(0) {
        }

    template <class TV> void FaceArray<TV>::
        Create() {
            std::cout << "Creating FaceArray\n";
            if (!grid_)
                delete grid_;
            if (!data_)
                delete data_;
            grid_ = new T_GRID(size_,
                    T_RANGE::Unit_Box(),
                    true);
            assert(grid_);
            data_ = new  T_FACE_ARRAY(*grid_);
            assert(data_);
        }

    template <class TV> void FaceArray<TV>::
        Destroy() {
            printf("Destroying face array\n");
            delete grid_;
            delete data_;
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
            ::communication::PhysbamVectorInt2d *pb_fa_size = 
                pb_fa.mutable_size();
            ::physbam_pb::make_pb_object(&size_, pb_fa_size);
            if (grid_ != NULL) {
                ::communication::PhysbamGrid2d *pb_fa_grid =
                    pb_fa.mutable_grid();
                ::physbam_pb::make_pb_object(grid_, pb_fa_grid);
            }
            if (data_ != NULL) {
                ::communication::PhysbamFaceArray2d *pb_fa_data =
                    pb_fa.mutable_data();
                ::physbam_pb::make_pb_object(data_, pb_fa_data);
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

    /* DeSerialize should be called only after Create has been called. Create
       ensures that required memory has been allocated. */
    template <class TV> bool FaceArray<TV>::
        DeSerialize(const SerializedData &ser_data, Data **result) {
            assert(result);
            const char *buffer = ser_data.data_ptr();
            const int buff_size = ser_data.size();
            if (buff_size <= 0)
                return false;
            assert(buffer);
            ::communication::AppFaceArray2d pb_fa;
            pb_fa.ParseFromString((std::string)buffer);
            TV_INT size;
            ::physbam_pb::make_physbam_object(&size, pb_fa.size());
            *result = new FaceArray<TV>(size);
            if (*result == NULL)
                return false;
            FaceArray<TV> *des = (FaceArray<TV> *)(*result);
            des->Create();
            if (pb_fa.has_grid()) {
                ::physbam_pb::make_physbam_object(des->grid_, pb_fa.grid());
            }
            if (pb_fa.has_data())
                ::physbam_pb::make_physbam_object(des->data_, pb_fa.data());
            return true;
        }

    template <class TV> void FaceArray<TV>::
        Extend_Array(
                T_FACE_ARRAY *original_vel,
                T_FACE_ARRAY *extended_vel,
                T_BOUNDARY *boundary,
                int bandwidth,
                TF time,
                bool extrapolate) {
            TV_INT new_min = original_vel->domain_indices.min_corner;
            TV_INT new_max = original_vel->domain_indices.max_corner;
            T_BOX extended_range = T_BOX(new_min-bandwidth, new_max + bandwidth);
            extended_vel->Resize(extended_range, true, false, 0);
            if (extrapolate) {
                T_GRID extended_grid(
                        new_max - new_min + 1,
                        T_RANGE::Unit_Box(),
                        true);
                boundary->Fill_Ghost_Cells_Face(
                        extended_grid,
                        *original_vel,
                        *extended_vel,
                        time,
                        bandwidth); // this also copies the center values
            }
            else
                T_FACE_ARRAY::Put(*original_vel, *extended_vel);
        }

    template <class TV> void FaceArray<TV>::
        Glue_Face_Array(
                T_FACE_ARRAY* from,
                T_BOX& box) {
            TV_INT i, j;
            for (int axis = 1; axis <= 2; axis++) {
                int dx = axis == 1? 1 : 0;
                int dy = axis == 2? 1 : 0;
                for(j.x = 1, i.x = box.min_corner.x;
                        i.x <= (box.max_corner.x + dx); i.x++, j.x++) {
                    for(j.y = 1, i.y = box.min_corner.y;
                            i.y <= (box.max_corner.y + dy); i.y++, j.y++) {
                        (*(data_))(axis, i) = (*(from))(axis, j);
                    }
                }
            }
        }

    /* This needs the center region right now due to the way it is
     * implemented. */
    template <class TV> void FaceArray<TV>::
        Glue_Regions(
                T_FACE_ARRAY* result,
                std::vector<FaceArray * > parts,
                int bandwidth,
                int offset = 0) {
            TV_INT min_corner, max_corner;
            T_BOX box;
            TV_INT max_d = result->domain_indices.Maximum_Corner() - offset;
            TV_INT min_d = result->domain_indices.Minimum_Corner() + offset;
            int len_x = max_d.x - min_d.x + 1;
            int len_y = max_d.y - min_d.y + 1;
            if (parts[1]) {
                box = T_BOX(1-bandwidth, 0, 1-bandwidth, 0);
                parts[1]->Glue_Face_Array(result,  box);
            }
            if (parts[3]) {
                box = T_BOX(len_x+1, len_x+bandwidth, 1-bandwidth, 0);
                parts[3]->Glue_Face_Array(result, box);
            }
            if (parts[6]) {
                box = T_BOX(1-bandwidth, 0, len_y+1, len_y+bandwidth);
                parts[6]->Glue_Face_Array(result, box);
            }
            if (parts[8]) {
                box = T_BOX
                    (len_x+1, len_x+bandwidth, len_y+1, len_y+bandwidth);
                parts[8]->Glue_Face_Array(result, box);
            }
            if (parts[2]) {
                box = T_BOX(1, len_x, 1-bandwidth, 0);
                parts[2]->Glue_Face_Array(result, box);
            }
            if (parts[4]) {
                box = T_BOX(1-bandwidth, 0, 1, len_y);
                parts[4]->Glue_Face_Array(result, box);
            }
            if (parts[5]) {
                box = T_BOX(len_x+1, len_x+bandwidth, 1, len_y);
                parts[5]->Glue_Face_Array(result, box);
            }
            if (parts[7]) {
                box = T_BOX
                    (1, len_x, len_y+1, len_y+bandwidth);
                parts[7]->Glue_Face_Array(result, box);
            }
            assert(parts[0]);
            box = T_BOX(1, len_x, 1, len_y);
            parts[0]->Glue_Face_Array(result, box);
        }

    template <class TV> void FaceArray<TV>::
        Update_Face_Array(T_FACE_ARRAY* from, T_BOX& box) {
            TV_INT i, j;
            for (int axis = 1; axis <= 2; axis++) {
                int dx = axis == 1? 1 : 0;
                int dy = axis == 2? 1 : 0;
                for(j.x = 1, i.x = box.min_corner.x;
                        i.x <= (box.max_corner.x + dx); i.x++, j.x++) {
                    for(j.y = 1, i.y = box.min_corner.y;
                            i.y <= (box.max_corner.y + dy); i.y++, j.y++) {
                        (*(data_))(axis, j) = (*(from))(axis, j);
                    }
                }
            }
        }

    template <class TV> void FaceArray<TV>::
        Update_Regions(
                T_FACE_ARRAY* updated,
                std::vector<FaceArray* > parts,
                int bandwidth) {
            TV_INT min_corner, max_corner;
            T_BOX box;
            TV_INT max_d = updated->domain_indices.Maximum_Corner();
            TV_INT min_d = updated->domain_indices.Minimum_Corner();
            int len_x = max_d.x - min_d.x + 1;
            int len_y = max_d.y - min_d.y + 1;
            if (parts[1]) {
                box = T_BOX(1-bandwidth, 0, 1-bandwidth, 0);
                parts[1]->Update_Face_Array(updated, box);
            }
            if (parts[3]) {
                box = T_BOX(len_x+1, len_x+bandwidth, 1-bandwidth, 0);
                parts[3]->Update_Face_Array(updated, box);
            }
            if (parts[6]) {
                box = T_BOX(1-bandwidth, 0, len_y+1, len_y+bandwidth);
                parts[6]->Update_Face_Array(updated, box);
            }
            if (parts[8]) {
                box = T_BOX
                    (len_x+1, len_x+bandwidth, len_y+1, len_y+bandwidth);
                parts[8]->Update_Face_Array(updated, box);
            }
            if (parts[2]) {
                box = T_BOX(1, len_x, 1-bandwidth, 0);
                parts[2]->Update_Face_Array(updated, box);
            }
            if (parts[4]) {
                box = T_BOX(1-bandwidth, 0, 1, len_y);
                parts[4]->Update_Face_Array(updated, box);
            }
            if (parts[5]) {
                box = T_BOX(len_x+1, len_x+bandwidth, 1, len_y);
                parts[5]->Update_Face_Array(updated, box);
            }
            if (parts[7]) {
                box = T_BOX
                    (1, len_x, len_y+1, len_y+bandwidth);
                parts[7]->Update_Face_Array(updated, box);
            }
            assert(parts[0]);
            box = T_BOX(1, len_x, 1, len_y);
            parts[0]->Update_Face_Array(updated, box);
        }

    template class ::water_app_data::FaceArray<TVF2>;

} // namespace water_app_data
