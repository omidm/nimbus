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
#include "proto_files/app_messages_2d.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"
#include "string.h"

namespace water_app_data {

    namespace {
        typedef ::PhysBAM::VECTOR<float, 2> TVF2;
    } // namespace

    template <class TV> FaceArray<TV>::
        FaceArray(::nimbus::GeometricRegion &region, bool lor) :
            SimData(region),
            grid_(0),
            data_(0),
            size_(region.dx(), region.dy()),
            left_or_right(lor){};

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
            std::cout << "Cloning facearray with region "<<region()->ToNetworkData()<<std::endl;
            return (new FaceArray<TV>(*region(), left_or_right));
        }

    template <class TV> int FaceArray<TV>::
        get_debug_info() {
            return face_array_id;
        }

    /* We need to serialize only region and data. Information from these is
     * sufficient to reconstruct the object.
     */
    template <class TV> bool FaceArray<TV>::
        Serialize(SerializedData *ser_data) {
            assert(ser_data);
            ::communication::AppFaceArray2d pb_fa;
            make_pb_object(data(), region(), &pb_fa);
            std::string ser;
            bool success = pb_fa.SerializeToString(&ser);
            char *buf = new char[ser.length()];
            memcpy(buf, ser.c_str(), sizeof(char) * (ser.length()+1));
            if (success)
                ser_data->set_data_ptr(buf);
            ser_data->set_size(sizeof(char) * (ser.length()+1));
            return success;
        }

    /* DeSerialize should be called only after Create has been called. Create
       ensures that required memory has been allocated, and grid and size have
       been initialized correctly. We need to update only region and data. */
    template <class TV> bool FaceArray<TV>::
        DeSerialize(const SerializedData &ser_data, Data **result) {
            assert(result);
            const char *buffer = ser_data.data_ptr_raw();
            const int buff_size = ser_data.size();
            if (buff_size <= 0)
                return false;
            assert(buffer);
            ::communication::AppFaceArray2d pb_fa;
            std::string temp(buffer, buff_size/sizeof(char)-1);
            pb_fa.ParseFromString(temp);
            FaceArray<TV> *des = (FaceArray<TV> *)(*result);
            make_app_object(des->data(), des->region(), pb_fa);
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
                T_FACE_ARRAY* result,
                T_BOX *box) {
            TV_INT i, j;
            for (int axis = 1; axis <= 2; axis++) {
                int dx = axis == 1? 1 : 0;
                int dy = axis == 2? 1 : 0;
                for(j.x = 1, i.x = box->min_corner.x;
                        i.x <= (box->max_corner.x + dx); i.x++, j.x++) {
                    for(j.y = 1, i.y = box->min_corner.y;
                            i.y <= (box->max_corner.y + dy); i.y++, j.y++) {
                        (*(result))(axis, i) = (*(data_))(axis, j);
                        //result->Component(axis)(i) = data()->Component(axis)(j);
                    }
                }
            }
        }

    template <class TV> void FaceArray<TV>::
        Glue_Regions(
                T_FACE_ARRAY* result,
                std::vector<FaceArray * > &parts,
                ::nimbus::GeometricRegion &region,
                int o_x = 0,
                int o_y = 0) {
            T_BOX box;
            unsigned int parts_num = parts.size();
            for (unsigned int i = 0; i < parts_num; i++)
            {
                FaceArray *part = parts[i];
                ::nimbus::GeometricRegion *r = part->region();
                // allow copy from only covered regions
                if (region.Covers(r)) {
                    box = T_BOX(
                            r->x() + o_x,
                            r->x()+r->dx() + o_x -1,
                            r->y() + o_y,
                            r->y()+r->dy() + o_y -1);
                    part->Glue_Face_Array(result, &box);
                }
            }
        }

    template <class TV> void FaceArray<TV>::
        Update_Face_Array(T_FACE_ARRAY* from, T_BOX *box) {
            TV_INT i, j;
            for (int axis = 1; axis <= 2; axis++) {
                int dx = axis == 1? 1 : 0;
                int dy = axis == 2? 1 : 0;
                for(j.x = 1, i.x = box->min_corner.x;
                        i.x <= (box->max_corner.x + dx); i.x++, j.x++) {
                    for(j.y = 1, i.y = box->min_corner.y;
                            i.y <= (box->max_corner.y + dy); i.y++, j.y++) {
                        (*(data_))(axis, j) = (*(from))(axis, i);
                        //data()->Component(axis)(j) = from->Component(axis)(i);
                    }
                }
            }
        }

    template <class TV> void FaceArray<TV>::
        Update_Regions(
                T_FACE_ARRAY* updated,
                std::vector<FaceArray* > &parts,
                ::nimbus::GeometricRegion &region,
                int o_x = 0,
                int o_y = 0) {
            T_BOX box;
            unsigned int parts_num = parts.size();
            for (unsigned int i = 0; i < parts_num; i++)
            {
                FaceArray *part = parts[i];
                ::nimbus::GeometricRegion *r = part->region();
                // allow update to only covered regions
                if (region.Covers(r)) {
                    box = T_BOX(
                            r->x() + o_x,
                            r->x()+r->dx() + o_x -1,
                            r->y() + o_y,
                            r->y()+r->dy() + o_y -1);
                    part->Update_Face_Array(updated, &box);
                }
            }
        }

    template class ::water_app_data::FaceArray<TVF2>;

} // namespace water_app_data
