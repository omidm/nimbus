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
 * Implementation of functions for data classes used by application.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "shared/nimbus.h"
#include "physbam_include.h"
#include "water_app_data.h"

namespace water_app_data {

    template <class TV> FaceArray<TV>::
        FaceArray(int size) {
            id_debug = face_array_id;
            this->size_ = size;
            grid = NULL;
            data = NULL;
        }

    template <class TV> void FaceArray<TV>::
        create() {
            std::cout << "Creating FaceArray\n";
            grid = new T_GRID(TV_INT::All_Ones_Vector()*size_,
                    RANGE<TV>::Unit_Box(), true);
            assert(grid);
            data = new  T_FACE_ARRAY();
            assert(data);
            data->Resize(*grid);
        }

    template <class TV> ::nimbus::Data* FaceArray<TV>::
        clone() {
            std::cout << "Cloning facearray\n";
            return new FaceArray<TV>(size_);
        }

    template <class TV> int FaceArray<TV>::
        get_debug_info() {
            return id_debug;
        }

    namespace {
        typedef ::PhysBAM::VECTOR<float, 2> TVF2;
    } // namespace
    template class ::water_app_data::FaceArray<TVF2>;

} // namespace water_app_data
