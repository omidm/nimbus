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

#include "assert.h"
#include "physbam_deserialize_data_arrays_2d.h"
#include "physbam_deserialize_data_common_2d.h"
#define INITIALIZATION_VALUE 0

namespace physbam_pb {

    void make_physbam_object(FaceArray2 *phys_fa,
            const ::communication::PhysbamFaceArray2d &pb_fa) {

        assert(phys_fa);
        phys_fa->Clean_Memory();
        
        if (pb_fa.has_domain_indices())
            make_physbam_object(
                    &phys_fa->domain_indices,
                    pb_fa.domain_indices());
        
        if (pb_fa.has_buffer_size())
            phys_fa->buffer_size = pb_fa.buffer_size();
        assert(pb_fa.buffer_size() == pb_fa.values_size());
        
        /*
         * The following code is for initializing phys_fa, by copying Resize() method in /PhysBAM/PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h
         * Now we directly call Clean_Memory() method. If find any problem with Clean_Memory(), may switch back to Resize() method.
         * 
        assert(pb_fa.has_domain_indices()); //zhihao: didn't handle the case in which pb_fa doesn't have domain_indices
        RangeI2 domain = phys_fa->domain_indices;
        int dimension = phys_fa->dimension;
        VectorRangeI2 domains;
        VI2 sizes_new;        
        for (int i = 1; i <= dimension; i++) {
            domains(i) = domain;
            domains(i).max_corner(i)++;
            sizes_new(i) = (domains(i).Edge_Lengths() + 1).Product();            
        }
        
        int buffer_size = sizes_new.Sum();
        
        float* p = new float[buffer_size];
        float* p_start = p;
        for (int i = 1; i <= dimension; i++) {
        	ArrayView2 array_new(domains(i), p_start);
            array_new.Fill(INITIALIZATION_VALUE);
            p_start += sizes_new(i);
        }
        phys_fa->base_pointer = p;

        assert(buffer_size = pb_fa.buffer_size()); //zhihao: do we have to assert this? 
        */
        
        float *buff_values = phys_fa->base_pointer;
        for (int i = 0; i < pb_fa.values_size(); i++) {
            buff_values[i] = pb_fa.values(i);
        }
    }

} // namespace physbam_pb
