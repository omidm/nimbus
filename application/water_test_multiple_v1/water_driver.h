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
 * Include job and function definitions here. Also, any class or struct
 * definitions to group data used by a water simulation can be included
 * here.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_WATER_DRIVER_H_
#define NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_WATER_DRIVER_H_

/* Include relevant PhysBAM files here.
*/
#include "shared/nimbus.h"
#include "physbam_include.h"

#define driver_id 15

using namespace PhysBAM;

/* This class contains all constant parameters and policies, and functions to
 * operate on data supplied by Nimbus.
 */
template <class TV>
class WaterDriver : public LEVELSET_CALLBACKS<GRID<TV> >, public Data
{
    private:
        /* typedefs */
        typedef typename TV::SCALAR T;
        typedef typename TV::template REBIND<int>::TYPE TV_INT;

        typedef typename LEVELSET_POLICY<GRID<TV> >::
            LEVELSET T_LEVELSET;
        typedef typename ADVECTION_POLICY<GRID<TV> >::
            ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
        typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::
            ARRAYS_SCALAR T_ARRAYS_SCALAR;
        typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::
            FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
        typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::
            TYPE T_ARRAYS_BOOL;
        typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::
            TYPE T_FACE_ARRAYS_BOOL;

        int size_;
        int frame_done;
        T target_time;

    public:

        Data *face_velocities;
        Data *sim_data;

        int id_debug;

        WaterDriver(const STREAM_TYPE stream_type_input);
        virtual ~WaterDriver() {}

        virtual void create();
        virtual Data* clone();
        virtual int get_debug_info();

        /* water simulation parameters
        */
        //TODO(eventually): time and frame information should not be here
        STREAM_TYPE stream_type;
        int number_of_ghost_cells;
        T cfl, initial_time, time, dt, frame_rate;
        int first_frame, last_frame, current_frame;
        int output_number;
        std::string frame_title;

        /* I/O, logging
        */
        int write_substeps_level;
        bool write_output_files;
        std::string output_directory;

        /* water driver functions, these should be called from the execute
         * functions for the jobs
         */
        bool CheckProceed();
        void IncreaseTime();
        bool IsFrameDone();

        /* helper functions.
        */
        void Write_Substep
            (const std::string &title, const int substep, const int level = 0);
        void Write_Output_Files(const int frame);

        T Time_At_Frame(const int frame) const
        {
            return initial_time + (frame-first_frame)/frame_rate;
        }

        /* callback functions needed here
        */
        void Adjust_Particle_For_Domain_Boundaries(
                PARTICLE_LEVELSET_PARTICLES<TV> &particles,
                const int index,TV &V,
                const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
                const T dt, const T time);

        void Get_Levelset_Velocity(
                const GRID<TV> &grid,
                T_LEVELSET &levelset,
                ARRAY<T,FACE_INDEX<TV::dimension> > &V_levelset,
                const T time) const PHYSBAM_OVERRIDE;

        void Get_Levelset_Velocity(
                const GRID<TV> &grid,
                LEVELSET_MULTIPLE_UNIFORM<GRID<TV> > &levelset_multiple,
                ARRAY<T, FACE_INDEX<TV::dimension> > &V_levelset,
                const T time) const PHYSBAM_OVERRIDE {}
};

#endif  // NIMBUS_APPLICATION_WATER_TEST_MULTIPLE_V1_WATER_DRIVER_H_
