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

#ifndef NIMBUS_APPLICATION_WATER_TEST_WATER_DRIVER_H_
#define NIMBUS_APPLICATION_WATER_TEST_WATER_DRIVER_H_

/* Include relevant PhysBAM files here.
 */
#include "./physbam_include.h"
#include "./water_data_types.h"

using namespace PhysBAM;

/* This is more like WATER_EXAMPLE.h than WATER_DRIVER.h from the original
 * PhysBAM project Water, in the sense that it directly contains all the data
 * that methods in WATER_DRIVER.cpp are operating on, rather than accessing the
 * data through driver->example. However, each machine will launch its own copy
 * of WaterDriver, like in the PhysBAM project. This initialization should
 * happen in load() in WaterApp, after which the job and data maps should be
 * built.
 */
template <class TV>
class WaterDriver : public LEVELSET_CALLBACKS<GRID<TV> >,
    public RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>
{
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

    public:

    WaterDriver();
    virtual ~WaterDriver();

    /* water simulation data
     */
    Grid<TV> mac_grid(int size);
    MPIGrid<TV> mpi_grid(int size);
    FaceArray<TV> face_velocities(int size);
    FaceArrayGhost<TV> face_velocities_ghost(int size);
    NonAdvData<TV, T> sim_data(int size);

    /* water simulation parameters
     */
    int number_of_ghost_cells;
    T cfl, initial_time, time, frame_rate;
    int first_frame, last_frame, current_frame, output_number;
//    STREAM_TYPE stream_type;
    int write_substeps_level;
    bool write_output_files_flag;
    std::string frame_title, output_directory;

    /* water driver functions, these should be called from the execute
     * functions for the jobs
     */
    void initialize(bool distributed);
    void run_upto_advection();
    void run_advect();
    void run_after_advection();

    /* helper functions.
     */
    void Write_Output_Files(const int frame);
    void Read_Output_Files(const int frame);

    T Time_At_Frame(const int frame) const
    {
        return initial_time + (frame-first_frame)/frame_rate;
    }

//    bool Set_Kinematic_Positions(
//            FRAME<TV> &frame,
//            const T time,
//            const int id)
//    {
//        T range = 0.6;
//        frame.t = TV::All_Ones_Vector()*0.5;
//        if(time<=2)
//            frame.t(2) = time*range+(1-range)/2.;
//        return true;
//    }
//
//    void Get_Levelset_Velocity(
//            const GRID<TV> &grid,
//            T_LEVELSET& levelset,
//            ARRAY<T, FACE_INDEX<TV::dimension> > &V_levelset,
//            const T time) const PHYSBAM_OVERRIDE 
//    {
//        V_levelset = *face_velocities.data;
//    }
//
//    void Get_Levelset_Velocity(
//            const GRID<TV> &grid,
//            LEVELSET_MULTIPLE_UNIFORM<GRID<TV> > &levelset_multiple,
//            ARRAY<T,FACE_INDEX<TV::dimension> > &V_levelset,
//            const T time) const PHYSBAM_OVERRIDE {}
//
//    void Adjust_Particle_For_Domain_Boundaries(
//            PARTICLE_LEVELSET_PARTICLES<TV> &particles,
//            const int index,
//            TV &V,
//            const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
//            const T dt,
//            const T time);
//    void Initialize_Grid(TV_INT counts, RANGE<TV> range);
//    void Set_Boundary_Conditions(const T time);
//    void Adjust_Phi_With_Sources(const T time);
//    void Adjust_Phi_With_Objects(const T time);
//    void Extrapolate_Phi_Into_Objects(const T time);
//    void Initialize_Phi();
};

#endif  // NIMBUS_APPLICATION_WATER_TEST_WATER_DRIVER_H_
