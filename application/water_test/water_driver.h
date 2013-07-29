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
 */

#ifndef NIMBUS_APPLICATION_WATER_TEST_WATER_DRIVER_H_
#define NIMBUS_APPLICATION_WATER_TEST_WATER_DRIVER_H_

/* Include relevant PhysBAM files here.
 */
#include "./physbam_include.h"
#include "./water_data_types.h"

using namespace PhysBAM;    // NOLINT

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
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    public:

    WaterDriver();
    virtual ~WaterDriver();

    /* water simulation data
     */
    Grid<TV> mac_grid;
    MPIGrid<TV> mpi_grid;
    FaceArray<TV> face_velocities;
    FaceArrayGhost<TV> face_velocities_ghost;
    NonAdvData<TV, T> sim_data;

    /* water simulation parameters
     */
    int number_of_ghost_cells;
    T cfl, initial_time, time, frame_rate;
    int last_frame, current_frame, output_number;
    int write_substeps_level;
    bool write_output_files_flag;
    std::string frame_title, output_directory;

    /* water driver functions, these should be called from the execute
     * functions for the jobs
     */
    void initialize();
    void run_upto_advection();
    void run_advect();
    void run_after_advection();

    /* helper functions.
     */
    void time_at_frame(const int frame) const;
    void write_output_files(const int frame);
};

#endif  // NIMBUS_APPLICATION_WATER_TEST_WATER_DRIVER_H_
