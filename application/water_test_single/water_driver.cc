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

#include <iostream>
#include "./water_driver.h"
#include "./water_data_types.h"

using namespace PhysBAM;

#ifndef TEMPLATE_USE
#define TEMPLATE_USE
typedef VECTOR<float, 2> TVF2;
typedef VECTOR<float, 3> TVF3;
typedef float TF;
#endif  // TEMPLATE_USE

template <class TV> WaterDriver<TV> ::
WaterDriver(const STREAM_TYPE stream_type_input):
    stream_type(stream_type_input)
{
    frame_done = true;
    target_time = (T)0;

    // setup time
    initial_time = (T)0;
    first_frame = (T)0;
    time = (T)0;

    last_frame = 20;
    frame_rate = 24;
    current_frame = 0;

    // other parameters
    number_of_ghost_cells = 3;
    cfl = 0.9;

    // I/O & logging
    write_substeps_level = -1;
    write_output_files = true;
    output_directory = "output";
    frame_title = "";
    output_number = 0;

    // debugging information
    id_debug = driver_id;
};

template<class TV>
void Write_Substep_Helper
(void *writer, const std::string &title, int substep, int level)
{
    ((WaterDriver<TV> *)writer)->Write_Substep(title, substep, level);
};

template <class TV>
Data* WaterDriver<TV> :: Clone()
{
    std::cout << "Cloning waterdriver\n";
    return new WaterDriver<TV>(stream_type);
};

template <class TV>
void WaterDriver<TV> :: Create()
{
    std::cout << "Initialize water driver ...\n";

    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File
        (output_directory+"/common/log.txt", false);

    Initialize_Particles();
    Initialize_Read_Write_General_Structures();

    DEBUG_SUBSTEPS::Set_Substep_Writer((void *)this, &Write_Substep_Helper<TV>);

    std::cout << "Completed initializing water driver\n";
}

template <class TV>
int WaterDriver<TV> :: get_debug_info()
{
    return id_debug;
}

template <class TV> void WaterDriver<TV>::
Get_Levelset_Velocity(
        const GRID<TV> &grid,
        T_LEVELSET& levelset,
        ARRAY<T,FACE_INDEX<TV::dimension> > &V_levelset,
        const T time) const PHYSBAM_OVERRIDE
{
    FaceArray<TV> *fv = (FaceArray<TV> *)face_velocities;
    V_levelset = *fv->data;
}

template <class TV> void WaterDriver<TV>::
Adjust_Particle_For_Domain_Boundaries(
        PARTICLE_LEVELSET_PARTICLES<TV> &particles,
        const int index, TV &V,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt, const T time)
{
    NonAdvData<TV, T> *sd = (NonAdvData<TV, T> *)sim_data;
    sd->Adjust_Particle_For_Domain_Boundaries(particles, index, V,
            particle_type, dt, time);
}

template <class TV> void WaterDriver<TV>::
Write_Substep(
        const std::string &title,
        const int substep,
        const int level)
{
    if (level <= write_substeps_level)
    {
        frame_title = title;
        std::stringstream ss;
        ss << "Writing substep [" << frame_title << "]: output_number="
            <<output_number+1 << ", time=" << time << ", frame=" <<
            current_frame << ", substep=" << substep << std::endl;
        LOG::filecout(ss.str());
        Write_Output_Files(++output_number);
        frame_title="";
    }
}

template<class TV> void WaterDriver<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory +
            STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(output_directory + "/common");

    FILE_UTILITIES::Write_To_Text_File(output_directory +
            STRING_UTILITIES::string_sprintf("/%d/frame_title", frame),
            frame_title);

    if(frame == first_frame) 
        FILE_UTILITIES::Write_To_Text_File
            (output_directory + "/common/first_frame", frame, "\n");

    NonAdvData<TV, T> *sd = (NonAdvData<TV, T> *)sim_data;
    FaceArray<TV> *fv = (FaceArray<TV> *)face_velocities;
    sd->Write_Output_Files_EF(frame, fv, this);

    FILE_UTILITIES::Write_To_Text_File
        (output_directory + "/common/last_frame", frame, "\n");
}

template<class TV> bool WaterDriver<TV>::
CheckProceed()
{
    NonAdvData<TV, T> *sd = (NonAdvData<TV, T> *)sim_data;

    std::cout << "## Simulating frame: " << current_frame
        << ", time :" << time << "\n";

    if (frame_done)
    {
        frame_done = false;
        current_frame++;
        sd->current_frame = current_frame;
        if (current_frame > last_frame)
            return false;
        target_time = Time_At_Frame(current_frame);
    }

    LOG::Time("Calculate Dt");
    FaceArray<TV> *fv = (FaceArray<TV> *)face_velocities;

    sd->particle_levelset_evolution->Set_Number_Particles_Per_Cell(16);
    dt = cfl * sd->incompressible->CFL(*fv->data);
    dt = min(dt, sd->particle_levelset_evolution->CFL(false, false));
    if ( time + dt >= target_time)
    {
        dt = target_time-time;
        frame_done=true;
    }
    else if (time + 2*dt >= target_time)
    {
        dt = .5*(target_time - time);
    }

    sd->dt = dt;
    return true;
}

template<class TV> void WaterDriver<TV>::
IncreaseTime()
{
    time += dt;
    NonAdvData<TV, T> *sd = (NonAdvData<TV, T> *)sim_data;
    sd->time = time;
}

template<class TV> bool WaterDriver<TV>::
IsFrameDone()
{
    return frame_done;
}

template class WaterDriver<TVF2>;
template class WaterDriver<TVF3>;
