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
 * Helper functions in water_driver.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "./water_app.h"
#include "./water_data_types.h"
#include "./water_driver.h"
#include "./physbam_include.h"

template <class TV, class T> void NonAdvData<TV, T>::
Initialize_Phi()
{
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    ARRAY<T,TV_INT> &phi = particle_levelset_evolution->phi;
    for (typename GRID<TV>::CELL_ITERATOR iterator(*grid);
            iterator.Valid(); iterator.Next())
    {
        const TV &X = iterator.Location();
        phi(iterator.Cell_Index()) = X.y - (T)grid->min_dX*5;
    }
}

template <class TV, class T> void NonAdvData<TV, T>::
Set_Boundary_Conditions(
        WaterDriver<TV> *driver,
        const T time,
        FaceArray<TV> *face_velocities)
{
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    projection->elliptic_solver->psi_D.Fill(false);
    projection->elliptic_solver->psi_N.Fill(false);

    for(int axis = 1; axis <= TV::dimension; axis++)
    {
        for(int axis_side =1; axis_side <= 2; axis_side++)
        {
            int side = 2*(axis-1) + axis_side;
            TV_INT interior_cell_offset = axis_side==1?
                TV_INT() : -TV_INT::Axis_Vector(axis);
            TV_INT exterior_cell_offset = axis_side==1?
                -TV_INT::Axis_Vector(axis) : TV_INT();
            TV_INT boundary_face_offset = axis_side==1?
                TV_INT::Axis_Vector(axis) : -TV_INT::Axis_Vector(axis);
            if ((*domain_boundary)(axis)(axis_side))
            {
                for (typename GRID<TV>::FACE_ITERATOR
                        iterator(*grid, 1, GRID<TV>::BOUNDARY_REGION, side);
                        iterator.Valid(); iterator.Next())
                {
                    TV_INT face = iterator.Face_Index() + boundary_face_offset;
                    if (particle_levelset_evolution->phi(face + interior_cell_offset) <= 0)
                    {
                        if (face_velocities->data->Component(axis).Valid_Index(face))
                        {
                            projection->elliptic_solver->psi_N.Component(axis)(face) = true;
                            face_velocities->data->Component(axis)(face)=0;
                        }
                    }
                    else
                    {
                        TV_INT cell = face + exterior_cell_offset;
                        projection->elliptic_solver->psi_D(cell)=true;
                        projection->p(cell)=0;
                    }
                }
            }
            else for (typename GRID<TV>::FACE_ITERATOR
                    iterator(*grid, 1, GRID<TV>::BOUNDARY_REGION, side);
                    iterator.Valid(); iterator.Next())
            {
                TV_INT cell = iterator.Face_Index() + interior_cell_offset;
                projection->elliptic_solver->psi_D(cell) = true;
                projection->p(cell)=0;
            }
        }
    }
    for (typename GRID<TV>::FACE_ITERATOR iterator(*grid);
            iterator.Valid(); iterator.Next())
    {
        for(int i=1; i <= sources->m; i++)
        {
            if(time<=3 && (*sources)(i)->Lazy_Inside(iterator.Location()))
            {
                projection->elliptic_solver->
                    psi_N(iterator.Full_Index()) = true;
                if((TV::dimension==2 && iterator.Axis() == 1) ||
                        (TV::dimension==3 && iterator.Axis() == 3))
                    (*face_velocities->data)(iterator.Full_Index()) = -1;
                else
                    (*face_velocities->data)(iterator.Full_Index()) = 0;
            }
        }
    }
}

template <class TV, class T> void NonAdvData<TV, T>::
Write_Output_Files_EF(const int frame, FaceArray<TV>* face_velocities,
    WaterDriver<TV> * driver) {
  if(!driver->write_output_files) return;
  std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      driver->output_directory+"/"+f+"/mac_velocities",
      *(face_velocities->data));
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      driver->output_directory+"/common/grid", *grid);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      driver->output_directory+"/"+f+"/pressure", (*incompressible).projection.p);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      driver->output_directory+"/"+f+"/psi_N", (*projection).elliptic_solver->psi_N);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      driver->output_directory+"/"+f+"/psi_D", (*projection).elliptic_solver->psi_D);
  // T_PARTICLE_LEVELSET& particle_levelset = (*particle_levelset_evolution).particle_levelset;
//  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
//  FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
//  FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
//  FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
//  FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
//  FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
//#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
//#else
//  PHYSBAM_FATAL_ERROR("Cannot
//      read
//      doubles");
//#endif
//


}

#ifndef TEMPLATE_USE
#define TEMPLATE_USE
typedef VECTOR<float, 2> TVF2;
typedef float TF;
#endif  // TEMPLATE_USE

template class NonAdvData<TVF2, TF>;
