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
        phi(iterator.Cell_Index()) = X.y + X.x - (T)grid->min_dX*100;
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
Write_Output_Files_EF(const int frame, FaceArray<TV> *face_velocities,
    WaterDriver<TV> *driver) {
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
  T_PARTICLE_LEVELSET& particle_levelset = (*particle_levelset_evolution).particle_levelset;
  FILE_UTILITIES::Write_To_File(driver->stream_type, driver->output_directory+"/"+f+"/levelset",
      particle_levelset.levelset);
  FILE_UTILITIES::Write_To_File(driver->stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",
        driver->output_directory.c_str(),frame,"positive_particles"), particle_levelset.positive_particles);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      STRING_UTILITIES::string_sprintf("%s/%d/%s", driver->output_directory.c_str(),frame,"negative_particles"),
      particle_levelset.negative_particles);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      STRING_UTILITIES::string_sprintf("%s/%d/%s", driver->output_directory.c_str(),frame,"removed_positive_particles"), 
      particle_levelset.removed_positive_particles);
  FILE_UTILITIES::Write_To_File(driver->stream_type,
      STRING_UTILITIES::string_sprintf("%s/%d/%s", driver->output_directory.c_str(), frame, "removed_negative_particles"),
      particle_levelset.removed_negative_particles);
  FILE_UTILITIES::Write_To_Text_File(driver->output_directory+"/"+f+"/last_unique_particle_id",
      particle_levelset.last_unique_particle_id);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
  PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif

}

template <class TV, class T> void NonAdvData<TV, T>::
Adjust_Phi_With_Sources(const T time)
{
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    if (time>3)
        return;

    for (typename GRID<TV>::CELL_ITERATOR iterator(*grid);
            iterator.Valid(); iterator.Next())
    {
        TV_INT index=iterator.Cell_Index();
        for(int i = 1; i <= sources->m; i++)
        {
            particle_levelset_evolution->phi(index) =
                min(particle_levelset_evolution->phi(index),
                        (*sources)(i)->Extended_Phi( iterator.Location() ) );
        }
    }
}

template <class TV, class T> void NonAdvData<TV, T>::
Adjust_Particle_For_Domain_Boundaries(
        PARTICLE_LEVELSET_PARTICLES<TV> &particles,
        const int index, TV &V,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt, const T time)
{
    if (particle_type == PARTICLE_LEVELSET_POSITIVE ||
            particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE)
        return;

    TV &X = particles.X(index);
    TV X_new = X + dt*V;

    T max_collision_distance = particle_levelset_evolution->
        particle_levelset.Particle_Collision_Distance
        (particles.quantized_collision_distance(index));
    T min_collision_distance = particle_levelset_evolution->particle_levelset.
        min_collision_distance_factor * max_collision_distance;
    TV min_corner = grid->domain.Minimum_Corner(),
       max_corner = grid->domain.Maximum_Corner();

    for (int axis=1; axis <= GRID<TV>::dimension; axis++)
    {
        if((*domain_boundary)[axis][1] &&
                X_new[axis] < min_corner[axis] + max_collision_distance)
        {
            T collision_distance = X[axis] - min_corner[axis];
            if (collision_distance > max_collision_distance)
                collision_distance = X_new[axis] - min_corner[axis];
            collision_distance =
                max(min_collision_distance, collision_distance);
            X_new[axis] += max((T)0, 
                    min_corner[axis] - X_new[axis] + collision_distance);
            V[axis] = max((T)0, V[axis]);
            X = X_new -dt * V;
        }
        if ((*domain_boundary)[axis][2] &&
                X_new[axis] > max_corner[axis] - max_collision_distance)
        {
            T collision_distance = max_corner[axis] - X[axis];
            if (collision_distance > max_collision_distance)
                collision_distance = max_corner[axis] - X_new[axis];
            collision_distance =
                max(min_collision_distance, collision_distance);
            X_new[axis] -= max((T)0,
                    X_new[axis] - max_corner[axis] + collision_distance);
            V[axis] = min((T)0, V[axis]);
            X = X_new - dt * V;
        }
    }
}

#ifndef TEMPLATE_USE
#define TEMPLATE_USE
typedef VECTOR<float, 2> TVF2;
typedef float TF;
#endif  // TEMPLATE_USE

template class NonAdvData<TVF2, TF>;
