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
 * Data types used by the application jobs and functions.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_TEST_WATER_DATA_TYPES_H_
#define NIMBUS_APPLICATION_WATER_TEST_WATER_DATA_TYPES_H_

/* Include relevant PhysBAM files here.
*/
#include "shared/nimbus.h"
#include "./physbam_include.h"
#include "./water_driver.h"

#define face_array_id 20
#define face_array_ghost_id 25
#define non_adv_id 30

/* WATER_EXAMPLE is structured as follows (with the equivalent here shown in
 * brackets):
 * mac_grid (equivalent *corresponding data.grid)
 * face_velocities (*face_velocities.data)
 * boundary_scalar (*sim_data.boundary_scalar)
 * boundary (sim_data.boundary)
 * phi_boundary (sim_data.phi_boundary)
 * phi_boundary_water (*sim_data.phi_boundary_water)
 * domain_boundary (*sim_data.domain_boundary)
 * sources (*sim_data.sources)
 * particle_levelset_evolution (*sim_data.particle_levelset_evolution)
 * advection_scalar (*sim_data.advection_scalar)
 * collision_bodies_affecting_fluid (*sim_data.collision_bodies_affecting_fluid)
 * projection (*sim_data.projection)
 * incompressible (*sim_data.incompressible)
 */

using namespace PhysBAM;
using nimbus::Data;

template <class TV, class T>
class NonAdvData;

/* Face array for storing quantities like face velocities.
*/
template <class TV>
class FaceArray : public Data {
    private:
        int size_;
        typedef typename TV::SCALAR T;
        typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;

    public:

        void Advection (WaterDriver<TV> *driver,
             NonAdvData<TV, T> *sim_data,
             const T time_target);

        int id_debug;

        FaceArray(int size);
        virtual void create();
        virtual Data* clone();
        virtual int get_debug_info();

        // physbam structures and methods
        GRID<TV> *grid;
        ARRAY<T, FACE_INDEX<TV::dimension> > *data;
};

/* Ghost face array for storing scalar quantities.
*/
template <class TV>
class FaceArrayGhost : public Data {
    private:
        int size_;
    public:

        int id_debug;

        FaceArrayGhost(int size);
        virtual void create();
        virtual Data* clone();
        virtual int get_debug_info();

        // physbam structures and methods
        GRID<TV> *grid;
        typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS *data;
};

/* Add all other data used by water simulation here.  DO NOT add scalar
 * values. Scalar values can be passed around directly as parameters.
 */
template <class TV, class T>
class NonAdvData : public Data {

  typedef typename TV::template REBIND<int>::TYPE TV_INT;
  typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
  typedef typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
  typedef typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::BOUNDARY_PHI_WATER T_BOUNDARY_PHI_WATER;
  typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;


    typedef typename ADVECTION_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;




    private:
        int size_;
    public:

        int id_debug;

        NonAdvData(int size);
        virtual void create();
        virtual Data* clone();
        virtual int get_debug_info();

        bool initialize
            (WaterDriver<TV> *driver,
             FaceArray<TV> *face_velocities,
             const int frame);

        void BeforeAdvection (WaterDriver<TV> *driver,
             FaceArray<TV> *face_velocities,
             const T time_target); 

        void AfterAdvection (WaterDriver<TV> *driver,
             FaceArray<TV> *face_velocities,
             const T time_target);

        // physbam structures and methods

        // TODO(chinmayee): need to have this as a parameter
        int number_of_ghost_cells;
        T time, dt;
        int current_frame;

        GRID<TV> *grid;

        // boundary information
        BOUNDARY_UNIFORM<GRID<TV>, T>
            *boundary_scalar,
            *boundary,
            *phi_boundary;
        typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::
            BOUNDARY_PHI_WATER *phi_boundary_water;
        VECTOR<VECTOR<bool, 2>, TV::dimension> *domain_boundary;

        // sources
        ARRAY<IMPLICIT_OBJECT<TV>*> *sources;

        // fluid data
        PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> > *particle_levelset_evolution;
        ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>, T> *advection_scalar;

        // collision geometry
        typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::
            GRID_BASED_COLLISION_GEOMETRY *collision_bodies_affecting_fluid;

        // other containers
        PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > *projection;
        INCOMPRESSIBLE_UNIFORM<GRID<TV> > *incompressible;

        // helper methods
        void Initialize_Phi();
        void Set_Boundary_Conditions(
                WaterDriver<TV> *driver,
                const T time,
                FaceArray<TV> *face_velocities);
        void Adjust_Phi_With_Sources(const T time);
        void Adjust_Particle_For_Domain_Boundaries(
                PARTICLE_LEVELSET_PARTICLES<TV> &particles,
                const int index,
                TV &V,
                const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
                const T dt,
                const T time); 

        void Write_Output_Files_EF(const int frame,
                FaceArray<TV>* face_velocities,
                WaterDriver<TV> * driver);
        void Read_Output_Files_EF(const int frame) {}
};

#endif  // NIMBUS_APPLICATION_WATER_TEST_WATER_DATA_TYPES_H_
