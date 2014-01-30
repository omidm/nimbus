//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <set>
#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/data_include.h"
#include "data/physbam/translator_physbam.h"
#include "data/scalar_data.h"
#include "application/water_alternate_fine/water_example.h"
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include "shared/nimbus.h"
#include "worker/physical_data_instance.h"
using namespace PhysBAM;
//#####################################################################
// WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
WATER_EXAMPLE(const STREAM_TYPE stream_type_input) :
    stream_type(stream_type_input),
    initial_time(0),
    first_frame(0),
    last_frame(application::kLastFrame),
    frame_rate(24),
    write_substeps_level(-1),
    write_output_files(true),
    output_directory(application::kOutputDir),
    number_of_ghost_cells(application::kGhostNum),
    cfl(1),
    mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    projection(*new PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >(mac_grid,false,false,false,false,NULL/*thread_queue*/)),
    particle_levelset_evolution(mac_grid,number_of_ghost_cells),
    incompressible(mac_grid,projection),
    boundary(0),
    collision_bodies_affecting_fluid(mac_grid)
{
    Initialize_Particles();
    Initialize_Read_Write_General_Structures();
}
//#####################################################################
// ~WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
~WATER_EXAMPLE()
{
    delete &projection;
}
//#####################################################################
// Initialize_Phi
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Initialize_Phi()
{
    ARRAY<T,TV_INT>& phi=particle_levelset_evolution.phi;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        const TV &X=iterator.Location();
        phi(iterator.Cell_Index())=X.y-(T)0.35;}
}
//#####################################################################
// Initialize_Phi
//#####################################################################
// NOT REQUIRED IN NIMBUS
template<class TV> void WATER_EXAMPLE<TV>::
Initialize_Grid(TV_INT counts,RANGE<TV> domain)
{
    mac_grid.Initialize(counts,domain,true);
}
//#####################################################################
// Set_Boundary_Conditions
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time)
{
    projection.elliptic_solver->psi_D.Fill(false);
    projection.elliptic_solver->psi_N.Fill(false);
    for (int axis=1;axis<=TV::dimension;axis++) {
      for(int axis_side=1;axis_side<=2;axis_side++) {
        int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset =
            axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset =
            axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset =
            axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)) {
          for (typename GRID<TV>::FACE_ITERATOR iterator(
                  mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);
              iterator.Valid();
              iterator.Next()) {
                TV_INT face = iterator.Face_Index()+boundary_face_offset;
                // Access levelset here.
                if (particle_levelset_evolution.phi(
                        face+interior_cell_offset)<=0) {
                  // Access face velocity here.
                  if (face_velocities.Component(axis).Valid_Index(face)){
                    projection.elliptic_solver->psi_N.Component(axis)(face)
                        =true;
                    face_velocities.Component(axis)(face)=0;
                  }
                } else {
                  TV_INT cell=face+exterior_cell_offset;
                  projection.elliptic_solver->psi_D(cell)=true;
                  projection.p(cell)=0;
                }
          }
        } else {
          for (typename GRID<TV>::FACE_ITERATOR iterator(
                  mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);
              iterator.Valid();
              iterator.Next()) {
            TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            projection.elliptic_solver->psi_D(cell)=true;
            projection.p(cell)=0;
          }
        }
      }
    }
    for (typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);
         iterator.Valid();
         iterator.Next()) {
      for (int i=1;i<=sources.m;i++) {
        if (time<=3 && sources(i)->Lazy_Inside(iterator.Location())) {
          projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
          // Access velocity.
          if ((TV::dimension==2 && iterator.Axis()==1) ||
              (TV::dimension==3 && iterator.Axis()==3)) {
            face_velocities(iterator.Full_Index())=-1;
          } else {
            face_velocities(iterator.Full_Index())=0;
          }
        }
      }
    }
}
//#####################################################################
// Adjust_Phi_With_Sources
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Adjust_Phi_With_Sources(const T time)
{
    if(time>3) return;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        for(int i=1;i<=sources.m;i++) particle_levelset_evolution.phi(index)=min(particle_levelset_evolution.phi(index),sources(i)->Extended_Phi(iterator.Location()));}
}
//#####################################################################
// Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=particle_levelset_evolution.particle_levelset.min_collision_distance_factor*max_collision_distance;
    TV min_corner=mac_grid.domain.Minimum_Corner(),max_corner=mac_grid.domain.Maximum_Corner();
    for(int axis=1;axis<=GRID<TV>::dimension;axis++){
        if(domain_boundary[axis][1] && X_new[axis]<min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];
            if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(domain_boundary[axis][2] && X_new[axis]>max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];
            if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    std::string filename;
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Save_To_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame)
{
    PdiVector pdv;

    // mac velocities
    const std::string fvstring = std::string(APP_FACE_VEL);
    if (application::GetTranslatorData(job, fvstring, da, &pdv)
        && data_config.GetFlag(DataConfig::VELOCITY)) {
      translator.WriteFaceArray(
          &application::kDomainFaceVel, &pdv, &face_velocities);
    }
    application::DestroyTranslatorObjects(&pdv);

    // mac velocities ghost
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    if (application::GetTranslatorData(job, fvgstring, da, &pdv)
        && data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      translator.WriteFaceArray(
          &application::kDomainFaceVelGhost, &pdv, &face_velocities_ghost);
    }
    application::DestroyTranslatorObjects(&pdv);

    // pressure
    const std::string pstring = std::string(APP_PRESSURE);
    if (application::GetTranslatorData(job, pstring, da, &pdv))
        translator.WriteScalarArray(&application::kDomainPressure,
                                   &pdv,
                                   &incompressible.projection.p);
    application::DestroyTranslatorObjects(&pdv);

    // particle leveset quantities
    T_PARTICLE_LEVELSET& particle_levelset = particle_levelset_evolution.particle_levelset;

    // levelset
    const std::string lsstring = std::string(APP_PHI);
    if (application::GetTranslatorData(job, lsstring, da, &pdv)
        && data_config.GetFlag(DataConfig::LEVELSET)) {
      translator.WriteScalarArray(
          &application::kDomainPhi, &pdv, &particle_levelset.levelset.phi);
    }
    application::DestroyTranslatorObjects(&pdv);

    // positive particles
    const std::string ppstring = std::string(APP_POS_PARTICLES);
    if (application::GetTranslatorData(job, ppstring, da, &pdv)
        && data_config.GetFlag(DataConfig::POSITIVE_PARTICLE)) {
      translator.WriteParticles(
          &application::kDomainParticles, &pdv, particle_levelset, true);
    }
    application::DestroyTranslatorObjects(&pdv);

    // negative particles
    const std::string npstring = std::string(APP_NEG_PARTICLES);
    if (application::GetTranslatorData(job, npstring, da, &pdv)
        && data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE)) {
      translator.WriteParticles(
          &application::kDomainParticles, &pdv, particle_levelset, false);
    }
    application::DestroyTranslatorObjects(&pdv);

    // Removed positive particles.
    const std::string prpstring = std::string(APP_POS_REM_PARTICLES);
    if (application::GetTranslatorData(job, prpstring, da, &pdv)
        && data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
      translator.WriteRemovedParticles(
          &application::kDomainParticles, &pdv, particle_levelset, true);
    }
    application::DestroyTranslatorObjects(&pdv);

    // Removed negative particles.
    const std::string nrpstring = std::string(APP_NEG_REM_PARTICLES);
    if (application::GetTranslatorData(job, nrpstring, da, &pdv)
        && data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
      translator.WriteRemovedParticles(
          &application::kDomainParticles, &pdv, particle_levelset, false);
    }
    application::DestroyTranslatorObjects(&pdv);

    // last unique particle id
    const std::string lupistring = std::string(APP_LAST_UNIQUE_PARTICLE_ID);
    if (Data *d = application::GetFirstData(lupistring, da)) {
        nimbus::ScalarData<int> *sd = static_cast<nimbus::ScalarData<int> * >(d);
        sd->set_scalar(particle_levelset.last_unique_particle_id);
    }
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Load_From_Nimbus(const nimbus::Job *job, const nimbus::DataArray &da, const int frame)
{
    PdiVector pdv;

    // mac velocities
    const std::string fvstring = std::string(APP_FACE_VEL);
    if (application::GetTranslatorData(job, fvstring, da, &pdv)
        && data_config.GetFlag(DataConfig::VELOCITY)) {
      translator.ReadFaceArray(
          &application::kDomainFaceVel, &pdv, &face_velocities);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating velocity.\n");

    // mac velocities
    const std::string fvgstring = std::string(APP_FACE_VEL_GHOST);
    if (application::GetTranslatorData(job, fvgstring, da, &pdv)
        && data_config.GetFlag(DataConfig::VELOCITY_GHOST)) {
      translator.ReadFaceArray(
          &application::kDomainFaceVelGhost, &pdv, &face_velocities_ghost);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating ghost velocity.\n");

    // pressure
    const std::string pstring = std::string(APP_PRESSURE);
    if (application::GetTranslatorData(job, pstring, da, &pdv))
        translator.ReadScalarArray(&application::kDomainPressure,
                                   &pdv,
                                   &incompressible.projection.p);
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating pressure.\n");

    // particle leveset quantities
    T_PARTICLE_LEVELSET& particle_levelset = particle_levelset_evolution.particle_levelset;

    // levelset
    const std::string lsstring = std::string(APP_PHI);
    if (application::GetTranslatorData(job, lsstring, da, &pdv)
        && data_config.GetFlag(DataConfig::LEVELSET)) {
      translator.ReadScalarArray(
          &application::kDomainPhi, &pdv, &particle_levelset.levelset.phi);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating levelset.\n");

    // positive particles
    const std::string ppstring = std::string(APP_POS_PARTICLES);
    if (application::GetTranslatorData(job, ppstring, da, &pdv)
        && data_config.GetFlag(DataConfig::POSITIVE_PARTICLE)) {
      translator.ReadParticles(
          &application::kDomainParticles, &pdv, particle_levelset, true);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating positive particles.\n");

    // negative particles
    const std::string npstring = std::string(APP_NEG_PARTICLES);
    if (application::GetTranslatorData(job, npstring, da, &pdv)
        && data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE)) {
      translator.ReadParticles(
          &application::kDomainParticles, &pdv, particle_levelset, false);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating negative particles.\n");

    // Removed positive particles.
    const std::string prpstring = std::string(APP_POS_REM_PARTICLES);
    if (application::GetTranslatorData(job, prpstring, da, &pdv)
        && data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
      translator.ReadRemovedParticles(
          &application::kDomainParticles, &pdv, particle_levelset, true);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating remove positive particles.\n");

    // Removed negative particles.
    const std::string nrpstring = std::string(APP_NEG_REM_PARTICLES);
    if (application::GetTranslatorData(job, nrpstring, da, &pdv)
        && data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE)) {
      translator.ReadRemovedParticles(
          &application::kDomainParticles, &pdv, particle_levelset, false);
    }
    application::DestroyTranslatorObjects(&pdv);
    dbg(APP_LOG, "Finish translating remove negative particles.\n");

    // last unique particle id
    const std::string lupistring = std::string(APP_LAST_UNIQUE_PARTICLE_ID);
    if (Data *d = application::GetFirstData(lupistring, da)) {
        nimbus::ScalarData<int> *sd = static_cast<nimbus::ScalarData<int> * >(d);
        particle_levelset.last_unique_particle_id = sd->scalar();
    }
    dbg(APP_LOG, "Finish translating particle id.\n");
}
//#####################################################################
template class WATER_EXAMPLE<VECTOR<float,3> >;
