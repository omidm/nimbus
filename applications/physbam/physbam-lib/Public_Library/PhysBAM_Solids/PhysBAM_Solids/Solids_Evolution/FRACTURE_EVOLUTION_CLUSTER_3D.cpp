//#####################################################################
// Copyright 2006-2007, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_BODY_FRACTURE_OBJECT_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION_CLUSTER_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
using namespace PhysBAM;

//#####################################################################
// Function Process_Fracture
//#####################################################################
template<class TV> void FRACTURE_EVOLUTION_CLUSTER_3D<TV>::
Process_Cluster_Fracture(const T dt,const T time,SOLIDS_EVOLUTION<TV>* solids_evolution,SOLIDS_EVOLUTION_CALLBACKS<TV>* solids_evolution_callbacks,bool force)
{
    PHYSBAM_FATAL_ERROR("RIGID_BODY_CLUSTER_3D is obsolete; use RIGID_BODY_CLUSTER_BINDING instead.");
#if 0
    if(force){
        {std::stringstream ss;ss<<"Preprocessing CLUSTER fracture"<<std::endl;LOG::filecout(ss.str());}
        
        for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i)){
            if(RIGID_BODY_CLUSTER_3D<T>* rigid_body_cluster=dynamic_cast<RIGID_BODY_CLUSTER_3D<T>*>(&rigid_body_particles.Rigid_Body(i))){
                if(rigid_body_cluster->perform_cluster_breaks){
                    RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>& impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>&>(*rigid_body_cluster->impulse_accumulator);
                    if(max(impulse_accumulator.accumulated_impulse.linear.Magnitude(),impulse_accumulator.accumulated_impulse.angular.Magnitude())>rigid_body_cluster->fracture_threshold
                        || rigid_body_cluster->force_cluster_break_checks)
                        clusters_to_check.Append(rigid_body_cluster);}}}
        
        if(clusters_to_check.m>0){
            {std::stringstream ss;ss<<"---------------- Max impulse above threshold!!"<<std::endl;LOG::filecout(ss.str());}
            for(int c=1;c<=clusters_to_check.m;c++) clusters_to_check(c)->Save_Constituent_Bodies_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_NEW_STATE,time+dt);
            for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i)) rigid_body_particles.Rigid_Body(i).Restore_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_OLD_STATE);
            // remove all tagged clusters -- ones that had a greater impulse
            for(int c=1;c<=clusters_to_check.m;c++){
                clusters_to_check(c)->Save_Constituent_Bodies_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_OLD_STATE,time);
                rigid_body_particles.Remove_Cluster_Body(clusters_to_check(c)->Id_Number(),false);
                clusters_to_check(c)->Pre_Strain_Calculation();}
        
            LOG::Time("advancing rigid bodies in Postprocess_Fracture");
            // evolve bodies so that we can determine strain between them; set collision iterations to use less during strain testing
            solids_evolution->Get_Current_Kinematic_Keyframes(dt,time);
            int old_collision_iterations=solids_parameters.rigid_body_collision_parameters.collision_iterations;
            solids_evolution_callbacks->Pre_Advance_Cluster_Fracture(dt,time);
            solids_parameters.rigid_body_collision_parameters.collision_iterations=collision_iterations;
            solids_parameters.solids_evolution->Get_Current_Kinematic_Keyframes(dt,time);
            PHYSBAM_FATAL_ERROR("RIGID_DEFORMABLE_EVOLUTION_OLD is gone.");
            solids_parameters.rigid_body_collision_parameters.collision_iterations=old_collision_iterations;
            solids_evolution_callbacks->Post_Advance_Cluster_Fracture(dt,time);
        
            for(int i=1;i<=clusters_to_check.m;i++){
                clusters_to_check(i)->Post_Strain_Calculation();
                if(clusters_to_check(i)->Create_New_Clusters_Based_On_Strain(solids_parameters.collision_body_list,time,time-previous_time)) delete clusters_to_check(i);}
            solids_parameters.solid_body_collection.deformable_object.Update_Simulated_Particles();
            
            // move all bodies back to old position so that we can evolve using new clusters
            for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i)) rigid_body_particles.Rigid_Body(i).Restore_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_OLD_STATE);
            solids_parameters.solids_evolution->Get_Current_Kinematic_Keyframes(dt,time);
            rigid_body_particles.Reset_Impulse_Accumulators();
            PHYSBAM_FATAL_ERROR("RIGID_DEFORMABLE_EVOLUTION_OLD is gone.");
            for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i))
                rigid_body_particles.Rigid_Body(i).Save_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_NEW_STATE,time+dt);}
        clusters_to_check.Remove_All();
        previous_time=time;}
    else{
        {std::stringstream ss;ss<<"Preprocessing CLUSTER fracture"<<std::endl;LOG::filecout(ss.str());}
        for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i)){
            if(RIGID_BODY_CLUSTER_3D<T>* rigid_body_cluster=dynamic_cast<RIGID_BODY_CLUSTER_3D<T>*>(&rigid_body_particles.Rigid_Body(i))){
                if(rigid_body_cluster->perform_cluster_breaks){
                    RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>& impulse_accumulator=dynamic_cast<RIGID_BODY_IMPULSE_ACCUMULATOR<TV,3>&>(*rigid_body_cluster->impulse_accumulator);
                    if(max(impulse_accumulator.accumulated_impulse.linear.Magnitude(),impulse_accumulator.accumulated_impulse.angular.Magnitude())>rigid_body_cluster->fracture_threshold
                        || rigid_body_cluster->force_cluster_break_checks)
                        clusters_to_check.Append(rigid_body_cluster);}}}
        
        if(clusters_to_check.m>0){
            {std::stringstream ss;ss<<"---------------- Max impulse above threshold!!"<<std::endl;LOG::filecout(ss.str());}
        
            for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i)) rigid_body_particles.Rigid_Body(i).Save_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_OLD_STATE);
            // remove all tagged clusters -- ones that had a greater impulse
            for(int c=1;c<=clusters_to_check.m;c++){
                clusters_to_check(c)->Save_Constituent_Bodies_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_OLD_STATE,time);
                rigid_body_particles.Remove_Cluster_Body(clusters_to_check(c)->Id_Number(),false);
                clusters_to_check(c)->Pre_Strain_Calculation();}
        
            LOG::Time("advancing rigid bodies in Postprocess_Fracture");
            // evolve bodies so that we can determine strain between them; set collision iterations to use less during strain testing
            solids_parameters.solids_evolution->Get_Current_Kinematic_Keyframes(dt,time);
            int old_collision_iterations=solids_parameters.rigid_body_collision_parameters.collision_iterations;
            solids_evolution_callbacks->Pre_Advance_Cluster_Fracture(dt,time);
            solids_parameters.rigid_body_collision_parameters.collision_iterations=collision_iterations;
            // TODO: merge 
            PHYSBAM_FATAL_ERROR("Must implement this using new evolution");
            solids_parameters.rigid_body_collision_parameters.collision_iterations=old_collision_iterations;
            solids_evolution_callbacks->Post_Advance_Cluster_Fracture(dt,time);
        
            for(int i=1;i<=clusters_to_check.m;i++){clusters_to_check(i)->Post_Strain_Calculation();if(clusters_to_check(i)->Create_New_Clusters_Based_On_Strain(solids_parameters.collision_body_list,time,time-previous_time)) delete clusters_to_check(i);}
            
            // move all bodies back to old position so that we can evolve using new clusters
            for(int i(1);i<=rigid_body_particles.array_collection->Size();i++) if(rigid_body_particles.Is_Active(i)) rigid_body_particles.Rigid_Body(i).Restore_State(COLLISION_BODY<TV>::SOLIDS_EVOLUTION_RIGID_BODY_OLD_STATE);}
        clusters_to_check.Remove_All();
        previous_time=time;}
#endif
}
//#####################################################################
template class FRACTURE_EVOLUTION_CLUSTER_3D<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_EVOLUTION_CLUSTER_3D<VECTOR<double,3> >;
#endif

