//#####################################################################
// Copyright 2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_PATTERN
//##################################################################### 
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_PATTERN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_REGION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> FRACTURE_PATTERN<T>::
FRACTURE_PATTERN()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FRACTURE_PATTERN<T>::
~FRACTURE_PATTERN()
{}
//#####################################################################
// Function Intersect_With_Rigid_Body
//#####################################################################
template<class T> void FRACTURE_PATTERN<T>::
Intersect_With_Rigid_Body(const RIGID_BODY<TV>& body,const TV& point_of_impact,ARRAY<int>& added_bodies,const bool allow_refracture,const bool use_particle_optimization,const bool generate_object_tessellation)
{
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
    // Construct body region
    FRAME<TV> levelset_frame=FRAME<TV>();
    IMPLICIT_OBJECT<TV>* implicit_object=body.implicit_object->object_space_implicit_object;
    while(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* transform=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(implicit_object)){
        levelset_frame=levelset_frame**transform->transform;
        implicit_object=transform->object_space_implicit_object;}
    FRACTURE_REGION<T> body_region(body.simplicial_object,dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(implicit_object),use_particle_optimization);
    body_region.need_destroy_data=false;
    MATRIX<T,3> R_ed;
    body_region.levelset_T=body.Rotation().Rotate(levelset_frame.t)+body.X();
    body_region.levelset_RS=(body.Rotation()*levelset_frame.r).Rotation_Matrix();
    body_region.object_T=body.X();
    body_region.object_RS=body.Rotation().Rotation_Matrix();
    R_ed=body_region.levelset_RS*DIAGONAL_MATRIX<T,3>(body_region.implicit_object->levelset.grid.dX);
    body_region.extra_levelset_frame=levelset_frame;
    body_region.fracture_offset=body_region.implicit_object->levelset.grid.Closest_Node(body_region.levelset_RS.Inverse()*(point_of_impact-body_region.levelset_T))-TV_INT::All_Ones_Vector();
    T density=body.Mass()/body_region.Compute_Volume();

    for(int r=1;r<=regions.m;r++){
        FRACTURE_REGION<T>& region=*regions(r);
        region.levelset_RS=region.object_RS=R_ed*DIAGONAL_MATRIX<T,3>(region.implicit_object->levelset.grid.dX).Inverse();
        region.levelset_T=region.object_T=R_ed*(-region.implicit_object->levelset.grid.domain.min_corner/region.implicit_object->levelset.grid.dX+
            TV(body_region.fracture_offset-region.fracture_offset))+(body_region.levelset_RS*body_region.implicit_object->levelset.grid.domain.min_corner)+
            body_region.levelset_T;
        ARRAY<FRACTURE_REGION<T>*> new_regions=body_region.Intersect_With_Rigid_Body(*regions(r),use_particle_optimization,generate_object_tessellation);
        if(!new_regions.m) continue;
        for(int i=1;i<=new_regions.m;i++){
            RIGID_BODY<TV>* new_body=new RIGID_BODY<TV>(body.rigid_body_collection,true);
            // Initialize frame, mass, and inertia. coefficient of friction, coefficient of restitution
            new_body->Set_Coefficient_Of_Friction(body.coefficient_of_friction);new_body->Set_Coefficient_Of_Restitution(body.coefficient_of_restitution);
            T_WORLD_SPACE_INERTIA_TENSOR inertia;
            new_regions(i)->extra_levelset_frame=levelset_frame;
            new_regions(i)->Compute_Inertial_Properties(density,new_body->X(),new_body->Mass(),inertia);
            new_body->Diagonalize_Inertia_Tensor(inertia);
            new_body->Add_Structure(*new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new_regions(i)->implicit_object,true,new FRAME<TV>(new_body->Frame().Inverse_Times(levelset_frame))));
            new_body->Add_Structure(*new_regions(i)->triangulated_surface);
            for(int p=1;p<=new_body->simplicial_object->particles.array_collection->Size();p++)
                new_body->simplicial_object->particles.X(p)=new_body->Object_Space_Point(new_body->simplicial_object->particles.X(p));
            new_regions(i)->need_destroy_data=false;
            delete new_regions(i);
            
            new_body->Set_Frame(body.Frame()*new_body->Frame());
            new_body->Angular_Velocity()=body.Angular_Velocity();
            new_body->V()=body.V()+VECTOR<T,3>::Cross_Product(body.Angular_Velocity(),new_body->X()-body.X());
            new_body->Update_Angular_Momentum();
            new_body->simplicial_object->Initialize_Hierarchy();
            new_body->simplicial_object->Update_Bounding_Box();
            new_body->Update_Bounding_Box();
            if(allow_refracture) new_body->fracture_threshold=2*body.fracture_threshold;
            added_bodies.Append(body.rigid_body_collection.Add_Rigid_Body_And_Geometry(new_body));}}
}
//#####################################################################
// Function Intersect_With_Rigid_Body
//#####################################################################
template<class T> void FRACTURE_PATTERN<T>::
Intersect_With_Rigid_Body(const RIGID_BODY<VECTOR<T,1> >& body,const VECTOR<T,1>& point_of_impact,ARRAY<int>& added_bodies,const bool allow_refracture,const bool use_particle_optimization,const bool generate_object_tessellation)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Intersect_With_Rigid_Body
//#####################################################################
template<class T> void FRACTURE_PATTERN<T>::
Intersect_With_Rigid_Body(const RIGID_BODY<VECTOR<T,2> >& body,const VECTOR<T,2>& point_of_impact,ARRAY<int>& added_bodies,const bool allow_refracture,const bool use_particle_optimization,const bool generate_object_tessellation)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void FRACTURE_PATTERN<T>::
Read(TYPED_ISTREAM& input)
{
    int region_count;
    Read_Binary(input,region_count);
    regions.Resize(region_count);
    for(int r=1;r<=region_count;r++){
        regions(r)=new FRACTURE_REGION<T>(0,0,false);
        regions(r)->Read(input);
        if(use_particle_partitions) regions(r)->Initialize_Particle_Partition();}
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void FRACTURE_PATTERN<T>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,regions.m);
    for(int r=1;r<=regions.m;r++) regions(r)->Write(output);
}
//#####################################################################
template class FRACTURE_PATTERN<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_PATTERN<double>;
#endif
