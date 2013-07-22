//#####################################################################
// Copyright 2006, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/RIGID_BODY_FRACTURE_OBJECT_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/FRACTURE_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
template<class TV> void FRACTURE_EVOLUTION<TV>::
Postprocess_Solids_Substep(const T time,const int substep)
{
    if(perform_fracture && Fracture_Where_High_Stress()) fractured_after_rebuild_topology=true;
    if(Add_Scripted_Cuts(time)) fractured_after_rebuild_topology=true;
    {std::stringstream ss;ss<<"XXX Postprocess_Solids_Substep: fractured after rebuild topology "<<fractured_after_rebuild_topology<<std::endl;LOG::filecout(ss.str());}
    if(fractured_after_rebuild_topology && substep%substeps_before_rebuild==0){{std::stringstream ss;ss<<"XXX rebuild topology"<<std::endl;LOG::filecout(ss.str());}Rebuild_Topology();fractured_after_rebuild_topology=false;}

    if(push_out){  // experimental push out code
        {std::stringstream ss;ss << std::endl;LOG::filecout(ss.str());}{std::stringstream ss;ss << "AFTER FRACTURE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;LOG::filecout(ss.str());}
        bool found_intersections=true;int attempts=0;
        while(found_intersections){
            found_intersections=false;{std::stringstream ss;ss << "ATTEMPT........................." << std::endl;LOG::filecout(ss.str());}
            if(Adjust_Nodes_For_Segment_Triangle_Intersections() > 0) found_intersections=true;
            if(++attempts == 1) found_intersections=false;}}
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void FRACTURE_EVOLUTION<TV>::
Postprocess_Frame(const int frame)
{
}
//#####################################################################
template class FRACTURE_EVOLUTION<VECTOR<float,1> >;
template class FRACTURE_EVOLUTION<VECTOR<float,2> >;
template class FRACTURE_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_EVOLUTION<VECTOR<double,1> >;
template class FRACTURE_EVOLUTION<VECTOR<double,2> >;
template class FRACTURE_EVOLUTION<VECTOR<double,3> >;
#endif

