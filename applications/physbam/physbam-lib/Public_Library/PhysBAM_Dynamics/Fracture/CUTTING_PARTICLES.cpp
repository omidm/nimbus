#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
CUTTING_PARTICLES::
CUTTING_PARTICLES()
{
}
//#####################################################################
// Destructor
//#####################################################################
CUTTING_PARTICLES::
~CUTTING_PARTICLES()
{
}
//#####################################################################
// Function Preallocate
//#####################################################################
void CUTTING_PARTICLES::
Preallocate(const int size)
{
    tet_node_indices.Preallocate(size);
    intersection_indices.Preallocate(size);
    particle_ids_types.Preallocate(size);
}
//#####################################################################
// Function Add_Tet_Node_Id
//#####################################################################
void CUTTING_PARTICLES::
Add_Tet_Node_Id(const int tet_node_index)
{
    int index=tet_node_indices.Append(tet_node_index);
    intersection_indices.Append(0);
    particle_ids_types.Append(TET_NODE_ID);
    tet_node_to_particle_id.Insert(tet_node_index,index);
}
//#####################################################################
// Function Add_Intersection_Id
//#####################################################################
void CUTTING_PARTICLES::
Add_Intersection_Id(const int intersection_index)
{
    if(intersection_to_particle_id.Contains(intersection_index)) return;
    int index=intersection_indices.Append(intersection_index);
    tet_node_indices.Append(0);
    particle_ids_types.Append(INTERSECTION_ID);
    intersection_to_particle_id.Insert(intersection_index,index);
}
//#####################################################################
// Function Add_Tet_Node_And_Intersection_Id
//#####################################################################
void CUTTING_PARTICLES::
Add_Tet_Node_And_Intersection_Id(const int tet_node_index,const int intersection_index)
{
    int index1=tet_node_indices.Append(tet_node_index);
    int index2=intersection_indices.Append(intersection_index);
    assert(index1==index2);
    particle_ids_types.Append(TET_NODE_AND_INTERSECTION_ID);
    tet_node_to_particle_id.Insert(tet_node_index,index1);
    intersection_to_particle_id.Insert(intersection_index,index2);
}
//#####################################################################
// Function Particle_Id_From_Tet_Node
//#####################################################################
int CUTTING_PARTICLES::
Particle_Id_From_Tet_Node(const int tet_node) const
{
    int particle_id=0;
    if(!tet_node_to_particle_id.Get(tet_node,particle_id)) PHYSBAM_FATAL_ERROR();
    return particle_id;
}
//#####################################################################
// Function Particle_Id_From_Intersection
//#####################################################################
int CUTTING_PARTICLES::
Particle_Id_From_Intersection(const int intersection) const
{
    int particle_id=0;
    if(!intersection_to_particle_id.Get(intersection,particle_id)) PHYSBAM_FATAL_ERROR();
    return particle_id;
}
//#####################################################################
// Function Print
//#####################################################################
void CUTTING_PARTICLES::
Print() const
{
    std::stringstream ss;
    for(int i=1;i<=Number();i++) ss<<"Particle id "<<i<<(particle_ids_types(i)==TET_NODE_ID?"(tet node)":(particle_ids_types(i)==INTERSECTION_ID?"(intersection)":"(both)"))
        <<" has tet node id "<<tet_node_indices(i)<<", intersection id "<<intersection_indices(i)<<std::endl;
    LOG::filecout(ss.str());
}
