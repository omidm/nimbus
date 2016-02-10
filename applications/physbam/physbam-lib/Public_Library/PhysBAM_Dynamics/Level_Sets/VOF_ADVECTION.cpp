//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOF_ADVECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Level_Sets/VOF_ADVECTION.h>
#if 0
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#endif
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VOF_ADVECTION<TV>::
VOF_ADVECTION(PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset_input,T_FAST_LEVELSET_ADVECTION& levelset_advection_input,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >* fluids_parameters_input)
    :particle_levelset(particle_levelset_input),levelset_advection(levelset_advection_input),grid(particle_levelset.levelset.grid),object(object_mesh,object_particles),
    old_object(old_object_mesh,old_object_particles),preimage(object),full_cell_size(grid.Cell_Size()),length_scale(grid.Minimum_Edge_Length()),epsilon((T).01*length_scale),
    maximum_refinement_depth(0),minimum_refinement_depth(0),maximum_material_refinement_depth(3),minimum_material_refinement_depth(1),fluids_parameters(fluids_parameters_input),
    laplace(grid,p,particle_levelset.levelset,true,false,true),volume_of_material_initialized(false),geometry_initialized(false)
{
    mpi_grid=particle_levelset.mpi_grid;implicit_object=new PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>(particle_levelset);simplices_in_cells.Resize(grid.Domain_Indices(3));
    postimage_simplices_in_cells.Resize(grid.Domain_Indices(3));volume_of_material.Resize(grid.Domain_Indices(3));geometric_volume.Resize(grid.Domain_Indices(3));p.Resize(grid.Domain_Indices(2));marked_cells.Resize(grid.Domain_Indices(2));
    fixed_cells.Resize(grid.Domain_Indices(3));borders_non_fixed_material_cell.Resize(grid.Domain_Indices(4));borders_fixed_material_cell.Resize(grid.Domain_Indices(4));cell_preimage_material_volume.Resize(grid.Domain_Indices(3));
    volume_fluxes.Resize(grid);removed_negative_particle_volumes_in_cell.Resize(grid.Domain_Indices(1));
    laplace.mpi_grid=mpi_grid;laplace.Set_Up_Second_Order_Cut_Cell_Method();laplace.pcg.Set_Maximum_Iterations(100);laplace.pcg.Show_Results();laplace.pcg.Show_Residuals();
    debugging=false;Set_Runge_Kutta_Order_Particles();Use_Frozen_Velocity();V.Resize(grid);V_ghost.Resize(grid,particle_levelset.number_of_ghost_cells);
    object.particles.array_collection->Preallocate(10000);old_object.particles.array_collection->Preallocate(10000);cell_to_postimage_particle_map.Resize(grid.Domain_Indices(3));
    dirichlet_neighbor.Resize(grid.Domain_Indices(3));

    cells_valid.Resize(grid.Domain_Indices(3));
    // find domain edges which can be traversed (i.e. parallel boundaries or outflow walls)
    bool valid_edges[GRID<TV>::dimension][2]; // TODO: Add outflow walls below
    for(int axis=1;axis<=GRID<TV>::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++) valid_edges[axis-1][axis_side-1]=mpi_grid?mpi_grid->Neighbor(axis,axis_side):false;
    RANGE<TV_INT> domain=grid.Domain_Indices();TV_INT domain_max=domain.Maximum_Corner(),domain_min=domain.Minimum_Corner();
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();
        bool valid=true;
        for(int axis=1;axis<=GRID<TV>::dimension;axis++)
            if((cell(axis)>domain_max(axis) && !valid_edges[axis-1][1]) || (cell(axis)<domain_min(axis) && !valid_edges[axis-1][0])){valid=false;break;}
        cells_valid(cell)=valid;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> VOF_ADVECTION<TV>::
~VOF_ADVECTION()
{
    delete implicit_object;
}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class TV> typename TV::SCALAR VOF_ADVECTION<TV>::
Negative_Material(const TV_INT& cell_index,const bool force_full_refinement)
{
    // TODO: We need to use phi PLS here, not phi LS or RPLS.
    // refinement step
    const GRID<TV>& grid=particle_levelset.levelset.grid;
    static TV_INT phi_indices[GRID<TV>::number_of_nodes_per_cell];grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,phi_indices);

    cell_particle_X.Resize(GRID<TV>::number_of_nodes_per_cell);
    for(int i=1;i<=GRID<TV>::number_of_nodes_per_cell;i++) cell_particle_X(i)=grid.Node(phi_indices[i-1]);
    cell_refinement_simplices.Remove_All();Refined_Object_Initialization_Helper(cell_refinement_simplices);

    int last_node=cell_particle_X.m;
    cell_phis.Resize(last_node);
    // compute phis for extant nodes
    T minimum_phi=FLT_MAX,maximum_phi=-FLT_MAX;

    for(int i=1;i<=last_node;i++){T& phi=cell_phis(i);
        phi=grid_nodal_phis(phi_indices[i-1]);
        if(phi<minimum_phi) minimum_phi=phi;
        if(phi>maximum_phi) maximum_phi=phi;}

    if(minimum_phi*maximum_phi>0 && min(abs(minimum_phi),abs(maximum_phi))>length_scale) return minimum_phi<=0?full_cell_size:0;

    int unrefined_point_count=last_node;
    int last_parent_simplex=0,first_parent_simplex=1;
    for(int depth=0;depth<maximum_material_refinement_depth;depth++){
        first_parent_simplex=last_parent_simplex+1;
        last_parent_simplex=cell_refinement_simplices.m;
        for(int s=last_parent_simplex;s>=first_parent_simplex;s--){
            const T_ELEMENT simplex=cell_refinement_simplices(s);
            if(depth < minimum_material_refinement_depth || Refinement_Condition(cell_particle_X,cell_phis,simplex) || force_full_refinement){
                last_parent_simplex--;
                cell_refinement_simplices.Remove_Index_Lazy(s);
                Refine_Simplex(cell_refinement_simplices,cell_particle_X,simplex);}}
        cell_phis.Resize(cell_particle_X.m);
        // compute phis for extant nodes
        for(int i=last_node+1;i<=cell_phis.m;i++) cell_phis(i)=Phi(cell_particle_X(i));
        last_node=cell_phis.m;}

    if(maximum_material_refinement_depth>0 && cell_phis.m==unrefined_point_count) return minimum_phi<=0?full_cell_size:0;
    // compute material
    T negative_material=0;
    for(int i=1;i<=cell_refinement_simplices.m;i++) negative_material+=T_SIMPLEX::Negative_Material(cell_particle_X,cell_phis,cell_refinement_simplices(i));
    for(int i=1;i<=removed_negative_particle_volumes_in_cell(cell_index).m;i++) negative_material+=removed_negative_particle_volumes_in_cell(cell_index)(i).z;
    // TODO: negative material should be strictly nonnegative - ought to be able to remove the clamp from below, but ought to check this.  Also, precompute the particle contributions.
    //negative_material+=Get_Unmodified_Particle_Volume_In_Cell(particle_levelset.negative_particles,cell_index)+Get_Particle_Volume_In_Cell(particle_levelset.removed_negative_particles,cell_index);
    return clamp(negative_material,(T)0,full_cell_size);
}
//#####################################################################
// Function Precompute_Particle_Volumes_In_Cells
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> void VOF_ADVECTION<TV>::
Precompute_Particle_Volumes_In_Cells(T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE& particle_volumes_in_cell,T_ARRAYS_PARTICLES& particles,int sign)
{
    const GRID<TV>& grid=particle_levelset.levelset.grid;
    T volumes[GRID<TV>::number_of_cells_per_node];TV_INT cells[GRID<TV>::number_of_cells_per_node];
    T one_over_radius_multiplier=-(T)sign;
    // it would be nice to unify removed and escaped particles...for the moment (heinous) assume they're all escaped
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) particle_volumes_in_cell(iterator.Cell_Index()).Remove_All();
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(particles(block)){
            PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(block);grid.Cells_Neighboring_Node(block,cells);
            for(int k=1;k<=cell_particles.array_collection->Size();k++){
                if(one_over_radius_multiplier*particle_levelset.levelset.Phi(cell_particles.X(k))>cell_particles.radius(k)){
                    SPHERE<TV>(cell_particles.X(k),cell_particles.radius(k)).Sector_Volumes(grid.Node(block),volumes,epsilon);
                    for(int i=0;i<GRID<TV>::number_of_cells_per_node;i++)
                        particle_volumes_in_cell(cells[i]).Append(TRIPLE<TV_INT,int,T>(block,k,volumes[i]));}}}} // overadding rather than using if; may be bad
}
//#####################################################################
// Function Get_Outside_Particle_Volume_In_Cell
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> typename TV::SCALAR VOF_ADVECTION<TV>::
Get_Outside_Particle_Volume_In_Cell(T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE& particle_volumes_in_cell,T_ARRAYS_PARTICLES& particles,const TV_INT& cell_index)
{
    T particle_volume=0;
    for(int i=1;i<=particle_volumes_in_cell(cell_index).m;i++) particle_volume+=particle_volumes_in_cell(cell_index)(i).z;
    return particle_volume;
}
//#####################################################################
// Function Modify_Particle_Material_Volume_In_Cell
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> void VOF_ADVECTION<TV>::
Modify_Particle_Material_Volume_In_Cell(T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE& particle_volumes_in_cell,T_ARRAYS_PARTICLES& particles,const TV_INT& cell_index,const T volume_of_material_density)
{
    for(int i=1;i<=particle_volumes_in_cell(cell_index).m;i++){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(particle_volumes_in_cell(cell_index)(i).x);
        ARRAY_VIEW<T>* material_volume=cell_particles.array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);
        (*material_volume)(particle_volumes_in_cell(cell_index)(i).y)+=particle_volumes_in_cell(cell_index)(i).z*volume_of_material_density;}
}
//#####################################################################
// Function Create_Preimage_Particles_From_Old_Postimage_Simplices
//#####################################################################
template<class TV> typename TV::SCALAR VOF_ADVECTION<TV>::
Create_Preimage_Particles_From_Old_Postimage_Simplices(const TV_INT& cell_index,const T cell_volume_of_material)
{
    ARRAY<TV> simplex_centers;
    ARRAY<T_ELEMENT> simplices_in_cell,junk_simplices;
    ARRAY<int>& cell_postimage_simplices=postimage_simplices_in_cells(cell_index);
    for(int i=1;i<=cell_postimage_simplices.m;i++){
        T_ELEMENT& simplex=old_object_mesh.elements(cell_postimage_simplices(i));

        simplex_particles.Remove_All();
        T_ELEMENT local_simplex;
        for(int j=1;j<=T_ELEMENT::dimension;j++){simplex_particles.Append(old_object_particles.X(simplex[j]));local_simplex[j]=j;}

        // dice material simplex
        simplices_in_cell.Remove_All();junk_simplices.Remove_All();simplices_in_cell.Append(local_simplex);
        for(int axis=1;axis<=GRID<TV>::dimension;axis++) for(int node=0;node<=1;node++){
            TV_INT right_cell=cell_index+node*TV_INT::Axis_Vector(axis);
            // cut away outside bits
            T_HYPERPLANE cutting_surface(TV::Axis_Vector(axis),grid.Face(axis,right_cell));
            for(int t=simplices_in_cell.m;t>=1;t--){
                T_ELEMENT simplex=simplices_in_cell(t);simplices_in_cell.Remove_Index_Lazy(t);
                if(node==0) T_SIMPLEX::Cut_With_Hyperplane(simplex_particles,cutting_surface,simplex,junk_simplices,simplices_in_cell);
                else T_SIMPLEX::Cut_With_Hyperplane(simplex_particles,cutting_surface,simplex,simplices_in_cell,junk_simplices);}}
        for(int j=1;j<=simplices_in_cell.m;j++) simplex_centers.Append(ARRAYS_COMPUTATIONS::Average(simplex_particles.Subset(simplices_in_cell(j))));}

    if(!simplex_centers.m){
        if(!cell_to_postimage_particle_map(cell_index)) PHYSBAM_FATAL_ERROR();
        simplex_centers.Append_Elements(old_object_particles.X.Subset(*cell_to_postimage_particle_map(cell_index)));}

    const T volume_per_particle=cell_volume_of_material/simplex_centers.m;
    ARRAY_PLUS_SCALAR<int,IDENTITY_ARRAY<> > added_particles=object.particles.array_collection->Add_Elements(simplex_centers.m);
    for(int i=1;i<=simplex_centers.m;i++){
        object.particles.X(added_particles(i))=simplex_centers(i);
        preimage_particle_material_volumes.Append(volume_per_particle);
        material_particles.Append(true);
        fixed_particle_list.Append(false);}
    return cell_volume_of_material;
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Create_Geometry()
{
    const T total_material_volume=volume_of_material.Sum();
    T accumulated_material_volume=(T)0;
    // Start with an identify and remove?
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(particle_levelset.removed_negative_particles(block)){
            PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particle_levelset.removed_negative_particles(block);
            ARRAY_VIEW<T>* material_volume=cell_particles.array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);
            ARRAYS_COMPUTATIONS::Fill(*material_volume,(T)0);}}
    T phi_meshing_threshold=-10*grid.Maximum_Edge_Length(); // TODO: change this.  for now, no fixed cells

    old_phis.Exchange(phis);
    old_object.particles.array_collection->Initialize(*object.particles.array_collection);
    object.particles.array_collection->Delete_All_Elements();

    old_object.mesh.elements.Exchange(object.mesh.elements);
    old_object.mesh.number_nodes=old_object.particles.array_collection->Size();
    object.mesh.number_nodes=0;
    object.mesh.elements.Remove_All();
    object.mesh.Delete_Auxiliary_Structures();
    fixed_segment_list.Remove_All();fixed_particle_list.Remove_All();
    material_particles.Remove_All();

    Precompute_Cell_Particle_Influence();
    Precompute_Particle_Volumes_In_Cells(removed_negative_particle_volumes_in_cell,particle_levelset.removed_negative_particles,-1);

    TV_INT nodes_in_cell[GRID<TV>::number_of_nodes_per_cell];
    fixed_cells.Fill(false);borders_non_fixed_material_cell.Fill(false);borders_fixed_material_cell.Fill(false);
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(cells_valid(cell_index)){
            simplices_in_cells(cell_index).Remove_All();
            if(particle_levelset.levelset.phi(cell_index)<phi_meshing_threshold){
                fixed_cells(cell_index)=true;
                grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,nodes_in_cell);
                for(int i=0;i<GRID<TV>::number_of_nodes_per_cell;i++) borders_fixed_material_cell(nodes_in_cell[i])=true;}
            else if(volume_of_material(cell_index)>0){
                grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,nodes_in_cell);
                for(int i=0;i<GRID<TV>::number_of_nodes_per_cell;i++) borders_non_fixed_material_cell(nodes_in_cell[i])=true;}}}

    // start by adding common nodes to object, then red simplices.  Then use IT to refine.
    nodes_to_particles_map.Resize(grid.Node_Indices(3));
    for(NODE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
        if(borders_non_fixed_material_cell(node_index)){
            int particle_index=object.particles.array_collection->Add_Element();
            object.particles.X(particle_index)=iterator.Location();
            fixed_particle_list.Append(borders_fixed_material_cell(node_index));
            nodes_to_particles_map(node_index)=particle_index;}}

    for(int level=1;level<=level_simplex_cells.m;level++) level_simplex_cells(level).Remove_All();

    level_simplex_cells.Resize(1);
    // create initial red refinement
    cell_preimage_material_volume.Fill((T)0);
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(fixed_cells(cell_index)) cell_preimage_material_volume(cell_index)=volume_of_material(cell_index);
        else if (volume_of_material(cell_index)>0 && cells_valid(cell_index)) Create_Initial_Meshing_For_Cell(grid,object,cell_index,nodes_to_particles_map,level_simplex_cells(1));}

    accumulated_material_volume+=cell_preimage_material_volume.Sum();

    object.mesh.Set_Number_Nodes(object.particles.array_collection->Size());
    object.mesh.Initialize_Element_Edges();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after initial meshing",0,0);
    preimage.Initialize();
    fixed_segment_list.Resize(preimage.segment_mesh.elements.m);
    for(int i=1;i<=preimage.segment_mesh.elements.m;i++){
        const VECTOR<int,2>& segment=preimage.segment_mesh.elements(i);
        fixed_segment_list(i)=fixed_particle_list(segment[1]) && fixed_particle_list(segment[2]);}

    // red green refine
    int last_particle_index=0,last_segment_index=fixed_segment_list.m;
    phis.Resize(object.particles.array_collection->Size());
    for(int i=last_particle_index+1;i<=object.particles.array_collection->Size();i++) phis(i)=Phi(object.particles.X(i));
    last_particle_index=object.particles.array_collection->Size();
    for(int depth=0;depth<minimum_refinement_depth;depth++){
        simplex_list.Remove_All();
        if(depth<minimum_refinement_depth) for(int i=1;i<=preimage.meshes(depth+1)->elements.m;i++) simplex_list.Append((*preimage.leaf_number(depth+1))(i));
        else for(int i=1;i<=preimage.meshes(depth+1)->elements.m;i++){
            const T_ELEMENT& simplex=preimage.meshes(depth+1)->elements(i);
            if(Refinement_Condition(phis.Subset(simplex),object.particles.X.Subset(simplex))) simplex_list.Append((*preimage.leaf_number(depth+1))(i));}
        preimage.Refine_Simplex_List(simplex_list);

        preimage.Initialize_Segment_Index_From_Midpoint_Index(); // TODO: may not want to keep calling this

        phis.Resize(object.particles.array_collection->Size());fixed_particle_list.Resize(object.particles.array_collection->Size());
        for(int i=last_particle_index+1;i<=object.particles.array_collection->Size();i++){
            phis(i)=Phi(object.particles.X(i));
            fixed_particle_list(i)=fixed_segment_list((*preimage.segment_index_from_midpoint_index)(i));}
        last_particle_index=object.particles.array_collection->Size();
        
        // update fixed segment list
        fixed_segment_list.Resize(preimage.segment_mesh.elements.m);
        for(int s=last_segment_index+1;s<=preimage.segment_mesh.elements.m;s++)
            fixed_segment_list(s)=fixed_particle_list(preimage.segment_mesh.elements(s)(1)) && fixed_particle_list(preimage.segment_mesh.elements(s)(2));
        last_segment_index=preimage.segment_mesh.elements.m;

        // hah hah hah, sucks to be you if you have to debug this!
        level_simplex_cells.Resize(depth+2);
        level_simplex_cells(depth+2).Resize(preimage.meshes(depth+2)->elements.m);
        for(int i=1;i<=preimage.meshes(depth+1)->elements.m;i++){
            int child_count=0,child_simplex;
            while(child_count<1<<TV::dimension && (child_simplex=(*preimage.children(depth+1))(i)(++child_count)))
                level_simplex_cells(depth+2)(child_simplex)=level_simplex_cells(depth+1)(i);}}

    Cut_Simplicial_Object_With_Zero_Isocontour();

    // TODO: do we need this?
    object.mesh.Refresh_Auxiliary_Structures();

    // find cells with material volume and no simplices, and move their material into cells with simplices (in the minus-phi-normal direction)
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(volume_of_material(cell_index) && !simplices_in_cells(cell_index).m){
            TV normal=T_LEVELSET::Normal_At_Node(grid,particle_levelset.levelset.phi,cell_index);
            TV_INT target_cell_index=grid.Clamp_To_Cell(iterator.Location()-grid.dX*normal); // TODO: this will need to change for parallel
            assert(simplices_in_cells(target_cell_index).m);
            volume_of_material(target_cell_index)+=volume_of_material(cell_index);
            volume_of_material(cell_index)=(T)0;}}

    debugging=true;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before volume assignment",0,0);
    debugging=false;

    ARRAY<T> leaf_simplex_material_volumes;
    ARRAY<T> leaf_simplex_geometric_volumes;
    leaf_simplex_material_volumes.Resize(object.mesh.elements.m);
    leaf_simplex_geometric_volumes.Resize(object.mesh.elements.m);

    last_simplex_particle=object_particles.array_collection->Size(); // this is to determine which particles are anonymous preimage particles
    preimage_particle_material_volumes.Remove_All();
    // add material volume to simplices and particles
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(!fixed_cells(cell_index) && volume_of_material(cell_index) && cells_valid(cell_index)){
            // compute simplex geometric volumes in cell
            T cell_material_volume=(T)0;
            T geometric_volume_in_cell=(T)0;
            ARRAY<int>& simplices_in_cell=simplices_in_cells(cell_index);
            for(int i=1;i<=simplices_in_cell.m;i++){
                int simplex_index=simplices_in_cell(i);
                leaf_simplex_geometric_volumes(simplex_index)=T_SIMPLEX::Signed_Size(object.particles.X.Subset(object.mesh.elements(simplex_index)));
                geometric_volume_in_cell+=leaf_simplex_geometric_volumes(simplex_index);}

            T negative_particle_volume=Get_Outside_Particle_Volume_In_Cell(removed_negative_particle_volumes_in_cell,particle_levelset.removed_negative_particles,cell_index);
            if(simplices_in_cell.m==0 && negative_particle_volume==0){
                cell_material_volume+=Create_Preimage_Particles_From_Old_Postimage_Simplices(cell_index,volume_of_material(cell_index));
                accumulated_material_volume+=cell_material_volume;
                continue;}

            for(int i=fixed_particle_list.m+1;i<=object.particles.array_collection->Size();i++) fixed_particle_list.Append(false);
            
            // compute simplex negative material and label used nodes
            int number_of_simplices=simplices_in_cell.m;
            simplex_preimage_volume.Resize(number_of_simplices);lower_dimensional_preimage_volume.Resize(number_of_simplices);
            lower_dimensional_preimage.Resize(number_of_simplices);ARRAYS_COMPUTATIONS::Fill(lower_dimensional_preimage,false);
            for(int i=material_particles.m+1;i<=object.particles.array_collection->Size();i++) material_particles.Append(false);
            for(int i=1;i<=simplices_in_cell.m;i++){
                const T_ELEMENT& simplex=object.mesh.elements(simplices_in_cell(i));
                for(int j=1;j<=TV::dimension+1;j++) material_particles(simplex[j])=true;
                simplex_preimage_volume(i)=T_SIMPLEX::Signed_Size(object.particles.X.Subset(simplex));
                if(simplex_preimage_volume(i)<=0){lower_dimensional_preimage(i)=true;lower_dimensional_preimage_volume(i)=T_SIMPLEX::Half_Boundary_Measure(object.particles.X.Subset(simplex));}}

            simplex_preimage_material_volume.Resize(object.mesh.elements.m);
            const T total_preimage_volume=ARRAYS_COMPUTATIONS::Sum(simplex_preimage_volume)+negative_particle_volume;
            if(total_preimage_volume<=0){
                const T total_lower_dimensional_preimage_volume=ARRAYS_COMPUTATIONS::Sum(lower_dimensional_preimage_volume);
                // TODO: need to modify particles here too, just in case we have any (although that's unlikely)
                if(!total_lower_dimensional_preimage_volume){ // divide it up evenly...not that it's likely to matter in this case
                    T average_simplex_volume=number_of_simplices?volume_of_material(cell_index)/number_of_simplices:0;
                    for(int i=1;i<=number_of_simplices;i++) simplex_preimage_material_volume(simplices_in_cell(i))=average_simplex_volume;}
                else{
                    T volume_of_material_density=volume_of_material(cell_index)/total_lower_dimensional_preimage_volume;
                    for(int i=1;i<=number_of_simplices;i++) simplex_preimage_material_volume(simplices_in_cell(i))=lower_dimensional_preimage_volume(i)*volume_of_material_density;}
                for(int i=1;i<=number_of_simplices;i++) cell_material_volume+=simplex_preimage_material_volume(simplices_in_cell(i));}
            else{
                T volume_of_material_density=volume_of_material(cell_index)/total_preimage_volume;
                cell_material_volume+=volume_of_material_density*negative_particle_volume;
                Modify_Particle_Material_Volume_In_Cell(removed_negative_particle_volumes_in_cell,particle_levelset.removed_negative_particles,cell_index,volume_of_material_density);
                for(int i=1;i<=number_of_simplices;i++){simplex_preimage_material_volume(simplices_in_cell(i))=simplex_preimage_volume(i)*volume_of_material_density;
                    cell_material_volume+=simplex_preimage_material_volume(simplices_in_cell(i));}}
            if(abs(cell_material_volume-volume_of_material(cell_index))/volume_of_material(cell_index)>.1){
                std::stringstream ss;
                ss<<"Cell "<<cell_index<<": was "<<volume_of_material(cell_index)<<", now "<<cell_material_volume<<std::endl;}
            accumulated_material_volume+=cell_material_volume;}}

    std::stringstream ss;
    ss<<"Before: "<<total_material_volume<<" after: "<<accumulated_material_volume<<std::endl;
    LOG::filecout(ss.str());
    object.mesh.Set_Number_Nodes(object.particles.array_collection->Size());
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after volume assignment",0,0);
    object.mesh.Refresh_Auxiliary_Structures();
    Initialize_Face_Neighbors(object);
    
    // adjust removed particle radii based on assigned material volume.
    const T one_over_unit_sphere_size=(T)1/(T)unit_sphere_size<TV::dimension>::value;
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(particle_levelset.removed_negative_particles(block)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& cell_particles=*particle_levelset.removed_negative_particles(block);
            ARRAY_VIEW<T>* material_volume=cell_particles.array_collection->template Get_Array<T>(ATTRIBUTE_ID_MATERIAL_VOLUME);
            for(int k=1;k<=cell_particles.array_collection->Size();k++)
                cell_particles.radius(k)=max(particle_levelset.minimum_particle_radius,(T)pow<1,TV::dimension>((*material_volume)(k)*one_over_unit_sphere_size));}}

    preimage.Initialize(); // recreate the preimage after all the surgery has been done
}
//#####################################################################
// Function Cut_Simplicial_Object_With_Zero_Isocontour
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Cut_Simplicial_Object_With_Zero_Isocontour()
{
    int number_of_particles=object.particles.array_collection->Size();
    object.mesh.Initialize_Segment_Mesh();Initialize_Segment_Neighbors(object);

    // create list of cut segments and list of cut simplices
    ARRAY<int> inefficient_segment_list;
    ARRAY<VECTOR<int,2> > cut_segments;ARRAY<int> cut_simplices;
    for(int s=1;s<=object.mesh.segment_mesh->elements.m;s++){
        int i,j;object.mesh.segment_mesh->elements(s).Get(i,j);
        if(phis(i)*phis(j)>=0) continue;
        inefficient_segment_list.Append(preimage.segment_mesh.Simplex(VECTOR<int,2>(i,j)));
        cut_segments.Append(VECTOR<int,2>(i,j));cut_simplices.Append_Elements(Segment_Neighbors(object,s));}

    Initialize_Refinement_Target_Segments(inefficient_segment_list,preimage);
    preimage.Refine_Simplex_List(cut_simplices);
    preimage.Initialize_Segment_Index_From_Midpoint_Index(); // TODO: may not want to keep calling this
    Clean_Refinement_Target_Segments(preimage);

    material_particles.Resize(object.particles.array_collection->Size(),false);
    for(int i=1;i<=phis.m;i++) material_particles(i)=(phis(i)<=0);
    phis.Resize(object.particles.array_collection->Size(),false);fixed_particle_list.Resize(object.particles.array_collection->Size(),false);
    for(int p=number_of_particles+1;p<=object.particles.array_collection->Size();p++){
        VECTOR<int,2> segment=preimage.segment_mesh.elements((*preimage.segment_index_from_midpoint_index)(p));
        int i=segment(1),j=segment(2);
        if(phis(i)*phis(j)<0){phis(p)=0;material_particles(p)=true;}
        else{phis(p)=(T).5*(phis(i)+phis(j));material_particles(p)=(material_particles(i) && material_particles(j));}
        fixed_particle_list(p)=(fixed_particle_list(i) && fixed_particle_list(j));}

    // update level_simplex_cells and negative material simplices
    // TODO: do we need level_simplex_cells ?
    ARRAY<int> simplices_to_remove;
    level_simplex_cells.Resize(preimage.meshes.m);
    for(int level=1;level<=level_simplex_cells.m;level++){
        level_simplex_cells(level).Resize(preimage.meshes(level)->elements.m);
        for(int t=1;t<=level_simplex_cells(level).m;t++)
            if(preimage.Leaf(level,t)){
                if(level>1) level_simplex_cells(level)(t)=level_simplex_cells(level-1)((*preimage.parent(level))(t));
                int leaf=(*preimage.leaf_number(level))(t);
                if(material_particles.Subset(object.mesh.elements(leaf)).Number_True()!=T_ELEMENT::dimension) simplices_to_remove.Append(leaf);}}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("before removal",0,0);
    ARRAY<HASHTABLE<int,int> > level_simplex_maps;
    preimage.Remove_Simplex_List(simplices_to_remove,&level_simplex_maps);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after removal before adjustment",0,0);
    // adjust the positions of the relevant midpoints to lie on the interface
    for(int s=1;s<=cut_segments.m;s++){
        int i,j;cut_segments(s).Get(i,j);
        int midpoint=preimage.segment_midpoints(preimage.segment_mesh.Segment(i,j));
        object.particles.X(midpoint)=Iterative_Find_Interface(object.particles.X(i),object.particles.X(j));}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after removal after adjustment",0,0);

    for(int level=1;level<=level_simplex_cells.m;level++){
        // for each level, apply the map to the level_simplex_cells
        ARRAY<TV_INT>& simplex_cells=level_simplex_cells(level);
        HASHTABLE<int,int>& simplex_map=level_simplex_maps(level);
        for(int i=1;i<=simplex_cells.m;i++){
            TV_INT& cell_index=simplex_cells(i);
            int index=i;
            simplex_map.Get(i,index);
            if(index && preimage.Leaf(level,index)) simplices_in_cells(cell_index).Append((*preimage.leaf_number(level))(index));}}

    // TODO: something here maybe to move stuff up to top level if it has no parents.  Ideally, call this after we don't care about simplices in cells no mo'
}
//#####################################################################
// Function Iterative_Find_Interface
//#####################################################################
template<class TV> TV VOF_ADVECTION<TV>::
Iterative_Find_Interface(TV left,TV right,const int iterations) const
{
    T phi_left=Phi(left),phi_right=Phi(right);
    TV theta=LINEAR_INTERPOLATION<T,TV>::Linear(left,right,LEVELSET_UTILITIES<T>::Theta(phi_left,phi_right));
    int phi_left_sign=(phi_left<=0?-1:1),phi_right_sign=(phi_right<=0?-1:1);
    for(int i=1;i<=iterations;i++){
        T phi=Phi(theta);
        int phi_sign=(phi<=0?-1:1);
        if(phi_left_sign*phi_sign<0){
            right=theta;phi_right=phi;phi_right_sign=phi_sign;
            theta=LINEAR_INTERPOLATION<T,TV>::Linear(left,theta,LEVELSET_UTILITIES<T>::Theta(phi_left,phi));}
        else if(phi_right_sign*phi_sign<0){
            left=theta;phi_left=phi;phi_left_sign=phi_sign;
            theta=LINEAR_INTERPOLATION<T,TV>::Linear(theta,right,LEVELSET_UTILITIES<T>::Theta(phi,phi_right));}
        else break;}
    return theta;
}
//#####################################################################
// Function Advect_Material_Preimages
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Advect_Material_Preimages(const T_FACE_LOOKUP& face_velocities,const T dt,const T start_time,const bool use_analytic_velocities)
{
    old_particle_X=object.particles.X;
    // being lazy
    ARRAY_VIEW<TV> particle_X(object.particles.X);
    const int particle_number=object.particles.array_collection->Size();
    // TODO: store copy of old particle positions
    if(runge_kutta_order_particles == 2 && use_frozen_velocity){
        // move material nodes, second-order
        for(int i=1;i<=particle_number;i++) if(material_particles(i) && !fixed_particle_list(i)){ // TODO: avoid conditionals here by keeping a list of material particles rather than flags
            TV& node=particle_X(i);
            TV offset=dt*Velocity(face_velocities,node,start_time,use_analytic_velocities);
            TV second_offset=dt*Velocity(face_velocities,node+offset,start_time,use_analytic_velocities);
            node+=(T).5*(offset+second_offset); // TODO: clamp to domain walls (not parallel/outflow walls)
            Adjust_Node_For_Domain_Boundaries(node);}
    }
    else if(runge_kutta_order_particles == 2 || runge_kutta_order_particles == 3){
        T time=start_time;
        // TODO: this is inefficient.  Need to only keep the particles we care about.  Do this by building a list and changing the indices on the fly?
        RUNGEKUTTA<ARRAY_VIEW<TV> >* rungekutta=RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(particle_X,runge_kutta_order_particles,dt,time);
        for(int k=1;k<=runge_kutta_order_particles;k++){
            if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
            // TODO: ghost velocities should not be necessary unless we have periodic boundary conditions or we're on the edge or we're doing parallel...I guess they might be necessary
            particle_levelset.levelset.boundary->Fill_Ghost_Cells_Face(grid,V,V_ghost,time,particle_levelset.number_of_ghost_cells);
            T_FACE_LOOKUP V_lookup(V_ghost);
            
            for(int i=1;i<=particle_number;i++) if(material_particles(i) && !fixed_particle_list(i)){ // TODO: avoid conditionals here by keeping a list of material particles rather than flags
                TV& node=particle_X(i);
                TV velocity=Velocity(V_lookup,node,start_time,use_analytic_velocities);
                node+=dt*velocity; // TODO: clamp to domain walls (not parallel/outflow walls)
                Adjust_Node_For_Domain_Boundaries(node);}
            time=rungekutta->Main();}
        delete rungekutta;}

    // what the heck do we want to do for higher order stuff here?  Time integrate on cell faces, or what?
    // transfer material from one simplex to another first
    // first do loop over grid faces.  Two interior cells do the transfer; a cell adjacent to simplices lets the simplices handle it based on area (in fact, should be automatic;
    // only special case is to single count the face, rather than double

    // ALE section
    TV dt_face_size=dt*grid.Face_Sizes();
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        if(fixed_cells(iterator.First_Cell_Index()) && fixed_cells(iterator.Second_Cell_Index())){
            T material_flux=face_velocities.V_face.Component(iterator.Axis())(iterator.Face_Index())*dt_face_size(iterator.Axis());
            cell_preimage_material_volume(iterator.First_Cell_Index())-=material_flux;
            cell_preimage_material_volume(iterator.Second_Cell_Index())+=material_flux;}
    volume_of_material=cell_preimage_material_volume;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after apply face fluxes",0,0);

    if(!use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,start_time+dt);
    T_FACE_LOOKUP face_velocities_initial(face_velocities.V_face),face_velocities_final(use_frozen_velocity?face_velocities.V_face:V);
    
    const T maximum_refinement_fraction=(T)1/(1<<(maximum_refinement_depth+2));
    for(int f=1;f<=Faces(object);f++){
        const VECTOR<int,2>& face_neighbors=Face_Neighbors(object,f);
        const T_FACE_ELEMENT& face=Face(object,f);
        int fixed_nodes=0;
        for(int j=1;j<=T_FACE_ELEMENT::dimension;j++) if(fixed_particle_list(face[j])) fixed_nodes++;
        TV face_normal;
        TV_INT opposing_cell;
        TV old_center=ARRAYS_COMPUTATIONS::Average(old_particle_X.Subset(face)),new_center=ARRAYS_COMPUTATIONS::Average(particle_X.Subset(face));
        TV initial_levelset_center_velocity=Velocity(face_velocities_initial,old_center,start_time,use_analytic_velocities);
        TV final_levelset_center_velocity=Velocity(face_velocities_final,new_center,start_time+dt,use_analytic_velocities);
        if(fixed_nodes==T_FACE_ELEMENT::dimension && Find_Fixed_Cell_On_Point(new_center,opposing_cell,face_normal,maximum_refinement_fraction)){
            T material_flux=TV::Dot_Product(face_normal,(T).5*(initial_levelset_center_velocity+final_levelset_center_velocity)*dt)*T_SIMPLEX_FACE::Size(old_particle_X.Subset(face));
            cell_preimage_material_volume(opposing_cell)-=material_flux;
            simplex_preimage_material_volume(face_neighbors[1]?face_neighbors[1]:face_neighbors[2])+=material_flux;}
        else if(face_neighbors[1] && face_neighbors[2]){
            TV center_position_change=new_center-old_center;
            T old_normal_velocity_product=TV::Dot_Product(T_SIMPLEX_FACE::Normal(old_particle_X.Subset(face)),initial_levelset_center_velocity*dt-center_position_change)*T_SIMPLEX_FACE::Size(old_particle_X.Subset(face));
            T new_normal_velocity_product=TV::Dot_Product(T_SIMPLEX_FACE::Normal(particle_X.Subset(face)),final_levelset_center_velocity*dt-center_position_change)*T_SIMPLEX_FACE::Size(object.particles.X.Subset(face));
            T material_flux=(T).5*(old_normal_velocity_product+new_normal_velocity_product);
            simplex_preimage_material_volume(face_neighbors[1])-=material_flux;
            simplex_preimage_material_volume(face_neighbors[2])+=material_flux;}}
    volume_of_material=cell_preimage_material_volume;
}
//#####################################################################
// Function Rasterize_Material_Postimages
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Rasterize_Material_Postimages()
{
    const GRID<TV>& grid=particle_levelset.levelset.grid;const T cell_size=grid.Cell_Size();
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(fixed_cells(cell_index)) geometric_volume(cell_index)=cell_size;
        else geometric_volume(cell_index)=0;}
    volume_of_material=cell_preimage_material_volume;
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()) postimage_simplices_in_cells(iterator.Cell_Index()).Remove_All();

    PHYSBAM_DEBUG_WRITE_SUBSTEP("rasterize",0,0);
    for(int i=1;i<=object.mesh.elements.m;i++){
        // prepare to dice material simplices with grid
        T_ELEMENT simplex=object.mesh.elements(i);
        simplex_particles.Remove_All();
        T_ELEMENT local_simplex;
        for(int j=1;j<=T_ELEMENT::dimension;j++){simplex_particles.Append(object.particles.X(simplex[j]));local_simplex[j]=j;}
        TV minimum_node=simplex_particles(1),maximum_node=simplex_particles(1);
        for(int j=2;j<=T_ELEMENT::dimension;j++){minimum_node=TV::Componentwise_Min(minimum_node,simplex_particles(j));maximum_node=TV::Componentwise_Max(maximum_node,simplex_particles(j));}
        RANGE<TV_INT> grid_cells(grid.Clamp_To_Cell(minimum_node),grid.Clamp_To_Cell(maximum_node)); // TOOD: for parallel, this will need to allow border regions.  Also for outflow walls.
        RANGE<TV_INT> grid_nodes=grid_cells+RANGE<TV_INT>::Unit_Box();
        const TV_INT &low_node=grid_nodes.Minimum_Corner(),&high_node=grid_nodes.Maximum_Corner();
        TV_INT cell_base=low_node-TV_INT::All_Ones_Vector();
        GRID<TV> box_grid(high_node-low_node,RANGE<TV>::Centered_Box(),true);

        for(int j=1;j<=simplex_lists.array.Size();j++) simplex_lists.array(j).Remove_All(); // TODO: this will likely be more expensive than necessary
        simplex_lists.Resize_In_Place(grid_cells);
        simplex_lists(low_node).Remove_All();
        simplex_lists(low_node).Append(local_simplex);

        bool lower_dimensional_postimage=false;
        if(Signed_Size(local_simplex,simplex_particles)==0) lower_dimensional_postimage=true;

        // dice material simplices
        for(int axis=1;axis<=GRID<TV>::dimension;axis++){
            for(int node=1;node<high_node(axis)-low_node(axis);node++){ // only need to cut in-between nodes
                for(FACE_ITERATOR iterator(box_grid,0,GRID<TV>::INTERIOR_REGION,0,axis);iterator.Valid();iterator.Next()){ // TODO: is this set of loops doing more work than necessary?
                    // for all tets in left, cut them and place pieces into right/left, then move to next cell
                    TV_INT left_cell=cell_base+iterator.First_Cell_Index(),right_cell=cell_base+iterator.Second_Cell_Index();
                    T_HYPERPLANE cutting_surface(TV::Axis_Vector(axis),grid.Face(axis,right_cell));
                    ARRAY<T_ELEMENT > &left_list=simplex_lists(left_cell),&right_list=simplex_lists(right_cell);
                    for(int t=left_list.m;t>=1;t--){
                        T_ELEMENT simplex=left_list(t);left_list.Remove_Index_Lazy(t);
                        T_SIMPLEX::Cut_With_Hyperplane(simplex_particles,cutting_surface,simplex,left_list,right_list);}}}}

        // compute cut material simplex volumes FIRST, because for sliver simplices they are likely NOT to sum to the volume of their parent
        cut_material_simplex_postimage_volumes.Resize_In_Place(grid_cells);
        T simplex_postimage_volume=0;
        for(CELL_ITERATOR iterator(grid,grid_cells);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            ARRAY<T_ELEMENT>& simplices=simplex_lists(cell);
            ARRAY<T>& cut_material_simplex_postimage_volumes_cell=cut_material_simplex_postimage_volumes(cell);
            cut_material_simplex_postimage_volumes_cell.Resize(simplices.m);
            for(int t=1;t<=simplices.m;t++){
                if(!lower_dimensional_postimage) cut_material_simplex_postimage_volumes_cell(t)=Signed_Size(simplices(t),simplex_particles);
                else cut_material_simplex_postimage_volumes_cell(t)=Half_Boundary_Measure(simplices(t),simplex_particles);
                simplex_postimage_volume+=cut_material_simplex_postimage_volumes_cell(t);}}
    
        // if stuff gets wedged into a corner, for instance:
        if(simplex_postimage_volume==0) for(CELL_ITERATOR iterator(grid,grid_cells);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            ARRAY<T_ELEMENT>& simplices=simplex_lists(cell);
            ARRAY<T>& cut_material_simplex_postimage_volumes_cell=cut_material_simplex_postimage_volumes(cell);
            for(int t=1;t<=simplices.m;t++) cut_material_simplex_postimage_volumes_cell(t)=(T)1;
            simplex_postimage_volume+=(T)simplices.m;}
        
        T volume_of_material_density=simplex_preimage_material_volume(i)/simplex_postimage_volume;

        // now, get volumes and masses
        for(CELL_ITERATOR iterator(grid,grid_cells);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            ARRAY<T_ELEMENT >& simplices=simplex_lists(cell);
            if(simplices.m){
                postimage_simplices_in_cells(cell).Append(i);
                ARRAY<T>& cut_material_simplex_postimage_volumes_cell=cut_material_simplex_postimage_volumes(cell);
                RANGE<TV> cell_domain=grid.Cell_Domain(cell);
                for(int t=1;t<=simplices.m;t++){
                    // TODO: add something here to flip/zero sign of simplex volumes if they don't match parent sign (above)
                    volume_of_material(cell)+=cut_material_simplex_postimage_volumes_cell(t)*volume_of_material_density;
                    geometric_volume(cell)+=cut_material_simplex_postimage_volumes_cell(t);}}}}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("rasterize before particles",0,0);

    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){delete cell_to_postimage_particle_map(iterator.Cell_Index());
        cell_to_postimage_particle_map(iterator.Cell_Index())=0;}
#if 0
    for(int i=last_simplex_particle+1;i<=object_particles.array_collection->Size();i++){
        const TV_INT cell=grid.Clamp_To_Cell(object_particles.X(i),3);
        const T material_volume=preimage_particle_material_volumes(i-last_simplex_particle);
        volume_of_material(cell)+=material_volume;
        geometric_volume(cell)+=material_volume;
        if(!cell_to_postimage_particle_map(cell)) cell_to_postimage_particle_map(cell)=new ARRAY<int>();
        (*cell_to_postimage_particle_map(cell)).Append(i);
    }
#endif
    PHYSBAM_DEBUG_WRITE_SUBSTEP("rasterize after particles",0,0);
}
//#####################################################################
// Function Adjust_Node_For_Domain_Boundaries
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Adjust_Node_For_Domain_Boundaries(TV& node)
{
    assert(fluids_parameters);
    TV min_corner=fluids_parameters->grid->domain.Minimum_Corner(),max_corner=fluids_parameters->grid->domain.Maximum_Corner();
    for(int axis=1;axis<=GRID<TV>::dimension;axis++){
        if(fluids_parameters->domain_walls(axis)(1) && node[axis] < min_corner[axis]) node[axis]=min_corner[axis];
        if(fluids_parameters->domain_walls(axis)(2) && node[axis] > max_corner[axis]) node[axis]=max_corner[axis];}
}
//#####################################################################
// Function Make_Approximately_Incompressible
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Make_Approximately_Incompressible(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(!fluids_parameters->analytic_test){
        Rasterize_Material_Postimages();
        const MPI_UNIFORM_GRID<GRID<TV> >* mpi_grid=particle_levelset.mpi_grid;
        T_FACE_ARRAYS_BOOL& psi_N=laplace.psi_N;T_ARRAYS_BOOL& psi_D=laplace.psi_D;
        const T cell_size=grid.Cell_Size(),one_over_cell_size=1/cell_size;

        LOG::Time("entering iteration of Make_Incompressible");
        laplace.f.Fill((T)0);

        if(mpi_grid) mpi_grid->Exchange_Boundary_Cell_Data(volume_of_material,3,false); // TODO: do we need this here?

        psi_N.Fill(false);psi_D.Fill(false);
        fluids_parameters->Set_Domain_Boundary_Conditions(laplace,volume_fluxes,time);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(particle_levelset.levelset.phi(cell_index)>0){psi_D(cell_index)=true;p(cell_index)=0;}
            else psi_D(cell_index)=false;}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before second volume step",0,0);

        dirichlet_neighbor.Fill(false);
        // Set up boundary conditions based on whether there are any
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(psi_D(cell_index)) for(int axis=1;axis<=TV::dimension;axis++){
                dirichlet_neighbor(cell_index+TV_INT::Axis_Vector(axis))=true;
                dirichlet_neighbor(cell_index-TV_INT::Axis_Vector(axis))=true;}}
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(!dirichlet_neighbor(cell_index) || volume_of_material(cell_index)>cell_size)
                laplace.f(cell_index)=(volume_of_material(cell_index)-cell_size)*one_over_cell_size;}

        laplace.Solve(time,false);

        const TV dx=grid.dX,one_over_dx=grid.one_over_dX;
        volume_fluxes.Fill((T)0);

        // COMPENSATE FOR COLLISIONS
        if(!laplace.second_order_cut_cell_method){
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(!psi_N(axis,face) && !(psi_D(first_cell) && psi_D(second_cell))){
                    volume_fluxes(axis,face)=(p(second_cell)-p(first_cell))*one_over_dx[axis];}}}
        else{
            T_ARRAYS_SCALAR& phi=particle_levelset.levelset.phi;
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
                if(!psi_N(axis,face_index) && !(psi_D(first_cell) && psi_D(second_cell))){
                    if(psi_D(first_cell) && !psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(second_cell),phi(first_cell)))
                        volume_fluxes(axis,face_index)=(p(second_cell)-laplace.u_interface(axis,face_index))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(second_cell),phi(first_cell)),laplace.second_order_cut_cell_threshold)*dx[axis]);
                    else if(!psi_D(first_cell) && psi_D(second_cell) && LEVELSET_UTILITIES<T>::Interface(phi(first_cell),phi(second_cell)))
                        volume_fluxes(axis,face_index)=(laplace.u_interface(axis,face_index)-p(first_cell))
                            /(max(LEVELSET_UTILITIES<T>::Theta(phi(first_cell),phi(second_cell)),laplace.second_order_cut_cell_threshold)*dx[axis]);
                    else volume_fluxes(axis,face_index)=(p(second_cell)-p(first_cell))*one_over_dx[axis];}}}

        bool store_use_frozen_velocity=use_frozen_velocity;
        int store_runge_kutta_order_particles=runge_kutta_order_particles;
        Set_Runge_Kutta_Order_Particles(2);
        Use_Frozen_Velocity(true);
        T_FACE_ARRAYS_SCALAR volume_fluxes_ghost;volume_fluxes_ghost.Resize(grid,particle_levelset.number_of_ghost_cells,false);
        particle_levelset.levelset.boundary->Fill_Ghost_Cells_Face(grid,volume_fluxes,volume_fluxes_ghost,time,particle_levelset.number_of_ghost_cells);

        T_FACE_LOOKUP volume_flux_lookup(volume_fluxes_ghost);

        Advect_Material_Preimages(volume_flux_lookup,(T)1,(T)0,false);
        Set_Runge_Kutta_Order_Particles(store_runge_kutta_order_particles);
        Use_Frozen_Velocity(store_use_frozen_velocity);

        Second_Order_Runge_Kutta_Step_Particles(volume_flux_lookup,particle_levelset.negative_particles,particle_levelset,PARTICLE_LEVELSET_NEGATIVE,(T)1,time);
        Second_Order_Runge_Kutta_Step_Particles(volume_flux_lookup,particle_levelset.positive_particles,particle_levelset,PARTICLE_LEVELSET_POSITIVE,(T)1,time);
        // TODO: need to correct the removed particle velocities.  However, it's unlikely that removed particles will be affected by this.
        // NOTE: false.  If removed particles overfill a cell, they can certainly be affected.  Also they can be near cells which are overfilled.
        Second_Order_Runge_Kutta_Step_Particles(volume_flux_lookup,particle_levelset.removed_negative_particles,particle_levelset,PARTICLE_LEVELSET_NEGATIVE,(T)1,time);
        Second_Order_Runge_Kutta_Step_Particles(volume_flux_lookup,particle_levelset.removed_positive_particles,particle_levelset,PARTICLE_LEVELSET_POSITIVE,(T)1,time);
        levelset_advection.Euler_Step(volume_fluxes_ghost,(T)1,time,particle_levelset.number_of_ghost_cells);
        volume_fluxes*=1/dt;
        face_velocities+=volume_fluxes;

        PHYSBAM_DEBUG_WRITE_SUBSTEP("after second volume step",0,0);

        particle_levelset.Exchange_Overlap_Particles();
        if(particle_levelset.mpi_grid){
            MPI_UNIFORM_GRID<GRID<TV> >& mpi_grid=*particle_levelset.mpi_grid;
            mpi_grid.Exchange_Boundary_Cell_Data(volume_of_material,particle_levelset.number_of_ghost_cells,true);}}
}
//#####################################################################
// Function Refine_Or_Coarsen_Geometry
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Refine_Or_Coarsen_Geometry()
{
    ARRAY<int> triangle_list;
    ARRAY<ARRAY<T> > level_material_volume(preimage.meshes.m);
    SEGMENT_MESH& segment_mesh=object.mesh.Get_Segment_Mesh();
    ARRAY<bool> excessive_edge_length(segment_mesh.elements.m);
    int number_of_particles=object.particles.array_collection->Size();
    const T threshold_length_squared=sqr((T)3.1*grid.Minimum_Edge_Length());
    for(int i=1;i<=preimage.meshes.m;i++) level_material_volume(i).Resize(preimage.meshes(i)->elements.m);
    for(int i=1;i<=excessive_edge_length.m;i++){
        int j,k;segment_mesh.elements(i).Get(j,k);
        excessive_edge_length(i)=((object.particles.X(j)-object.particles.X(k)).Magnitude_Squared()>threshold_length_squared);}

    for(int level=1;level<=preimage.meshes.m;level++)for(int t=1;t<=(*preimage.meshes(level)).elements.m;t++){
        int leaf_number=(*preimage.leaf_number(level))(t);
        assert(!level_material_volume(level)(t));
        if(leaf_number){
            if(level<maximum_refinement_depth){
                bool excessive_edge=false;
                for(int j=1;j<=(*object.mesh.element_edges)(leaf_number).Size();j++) excessive_edge|=excessive_edge_length((*object.mesh.element_edges)(leaf_number)(j));
                if(excessive_edge) triangle_list.Append(leaf_number);}
            if(!preimage.Red(level,t)){
                int parent=(*preimage.parent(level))(t);T parent_material_volume=(T)0;
                for(int j=1;j<=T_RED_GREEN_SIMPLICES::number_of_red_children && (*preimage.children(level-1))(parent)(j);j++){
                    int child=(*preimage.children(level-1))(parent)(j);
                    int child_leaf_number=(*preimage.leaf_number(level))(child);
                    parent_material_volume+=simplex_preimage_material_volume(child_leaf_number);
                    level_material_volume(level)(child)=(T)0;
                    simplex_preimage_material_volume(child_leaf_number)=(T)0;}
                level_material_volume(level-1)(parent)+=parent_material_volume;} // NOTE: parents of greenies get processed more than once
            else level_material_volume(level)(t)=simplex_preimage_material_volume(leaf_number);}}

    if(triangle_list.m){preimage.Refine_Simplex_List(triangle_list);preimage.Initialize_Segment_Index_From_Midpoint_Index();} // TODO: may not want to keep calling this
    
    material_particles.Resize(object.particles.array_collection->Size());fixed_particle_list.Resize(object.particles.array_collection->Size());phis.Resize(object.particles.array_collection->Size());
    for(int p=number_of_particles+1;p<=object.particles.array_collection->Size();p++){
        VECTOR<int,2> segment=preimage.segment_mesh.elements((*preimage.segment_index_from_midpoint_index)(p));int i,j;segment.Get(i,j);
        if(phis(i)*phis(j)<0) phis(p)=0;
        else phis(p)=(T).5*(phis(i)+phis(j));
        material_particles(p)=true;
        fixed_particle_list(p)=(fixed_particle_list(i) && fixed_particle_list(j));}

    simplex_preimage_material_volume.Resize(object.mesh.elements.m);
    level_material_volume.Resize(preimage.meshes.m);
    for(int level=1;level<=level_material_volume.m;level++) level_material_volume(level).Resize((*preimage.meshes(level)).elements.m);
    
    ARRAYS_COMPUTATIONS::Fill(simplex_preimage_material_volume,(T)0);
    for(int level=1;level<=level_material_volume.m;level++) for(int i=1;i<=level_material_volume(level).m;i++) if(level_material_volume(level)(i)){
        if(preimage.Leaf(level,i)) simplex_preimage_material_volume((*preimage.leaf_number(level))(i))=level_material_volume(level)(i);
        else{
            T child_volume[T_RED_GREEN_SIMPLICES::number_of_red_children];T parent_volume=0;
            for(int j=0;j<T_RED_GREEN_SIMPLICES::number_of_red_children;j++)child_volume[j]=(T)0;
            T actual_parent_volume=T_SIMPLEX::Signed_Size(object.particles.X.Subset((*preimage.meshes(level)).elements(i)));
            for(int j=1;j<=T_RED_GREEN_SIMPLICES::number_of_red_children && (*preimage.children(level))(i)(j);j++){
                int child=(*preimage.children(level))(i)(j);
                child_volume[j-1]=T_SIMPLEX::Signed_Size(object.particles.X.Subset((*preimage.meshes(level+1)).elements(child)));
                if(child_volume[j-1]*actual_parent_volume<0) child_volume[j-1]*=-1; // should only happen for slivers...
                parent_volume+=child_volume[j-1];}
            T one_over_parent_volume;
            T parent_material_volume=level_material_volume(level)(i);
            if(parent_volume==0) one_over_parent_volume=(T)1;
            else one_over_parent_volume=(T)1/parent_volume;
            for(int j=1;j<=T_RED_GREEN_SIMPLICES::number_of_red_children && (*preimage.children(level))(i)(j);j++){
                int child=(*preimage.children(level))(i)(j);
                T child_material_volume=parent_material_volume*child_volume[j-1]*one_over_parent_volume;
                level_material_volume(level+1)(child)=child_material_volume;}}
        }
}
//#####################################################################
// Function Adjust_Levelset_With_Material_Volumes
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Adjust_Levelset_With_Material_Volumes()
{
    // COMPENSATE FOR DISSIPATION
    //T tolerance=(T)5e-2*grid.Cell_Size(),dtau=(T).1,max_error=0;iteration=1,max_iterations=200;
    //T tolerance=(T)5e-2*grid.Cell_Size(),dtau=(T).1,max_error=0;int iteration=1,max_iterations=10;
    const T dtau=(T).5;
    const T one_over_mean_face_size=(T).5/grid.Face_Size(1); // TODO: make this general and/or base the estimate of face size on the actual interface in cell
    const T max_curvature=1/grid.Minimum_Edge_Length();
    const T half_band_width=particle_levelset.half_band_width;
    Set_Full_Cell_Size(grid.Cell_Size());

    Precompute_Particle_Volumes_In_Cells(removed_negative_particle_volumes_in_cell,particle_levelset.removed_negative_particles,-1);
    Set_Up_For_Refinement();
    particle_levelset.levelset.Compute_Curvature();
    T_ARRAYS_VECTOR gradient;particle_levelset.levelset.Compute_Gradient(gradient);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(abs(particle_levelset.levelset.phi(cell_index))<=half_band_width && abs((*particle_levelset.levelset.curvature)(cell_index))<max_curvature){
            T levelset_volume_of_material_in_cell=Negative_Material(cell_index);
            const T distance_to_move_interface=one_over_mean_face_size*(volume_of_material(cell_index)-levelset_volume_of_material_in_cell);
            const T phi_change_magnitude=dtau*distance_to_move_interface*gradient(cell_index).Magnitude();
            particle_levelset.levelset.phi(cell_index)-=phi_change_magnitude;}}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after level set dissipation correction",0,0);
    particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset.levelset.phi,particle_levelset.levelset.phi,0,0,particle_levelset.number_of_ghost_cells); // TODO: use real dt/time
    LOG::Time("after particles update levelset");
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Volume of material after conservation",0,0);
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step_Particles
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> void VOF_ADVECTION<TV>::
Second_Order_Runge_Kutta_Step_Particles(const T_FACE_LOOKUP& V_lookup,T_ARRAYS_PARTICLES& particles,T_PARTICLE_LEVELSET& particle_levelset,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    T_LINEAR_INTERPOLATION_SCALAR linear_interpolation; // use for second step since particle may not be in initial block

    for(NODE_ITERATOR iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(block_index);
        BLOCK_UNIFORM<GRID<TV> > block(particle_levelset.levelset.grid,block_index);
        T_LINEAR_INTERPOLATION_MAC_HELPER linear_interpolation_mac_helper(block,V_lookup.V_face);
        for(int k=1;k<=cell_particles.array_collection->Size();k++){
            TV velocity=linear_interpolation_mac_helper.Interpolate_Face(cell_particles.X(k));
            TV X_new=cell_particles.X(k)+dt*velocity;
            velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(particle_levelset.levelset.grid,V_lookup,X_new));
            particle_levelset.levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
            cell_particles.X(k)+=dt*velocity;}}}
    particle_levelset.Update_Particle_Cells(particles);
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step_Particles
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> void VOF_ADVECTION<TV>::
Second_Order_Runge_Kutta_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const T_ARRAYS_VECTOR& center_velocities,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    T_LINEAR_INTERPOLATION_VECTOR linear_interpolation; // use for second step since particle may not be in initial block

    for(NODE_ITERATOR iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(block_index);
        for(int k=1;k<=cell_particles.array_collection->Size();k++){
            TV velocity=linear_interpolation.Clamped_To_Array(grid,center_velocities,cell_particles.X(k));
            TV X_new=cell_particles.X(k)+dt*velocity;
            velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array(grid,center_velocities,X_new));
            particle_levelset.levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,dt,time);
            cell_particles.X(k)+=dt*velocity;}}}
    particle_levelset.Update_Particle_Cells(particles);
}
//#####################################################################
// Function Update_Advection_Equation_Cell_Lookup_Postimage
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Perform_Conservative_Advection(const T_FACE_LOOKUP& face_velocities,const T dt,const T time)
{
    T_FAST_LEVELSET& levelset=particle_levelset.levelset;
    const GRID<TV>& grid=levelset.grid;
    if(!volume_of_material_initialized){levelset_advection.Negative_Material(volume_of_material);volume_of_material_initialized=true;Set_Full_Cell_Size(grid.Cell_Size());}
    levelset.boundary->Fill_Ghost_Cells(grid,levelset.phi,levelset.phi,dt,time,particle_levelset.number_of_ghost_cells);

    // rasterize postimage materials
    if(!geometry_initialized){Create_Geometry();geometry_initialized=true;}
    Advect_Material_Preimages(face_velocities,dt,time,fluids_parameters->analytic_test);
    
    if(particle_levelset.mpi_grid){
        // TODO: check what boundary width we really need here
        particle_levelset.mpi_grid->Exchange_Boundary_Cell_Data(volume_of_material,particle_levelset.number_of_ghost_cells,true);}
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& input_directory,const std::string& f)
{
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/vof_object."+f,object);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/geometry_material_volume."+f,simplex_preimage_material_volume);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/cell_postimage_simplices."+f,postimage_simplices_in_cells);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/cell_postimage_phis."+f,phis);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/material_particles."+f,material_particles);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/fixed_particles."+f,fixed_particle_list);
    geometry_initialized=true;
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& output_directory,const std::string& f)
{
#if 0
    if(!debugging && simplex_preimage_material_volume.m==object.mesh.elements.m){
        ARRAY<OPENGL_COLOR > color_map(object.mesh.elements.m);
        for(int i=1;i<=object.mesh.elements.m;i++){
            T size=abs(T_SIMPLEX::Signed_Size(object.particles.X.Subset(object.mesh.elements(i))));
            if(size){
                T ratio=simplex_preimage_material_volume(i)/size;
                color_map(i)=OPENGL_COLOR(min((T)1,max((T)0,ratio-(T)1)),pow((T)1-min((T)1,sqr((T)1-ratio)),(T)3),min((T)1,max((T)0,(T)1-ratio)));}
            else color_map(i)=OPENGL_COLOR::Green();}
        int max_count=0;TV_INT max_index;
        for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(simplices_in_cells(cell_index).m>max_count){
                max_count=simplices_in_cells(cell_index).m;max_index=cell_index;}}
                FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/vof_colors."+f,color_map);}
    Write_VOF_Object(stream_type,output_directory+"/vof_object."+f,object); // TODO: fix this for restarts in 3D
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/geometry_material_volume."+f,simplex_preimage_material_volume);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/cell_postimage_simplices."+f,postimage_simplices_in_cells);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/cell_postimage_phis."+f,phis);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/material_particles."+f,material_particles);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/fixed_particles."+f,fixed_particle_list);
#endif
}
//#####################################################################
// Function Create_Initial_Meshing_For_Cell
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Create_Initial_Meshing_For_Cell(GRID<TV>& grid,TRIANGULATED_AREA<T>& object,const VECTOR<int,2>& cell_index,const ARRAY<int,VECTOR<int,2> >& nodes_to_particles_map,ARRAY<VECTOR<int,2> >& simplex_cells)
{
    VECTOR<int,2> nodes_in_cell[4];grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,nodes_in_cell);
    int particle_index=object.particles.array_collection->Add_Element();
    VECTOR<T,2> center;
    for(int i=0;i<4;i++) center+=object.particles.X(nodes_to_particles_map(nodes_in_cell[i]));
    center*=(T).25;object.particles.X(particle_index)=center;
    for(int i=1;i<=4;i++) simplex_cells.Append(cell_index);
    fixed_particle_list.Append(false);
    object.mesh.elements.Append(VECTOR<int,3>(nodes_to_particles_map(nodes_in_cell[0]),nodes_to_particles_map(nodes_in_cell[1]),particle_index));
    object.mesh.elements.Append(VECTOR<int,3>(nodes_to_particles_map(nodes_in_cell[1]),nodes_to_particles_map(nodes_in_cell[3]),particle_index));
    object.mesh.elements.Append(VECTOR<int,3>(nodes_to_particles_map(nodes_in_cell[0]),particle_index,nodes_to_particles_map(nodes_in_cell[2])));
    object.mesh.elements.Append(VECTOR<int,3>(particle_index,nodes_to_particles_map(nodes_in_cell[3]),nodes_to_particles_map(nodes_in_cell[2])));
}
//#####################################################################
// Function Create_Initial_Meshing_For_Cell
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Create_Initial_Meshing_For_Cell(GRID<TV>& grid,TETRAHEDRALIZED_VOLUME<T>& object,const VECTOR<int,3>& cell_index,const ARRAY<int,VECTOR<int,3> >& nodes_to_particles_map,ARRAY<VECTOR<int,3> >& simplex_cells)
{
    VECTOR<int,3> nodes_in_cell[8];grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell_index,nodes_in_cell);
    for(int i=1;i<=5;i++) simplex_cells.Append(cell_index);
    if(cell_index.Sum()%2){
        // corner tets
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[1]),nodes_to_particles_map(nodes_in_cell[5]),nodes_to_particles_map(nodes_in_cell[3]),nodes_to_particles_map(nodes_in_cell[0])));
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[2]),nodes_to_particles_map(nodes_in_cell[6]),nodes_to_particles_map(nodes_in_cell[0]),nodes_to_particles_map(nodes_in_cell[3])));
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[4]),nodes_to_particles_map(nodes_in_cell[5]),nodes_to_particles_map(nodes_in_cell[0]),nodes_to_particles_map(nodes_in_cell[6])));
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[7]),nodes_to_particles_map(nodes_in_cell[5]),nodes_to_particles_map(nodes_in_cell[6]),nodes_to_particles_map(nodes_in_cell[3])));
        // inside tet
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[0]),nodes_to_particles_map(nodes_in_cell[6]),nodes_to_particles_map(nodes_in_cell[5]),nodes_to_particles_map(nodes_in_cell[3])));}
    else{
        // corner tets
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[0]),nodes_to_particles_map(nodes_in_cell[2]),nodes_to_particles_map(nodes_in_cell[4]),nodes_to_particles_map(nodes_in_cell[1])));
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[3]),nodes_to_particles_map(nodes_in_cell[7]),nodes_to_particles_map(nodes_in_cell[2]),nodes_to_particles_map(nodes_in_cell[1])));
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[5]),nodes_to_particles_map(nodes_in_cell[4]),nodes_to_particles_map(nodes_in_cell[7]),nodes_to_particles_map(nodes_in_cell[1])));
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[6]),nodes_to_particles_map(nodes_in_cell[4]),nodes_to_particles_map(nodes_in_cell[2]),nodes_to_particles_map(nodes_in_cell[7])));

        // inside tet
        object.mesh.elements.Append(VECTOR<int,4>(nodes_to_particles_map(nodes_in_cell[1]),nodes_to_particles_map(nodes_in_cell[4]),nodes_to_particles_map(nodes_in_cell[7]),nodes_to_particles_map(nodes_in_cell[2])));}
}
//#####################################################################
// Function Refine_Simplex
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Refine_Simplex(ARRAY<VECTOR<int,3> >& tris,ARRAY<VECTOR<T,2> >& particle_X,const VECTOR<int,3>& indices)
{
    const int x1=indices[1],x2=indices[2],x3=indices[3];
    const VECTOR<T,2> X1=particle_X(x1),X2=particle_X(x2),X3=particle_X(x3);
    int first_midpoint=particle_X.Append((T).5*(X1+X2));
    int second_midpoint=particle_X.Append((T).5*(X2+X3));
    int third_midpoint=particle_X.Append((T).5*(X3+X1));
    tris.Append(VECTOR<int,3>(x1,first_midpoint,third_midpoint));
    tris.Append(VECTOR<int,3>(x2,second_midpoint,first_midpoint));
    tris.Append(VECTOR<int,3>(x3,third_midpoint,second_midpoint));
    tris.Append(VECTOR<int,3>(first_midpoint,second_midpoint,third_midpoint));
}
//#####################################################################
// Function Refine_Simplex
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Refine_Simplex(ARRAY<VECTOR<int,4> >& tets,ARRAY<VECTOR<T,3> >& particle_X,const VECTOR<int,4>& indices)
{
    const int x1=indices[1],x2=indices[2],x3=indices[3],x4=indices[4];
    const VECTOR<T,3> X1=particle_X(x1),X2=particle_X(x2),X3=particle_X(x3),X4=particle_X(x4);
    int first_midpoint=particle_X.Append((T).5*(X1+X2));
    int second_midpoint=particle_X.Append((T).5*(X1+X3));
    int third_midpoint=particle_X.Append((T).5*(X1+X4));
    int fourth_midpoint=particle_X.Append((T).5*(X2+X3));
    int fifth_midpoint=particle_X.Append((T).5*(X2+X4));
    int sixth_midpoint=particle_X.Append((T).5*(X3+X4));
    tets.Append(VECTOR<int,4>(x1,first_midpoint,second_midpoint,third_midpoint));
    tets.Append(VECTOR<int,4>(x2,fourth_midpoint,first_midpoint,fifth_midpoint));
    tets.Append(VECTOR<int,4>(x3,second_midpoint,fourth_midpoint,sixth_midpoint));
    tets.Append(VECTOR<int,4>(x4,third_midpoint,sixth_midpoint,fifth_midpoint));
    // pick our octagon diagonal as first_midpoint,sixth_midpoint
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,third_midpoint,fifth_midpoint));
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,fifth_midpoint,fourth_midpoint));
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,fourth_midpoint,second_midpoint));
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,second_midpoint,third_midpoint));
}
//#####################################################################
// Function Find_Fixed_Cell_On_Point
//#####################################################################
template<class TV> bool VOF_ADVECTION<TV>::
Find_Fixed_Cell_On_Point(const TV& point,TV_INT& opposing_cell,TV& face_normal,const T length_scale) const
{
    TV_INT cell=grid.Cell(point,3);
    for(int axis=1;axis<=GRID<TV>::dimension;axis++){
        if(grid.Face_Domain(axis,cell,grid.Minimum_Edge_Length()*length_scale).Lazy_Inside(point)){
            if(fixed_cells(cell)) opposing_cell=cell;else opposing_cell=cell-TV_INT::Axis_Vector(axis);
            face_normal=(grid.Face(axis,cell)-grid.Center(opposing_cell)).Normalized();
            return fixed_cells(opposing_cell);}
        if(grid.Face_Domain(axis,cell+TV_INT::Axis_Vector(axis),grid.Minimum_Edge_Length()*length_scale).Lazy_Inside(point)){
            if(fixed_cells(cell)) opposing_cell=cell;else opposing_cell=cell+TV_INT::Axis_Vector(axis);
            face_normal=(grid.Face(axis,cell+TV_INT::Axis_Vector(axis))-grid.Center(opposing_cell)).Normalized();
            return fixed_cells(opposing_cell);}}
    return false;
}
//#####################################################################
// Function Set_Up_For_Refinement
//#####################################################################
template<class TV> void VOF_ADVECTION<TV>::
Set_Up_For_Refinement()
{
    PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>* pls_implicit_object=dynamic_cast<PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>*>(implicit_object);
    if(pls_implicit_object && fluids_parameters) pls_implicit_object->Precompute_Cell_Particle_Influence();
    // also precompute nodal phis
    T_LEVELSET& levelset=implicit_object->levelset;
    grid_nodal_phis.Resize(levelset.grid.Node_Indices(particle_levelset.number_of_ghost_cells));
    for(NODE_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next())
        grid_nodal_phis(iterator.Node_Index())=Phi(iterator.Location());
}
//#####################################################################
template class VOF_ADVECTION<VECTOR<float,2> >;
template class VOF_ADVECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VOF_ADVECTION<VECTOR<double,2> >;
template class VOF_ADVECTION<VECTOR<double,3> >;
#endif
