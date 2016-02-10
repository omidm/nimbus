//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_THINSHELL_ADVECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Computations/BOX_BOX_INTERSECTION_AREA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/EULER_THINSHELL_ADVECTION.h>

namespace PhysBAM {
//#####################################################################
// Construct_Sparce_Indices
//#####################################################################
namespace{
template<class T,int d>
int Construct_Sparce_Indices(const GRID<VECTOR<T,d> >& grid, const ARRAY<bool,VECTOR<int,d> >& psi, const ARRAY<bool,VECTOR<int,d> >& near_interface, const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,
                              HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices)
{
    typedef VECTOR<int,d> TV_INT;typedef typename GRID<VECTOR<T,d> >::CELL_ITERATOR CELL_ITERATOR;
    int num_indices=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const TV_INT cell_index=iterator.Cell_Index();
        if(near_interface(cell_index) && psi(cell_index) && !cut_cells(cell_index))
            indices.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,0))=++num_indices;
        if(cut_cells(cell_index)) for(int poly=1;poly<=cut_cells(cell_index)->geometry.Size();++poly){
            if(cut_cells(cell_index)->visibility(poly).Size() && psi(cut_cells(cell_index)->visibility(poly)(1))) indices.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,poly))=++num_indices;}}

    return num_indices;
}
//#####################################################################
// Output_Stuff
//#####################################################################
template<class T,int d>
void Output_Stuff(HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices,ARRAY<T>& volume,ARRAY<VECTOR<T,d+2> >& state)
{
    std::stringstream ss;
    for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices);iter.Valid();iter.Next())
        ss<<"Stuff in "<<iter.Data()<<" = "<<volume(iter.Data())<<" * "<<state(iter.Data())<<" => "<<volume(iter.Data())*state(iter.Data())<<std::endl;
    LOG::filecout(ss.str());
}
//#####################################################################
// Fill_Volume
//#####################################################################
template<class T,int d>
void Fill_Volume(const GRID<VECTOR<T,d> >& grid,HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,ARRAY<T>& volume)
{
    for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices);iter.Valid();iter.Next())
        if(!iter.Key().y) volume(iter.Data())=grid.Cell_Size();
        else volume(iter.Data())=cut_cells(iter.Key().x)->geometry(iter.Key().y).Area();
}
//#####################################################################
// Fill_State
//#####################################################################
template<class T,int d>
void Fill_State(const GRID<VECTOR<T,d> >& grid,HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,const ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U,ARRAY<VECTOR<T,d+2> >& state)
{
    for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices);iter.Valid();iter.Next()){
        if(!iter.Key().y){
            state(iter.Data())=U(iter.Key().x);}
        else{
            if(cut_cells(iter.Key().x)->dominant_element==iter.Key().y){
                state(iter.Data())=U(iter.Key().x);}
            else{const int num_neighbors=cut_cells(iter.Key().x)->visibility(iter.Key().y).Size();
                if(!num_neighbors){std::stringstream ss;ss<<"ERROR: Trying to fill cell "<<iter.Key().x<<", poly "<<iter.Key().y<<" with data, but no visible neighbors exist!"<<std::endl;LOG::filecout(ss.str());continue;}
                for(int i=1;i<=num_neighbors;++i) state(iter.Data())+=U(cut_cells(iter.Key().x)->visibility(iter.Key().y)(i));
                state(iter.Data()) *= (T)1/num_neighbors;}}}
}
//#####################################################################
// Restore_State
//#####################################################################
template<class T,int d>
void Restore_State(const GRID<VECTOR<T,d> >& grid,HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,const ARRAY<bool,VECTOR<int,d> >& near_interface,
                   const ARRAY<bool,VECTOR<int,d> >& psi,const ARRAY<T,VECTOR<int,d> >& cell_volumes,const ARRAY<VECTOR<T,d+2> >& state,const ARRAY<T>& volume,ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U)
{
    typedef VECTOR<int,d> TV_INT;typedef typename GRID<VECTOR<T,d> >::CELL_ITERATOR CELL_ITERATOR;
    for(typename HASHTABLE<PAIR<TV_INT,int>,int>::ITERATOR iter(indices);iter.Valid();iter.Next())
        if(near_interface(iter.Key().x) && psi(iter.Key().x)) U(iter.Key().x)=VECTOR<T,d+2>();

    for(typename HASHTABLE<PAIR<TV_INT,int>,int>::ITERATOR iter(indices);iter.Valid();iter.Next())
        if(!iter.Key().y) U(iter.Key().x)+=state(iter.Data());
        else{
            if(cut_cells(iter.Key().x)->dominant_element==iter.Key().y){
                if(!psi(iter.Key().x)) {std::stringstream ss;ss<<"ERROR: Trying to distribute "<<state(iter.Data())<<" from "<<iter.Data()<<" to cell "<<iter.Key().x<<" Poly "<<iter.Key().y<<std::endl;LOG::filecout(ss.str());}
                U(iter.Key().x)+=state(iter.Data());}
            else{const int num_neighbors=cut_cells(iter.Key().x)->visibility(iter.Key().y).Size();
                if(!num_neighbors){std::stringstream ss;ss<<"WARNING: Trying to restore cell "<<iter.Key().x<<", poly "<<iter.Key().y<<" with data, but no visible neighbors exist!"<<std::endl;LOG::filecout(ss.str());continue;}
                const T one_over_num_neighbors=(T)1/num_neighbors;
                for(int i=1;i<=num_neighbors;++i){
                    if(!psi(cut_cells(iter.Key().x)->visibility(iter.Key().y)(i))){
                        std::stringstream ss;
                        ss<<"ERROR: Trying to distribute "<<state(iter.Data())<<" * "<<one_over_num_neighbors<<" from "<<iter.Data()<<" to cell "
                            <<cut_cells(iter.Key().x)->visibility(iter.Key().y)(i)<<std::endl;
                        LOG::filecout(ss.str());}
                    U(cut_cells(iter.Key().x)->visibility(iter.Key().y)(i))+=state(iter.Data())*one_over_num_neighbors;}}}

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const TV_INT cell_index=iterator.Cell_Index();
        if(near_interface(cell_index) && psi(cell_index)){assert(cell_volumes(cell_index)>(T)1e-4*grid.Cell_Size());
            if(cell_volumes(cell_index) < 1e-12){
#if 0
                std::stringstream ss;
                ss<<"Trying to divide by really small number; "<<cell_volumes(cell_index)<<"\t"<<cut_cells(cell_index)->dominant_element<<"\t"<<psi(cell_index)<<std::endl;
                LOG::filecout(ss.str());
#endif
                CUT_CELLS<T,d>::Print_Debug_Information(cut_cells(cell_index)->geometry);}
            U(cell_index) /= cell_volumes(cell_index);}}
}
//#####################################################################
// Compute_Velocity_Field
//#####################################################################
template<class T,int d>
void Compute_Velocity_Field(HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices_n,  const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_n,  const ARRAY<VECTOR<T,d+2> >& state_n,  
                            HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices_np1,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_np1,const ARRAY<bool,VECTOR<int,d> >& swept_cells,
                            const ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U_n,ARRAY<VECTOR<T,d> >& velocity_n,ARRAY<VECTOR<T,d> >& velocity_np1)
{
    typedef VECTOR<int,d> TV_INT;
    for(typename HASHTABLE<PAIR<TV_INT,int>,int>::ITERATOR iter(indices_n);  iter.Valid();iter.Next()){
        velocity_n(iter.Data())=EULER<GRID<VECTOR<T,d> > >::Get_Velocity(state_n(iter.Data()));}

    for(typename HASHTABLE<PAIR<TV_INT,int>,int>::ITERATOR iter(indices_np1);iter.Valid();iter.Next()){
        if(!iter.Key().y){
            velocity_np1(iter.Data())=EULER<GRID<VECTOR<T,d> > >::Get_Velocity(U_n(iter.Key().x));}
        else{
            if(cut_cells_np1(iter.Key().x)->dominant_element==iter.Key().y && !swept_cells(iter.Key().x)){
                velocity_np1(iter.Data())=EULER<GRID<VECTOR<T,d> > >::Get_Velocity(U_n(iter.Key().x));}
            else{int num_neighbors=cut_cells_np1(iter.Key().x)->visibility(iter.Key().y).Size();
                for(int i=num_neighbors;i>=1;--i)
                    if(swept_cells(cut_cells_np1(iter.Key().x)->visibility(iter.Key().y)(i))) --num_neighbors;
                    else{velocity_np1(iter.Data())+=EULER<GRID<VECTOR<T,d> > >::Get_Velocity(U_n(cut_cells_np1(iter.Key().x)->visibility(iter.Key().y)(i)));}
                if(!num_neighbors){std::stringstream ss;ss<<"ERROR: No visible neighbors for cell "<<iter.Key().x<<" cut cell "<<iter.Key().y<<std::endl;LOG::filecout(ss.str());continue;}
                assert(num_neighbors);velocity_np1(iter.Data()) *= (T)1/num_neighbors;}}}
}
};
//#####################################################################
// Compute_Hybrid_Boundary_Fluxes
//#####################################################################
template<class T,int d>
void Compute_Hybrid_Boundary_Fluxes(const GRID<VECTOR<T,d> >& grid,const T dt,const ARRAY<bool,VECTOR<int,d> >& near_interface_mask,const ARRAY<bool,VECTOR<int,d> >& psi,
                                    const ARRAY<VECTOR<T,d+2>,FACE_INDEX<d> >& flux_boundary_conditions,ARRAY<TRIPLE<VECTOR<int,d>,VECTOR<int,d>,VECTOR<T,d+2> > >& hybrid_flux_data)
{
    typedef VECTOR<int,d> TV_INT;typedef typename GRID<VECTOR<T,d> >::FACE_ITERATOR FACE_ITERATOR;
    const VECTOR<T,d> flux_face_size=dt*grid.Cell_Size()*grid.One_Over_DX();
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const TV_INT& first_cell_index=iterator.First_Cell_Index();const TV_INT& second_cell_index=iterator.Second_Cell_Index();
        if(!(near_interface_mask(first_cell_index) ^ near_interface_mask(second_cell_index))
           || (near_interface_mask(first_cell_index ) && psi.Valid_Index(second_cell_index) && !psi(second_cell_index))
           || (near_interface_mask(second_cell_index) && psi.Valid_Index(first_cell_index ) && !psi(first_cell_index ))) continue;
        hybrid_flux_data.Append(TRIPLE<TV_INT,TV_INT,VECTOR<T,d+2> >(first_cell_index,second_cell_index,flux_face_size(iterator.Axis())*flux_boundary_conditions(iterator.Full_Index())));}
}
//#####################################################################
// Advect_Near_Interface_Data
//#####################################################################
namespace {
template<class T> void
Add_Weight_To_Advection(const T weight,const int donor_index,const int receiver_index,ARRAY<T>& sigma,ARRAY<ARRAY<int> >& donors,ARRAY<ARRAY<int> >& receivers,ARRAY<PAIR<T,int> >& weights)
{
    int index=0;
    if(donor_index){sigma(donor_index)+=weight;
        for(int i=1;i<=donors(donor_index).Size();++i) if(weights(donors(donor_index)(i)).y==receiver_index) index=donors(donor_index)(i);}

    if(!index && receiver_index){
        index=weights.Append(PAIR<T,int>(weight,receiver_index));
        receivers(receiver_index).Append(index);
        if(donor_index) donors(donor_index).Append(index);}
    else if(index) weights(index).x += weight;

#if 0
    std::stringstream ss;
    ss<<"ADVECTION "<<index<<": Adding weight "<<weight<<" from "<<donor_index<<" to "<<receiver_index<<"; sigma("<<donor_index<<") = "<<(donor_index ? sigma(donor_index) : 0)<<std::endl;
    LOG::filecout(ss.str());
#endif
}
template<class T,int d> int
Get_Backward_Visible_Neighbor_Index(const VECTOR<int,d>& cell_index,const int cell_polygon,const VECTOR<int,d>& neighbor_index,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices,HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >& visibility_graph,int* poly)
{
    int donor=0;
    if(!cut_cells(neighbor_index)){
        if(neighbor_index==cell_index) indices.Get(PAIR<VECTOR<int,d>,int>(cell_index,0),donor);
        else if(visibility_graph.Contains(PAIR<VECTOR<int,d>,int>(neighbor_index,0)) && visibility_graph.Get(PAIR<VECTOR<int,d>,int>(neighbor_index,0)).Contains(PAIR<VECTOR<int,d>,int>(cell_index,cell_polygon))){
                indices.Get(PAIR<VECTOR<int,d>,int>(neighbor_index,0),donor);if(poly) *poly=0;}}
    else{
        for(int neighbor_poly=1;neighbor_poly<=cut_cells(neighbor_index)->geometry.Size();++neighbor_poly){
            if(visibility_graph.Contains(PAIR<VECTOR<int,d>,int>(neighbor_index,neighbor_poly)) && visibility_graph.Get(PAIR<VECTOR<int,d>,int>(neighbor_index,neighbor_poly)).Contains(PAIR<VECTOR<int,d>,int>(cell_index,cell_polygon))){
                indices.Get(PAIR<VECTOR<int,d>,int>(neighbor_index,neighbor_poly),donor);if(poly) *poly=neighbor_poly;}}}
    return donor;
}
template<class T,int d> int
Get_Forward_Visible_Neighbor_Index(const VECTOR<int,d>& cell_index,const int cell_polygon,const VECTOR<int,d>& neighbor_index,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells,HASHTABLE<PAIR<VECTOR<int,d>,int>,int>& indices,HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >& visibility_graph)
{
    int donor=0;
    if(!cut_cells(neighbor_index)){
        if(neighbor_index==cell_index) indices.Get(PAIR<VECTOR<int,d>,int>(cell_index,0),donor);
        else if(visibility_graph.Contains(PAIR<VECTOR<int,d>,int>(cell_index,cell_polygon)) && visibility_graph.Get(PAIR<VECTOR<int,d>,int>(cell_index,cell_polygon)).Contains(PAIR<VECTOR<int,d>,int>(neighbor_index,0)))
            indices.Get(PAIR<VECTOR<int,d>,int>(neighbor_index,0),donor);}
    else{
        if(!visibility_graph.Contains(PAIR<VECTOR<int,d>,int>(cell_index,cell_polygon))) return 0;
        const ARRAY<PAIR<VECTOR<int,d>,int> >& neighbor_polygons(visibility_graph.Get(PAIR<VECTOR<int,d>,int>(cell_index,cell_polygon)));
        for(int i=1;i<=neighbor_polygons.Size();++i){
            if(neighbor_polygons(i).x==neighbor_index){indices.Get(PAIR<VECTOR<int,d>,int>(neighbor_polygons(i)),donor);}}}
    return donor;
}
};
template<class T,int d> void Advect_Near_Interface_Data(
        const GRID<VECTOR<T,d> >& grid,const T collision_thickness,const T dt,const ARRAY<bool,VECTOR<int,d> >& near_interface_mask,
        const ARRAY<bool,VECTOR<int,d> >& swept_cells,const ARRAY<VECTOR<T,d+2>,FACE_INDEX<d> >& flux_boundary_conditions,const ARRAY<T,VECTOR<int,d> >& cell_volumes_np1,
        HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >& visibility_graph,
        const ARRAY<bool,VECTOR<int,d> >& psi_n,  const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_n,  const ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U_n,
        const ARRAY<bool,VECTOR<int,d> >& psi_np1,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_np1,      ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >& U_np1)
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;typedef VECTOR<T,d+2> TV_DIMENSION;typedef GRID<VECTOR<T,d> > T_GRID;typedef typename GRID<VECTOR<T,d> >::CELL_ITERATOR CELL_ITERATOR;

    HASHTABLE<PAIR<TV_INT,int>,int> indices_n,indices_np1;
    int num_indices_n,num_indices_np1;
    {num_indices_n  =Construct_Sparce_Indices<T,d>(grid,psi_n  ,near_interface_mask,cut_cells_n  ,indices_n  );}
    {num_indices_np1=Construct_Sparce_Indices<T,d>(grid,psi_np1,near_interface_mask,cut_cells_np1,indices_np1);}

    ARRAY<T> volumes_n(num_indices_n),volumes_np1(num_indices_np1);
    {Fill_Volume<T,d>(grid,indices_n,  cut_cells_n,  volumes_n);Fill_Volume<T,d>(grid,indices_np1,cut_cells_np1,volumes_np1);}

    ARRAY<TV_DIMENSION> state_n(num_indices_n),state_np1(num_indices_np1);
    {ARRAYS_COMPUTATIONS::Fill(state_n,  TV_DIMENSION());ARRAYS_COMPUTATIONS::Fill(state_np1,TV_DIMENSION());Fill_State<T,d>(grid,indices_n,cut_cells_n,U_n,state_n);}

    ARRAY<TV> velocity_n(num_indices_n),velocity_np1(num_indices_np1);
    {Compute_Velocity_Field<T,d>(indices_n,cut_cells_n,state_n,indices_np1,cut_cells_np1,swept_cells,U_n,velocity_n,velocity_np1);}

    TV_DIMENSION hybrid_flux_sum=TV_DIMENSION();
    ARRAY<TRIPLE<TV_INT,TV_INT,TV_DIMENSION> > hybrid_flux_data;
    {Compute_Hybrid_Boundary_Fluxes<T,d>(grid,dt,near_interface_mask,psi_np1,flux_boundary_conditions,hybrid_flux_data);}

    for(int variable_index=1;variable_index<=d+2;++variable_index){
        ARRAY<T> sigma(num_indices_n),cell_stuff(num_indices_n);ARRAY<PAIR<T,int> > weights;
        ARRAY<ARRAY<int> > donors(num_indices_n),receivers(num_indices_np1);

        for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices_n);iter.Valid();iter.Next())
            cell_stuff(iter.Data()) = state_n(iter.Data())(variable_index)*volumes_n(iter.Data());

        {T accumulated_material=0;
        for(int i=1;i<=hybrid_flux_data.Size();++i){
            if(hybrid_flux_data(i).z(variable_index) >= 0){
                int first_index=0;  indices_n.Get(PAIR<TV_INT,int>(hybrid_flux_data(i).x,  (cut_cells_n.Valid_Index(hybrid_flux_data(i).x)   && cut_cells_n(hybrid_flux_data(i).x)   ? cut_cells_n(hybrid_flux_data(i).x)->dominant_element   : 0)),first_index);
                int second_index=0; indices_np1.Get(PAIR<TV_INT,int>(hybrid_flux_data(i).y,(cut_cells_np1.Valid_Index(hybrid_flux_data(i).y) && cut_cells_np1(hybrid_flux_data(i).y) ? cut_cells_np1(hybrid_flux_data(i).y)->dominant_element : 0)),second_index);
                if(first_index){cell_stuff(first_index) -= hybrid_flux_data(i).z(variable_index);sigma(first_index) -= hybrid_flux_data(i).z(variable_index);}
                if(first_index && second_index) {std::stringstream ss;ss<<"Weird: "<<first_index<<" -> "<<second_index<<"; "<<hybrid_flux_data(i).z(variable_index)<<std::endl;LOG::filecout(ss.str());}
                if(first_index) accumulated_material -= hybrid_flux_data(i).z(variable_index);
                if(second_index) accumulated_material += hybrid_flux_data(i).z(variable_index);
                Add_Weight_To_Advection(hybrid_flux_data(i).z(variable_index), first_index, second_index, sigma, donors, receivers, weights);}
            else{
                int first_index=0;  indices_n.Get(PAIR<TV_INT,int>(hybrid_flux_data(i).y,  (cut_cells_n.Valid_Index(hybrid_flux_data(i).y)   && cut_cells_n(hybrid_flux_data(i).y)   ? cut_cells_n(hybrid_flux_data(i).y)->dominant_element   : 0)),first_index);
                int second_index=0; indices_np1.Get(PAIR<TV_INT,int>(hybrid_flux_data(i).x,(cut_cells_np1.Valid_Index(hybrid_flux_data(i).x) && cut_cells_np1(hybrid_flux_data(i).x) ? cut_cells_np1(hybrid_flux_data(i).x)->dominant_element : 0)),second_index);
                if(first_index){cell_stuff(first_index) += hybrid_flux_data(i).z(variable_index);sigma(first_index) += hybrid_flux_data(i).z(variable_index);}
                if(first_index && second_index) {std::stringstream ss;ss<<"Weird: "<<first_index<<" -> "<<second_index<<"; "<<hybrid_flux_data(i).z(variable_index)<<std::endl;LOG::filecout(ss.str());}
                if(first_index) accumulated_material += hybrid_flux_data(i).z(variable_index);
                if(second_index) accumulated_material -= hybrid_flux_data(i).z(variable_index);
                Add_Weight_To_Advection(abs(hybrid_flux_data(i).z(variable_index)), first_index, second_index, sigma, donors, receivers, weights);}}
        hybrid_flux_sum(variable_index) = accumulated_material;
        }

        {for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices_np1);iter.Valid();iter.Next()){
                const int receiver_index=iter.Data();
                int donor_index;T volume_overlap=0;
                const TV offset=-dt*velocity_np1(receiver_index);
                RANGE<TV> bounding_box;
                if(cut_cells_np1(iter.Key().x)) bounding_box=RANGE<TV>::Bounding_Box(cut_cells_np1(iter.Key().x)->geometry(iter.Key().y).X)+offset;
                else bounding_box=RANGE<TV>(grid.Node(iter.Key().x)+offset,grid.Node(iter.Key().x+TV_INT::All_Ones_Vector())+offset);
#if 0
                // LOG::cout<<"Looking backward from Cell "<<iter.Key().x<<", Polygon "<<iter.Key().y<<" ["<<iter.Data()<<"] by velocity "<<velocity_np1(receiver_index)<<": "<<bounding_box<<std::endl;
#endif
                RANGE<TV_INT> affected_cells(grid.Clamp_To_Cell(bounding_box.min_corner-collision_thickness),grid.Clamp_To_Cell(bounding_box.max_corner+collision_thickness));
                for(CELL_ITERATOR intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                    int poly=0;
                    donor_index=Get_Backward_Visible_Neighbor_Index<T,d>(iter.Key().x,iter.Key().y,intersecting_iter.Cell_Index(),cut_cells_n,indices_n,visibility_graph,&poly);
                    if(donor_index){
                        RANGE<TV> intersecting_bounding_box=(poly ? RANGE<TV>::Bounding_Box(cut_cells_n(intersecting_iter.Cell_Index())->geometry(poly).X) : intersecting_iter.Bounding_Box());
                        volume_overlap=INTERSECTION::Intersection_Area(bounding_box,intersecting_bounding_box);
                        Add_Weight_To_Advection(volume_overlap*state_n(donor_index)(variable_index), donor_index, receiver_index, sigma, donors, receivers, weights);}}}
                // MPI: Move w_ij's to correct MPI domain
                // MPI: (Re)-compute sigma_i
        }

        {for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices_n);iter.Valid();iter.Next()){
                const int donor_index=iter.Data();
                if(abs(sigma(donor_index)) > abs(cell_stuff(iter.Data()))){
                    const T scale=cell_stuff(iter.Data())/sigma(donor_index);
                    for(int i=1;i<=donors(donor_index).Size();++i) weights(donors(donor_index)(i)).x *= scale;}
                else{
                    const T remainder=cell_stuff(iter.Data())-sigma(donor_index);
                    int receiver_index;T volume_overlap=0;
                    const TV offset=dt*velocity_n(donor_index);
#if 0
                    std::stringstream ss;
                    ss<<"Looking forward from Cell "<<iter.Key().x<<", Polygon "<<iter.Key().y<<" by velocity "<<velocity_n(donor_index)<<std::endl;
                    LOG::filecout(ss.str());
#endif
                    RANGE<TV> bounding_box;
                    if(cut_cells_n(iter.Key().x)){bounding_box=RANGE<TV>::Bounding_Box(cut_cells_n(iter.Key().x)->geometry(iter.Key().y).X)+offset;}
                    else{bounding_box=RANGE<TV>(grid.Node(iter.Key().x)+offset,grid.Node(iter.Key().x+TV_INT::All_Ones_Vector())+offset);}

                    ARRAY<PAIR<int,T> > forward_weights;T distributed_volume=0;
                    RANGE<TV_INT> affected_cells(grid.Clamp_To_Cell(bounding_box.min_corner-collision_thickness),grid.Clamp_To_Cell(bounding_box.max_corner+collision_thickness));
                    for(CELL_ITERATOR intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                        int poly=0;
                        receiver_index=Get_Forward_Visible_Neighbor_Index<T,d>(iter.Key().x,iter.Key().y,intersecting_iter.Cell_Index(),cut_cells_np1,indices_np1,visibility_graph);
                        if(receiver_index){
                            RANGE<TV> intersecting_bounding_box=(poly ? RANGE<TV>::Bounding_Box(cut_cells_np1(intersecting_iter.Cell_Index())->geometry(poly).X) : intersecting_iter.Bounding_Box());
                            volume_overlap=INTERSECTION::Intersection_Area(bounding_box,intersecting_bounding_box);
                            distributed_volume+=volume_overlap;
                            forward_weights.Append(PAIR<int,T>(receiver_index,volume_overlap));}}
                    if(distributed_volume < (T)1e-8*grid.Cell_Size()){
                        if(donors(donor_index).Size()){const T one_over_num_neighbors=(T)1/donors(donor_index).Size();
                            for(int i=1;i<=donors(donor_index).Size();++i) weights(donors(donor_index)(i)).x += remainder*one_over_num_neighbors;}
                        else{
                            std::stringstream ss;
                            ss<<"WARNING: Stuff in "<<donor_index<<" [Cell "<<iter.Key().x<<", Poly "<<iter.Key().y<<"] just can't be put ANYWHERE; hacking together a solution...!"<<std::endl;
                            bool found_target_distributor=false;
                            if(cut_cells_n(iter.Key().x)) for(int i=1;!found_target_distributor && i<=cut_cells_n(iter.Key().x)->visibility(iter.Key().y).Size();++i){
                                const TV_INT& neighbor_index(cut_cells_n(iter.Key().x)->visibility(iter.Key().y)(i));
                                int neighbor_poly=0;if(cut_cells_n(neighbor_index)) neighbor_poly=cut_cells_n(neighbor_index)->dominant_element;
                                int pseudo_donor=0;indices_n.Get(PAIR<TV_INT,int>(neighbor_index,neighbor_poly),pseudo_donor);
                                if(pseudo_donor && donors(pseudo_donor).Size()){const T one_over_num_neighbors=(T)1/donors(pseudo_donor).Size();
                                    for(int i=1;i<=donors(pseudo_donor).Size();++i) weights(donors(pseudo_donor)(i)).x += remainder*one_over_num_neighbors;found_target_distributor=true;}}
                            else ss<<"ERROR: Cell "<<iter.Key().x<<" Poly "<<iter.Key().y<<" Has no temporal visible neighbors, and isn't even a cut cell..."<<std::endl;
                            if(!found_target_distributor) ss<<"ERROR:\tThe above warning could not be hacked together!"<<std::endl;
                            LOG::filecout(ss.str());}}
                    else{
                        const T dilation_factor=remainder/distributed_volume;
                        for(int i=1;i<=forward_weights.Size();++i) Add_Weight_To_Advection(dilation_factor*forward_weights(i).y, donor_index, forward_weights(i).x, sigma, donors, receivers, weights);}}}
                // MPI: Move f_ij to proper receiving indices
        }

        for(int i=1;i<=weights.Size();++i){
            state_np1(weights(i).y)(variable_index) += weights(i).x;}
    }

    Restore_State<T,d>(grid,indices_np1,cut_cells_np1,near_interface_mask,psi_np1,cell_volumes_np1,state_np1,volumes_np1,U_np1);

#if 1 // Debug validation and verification
    TV_DIMENSION accumulated_material_n=TV_DIMENSION();
    for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices_n);iter.Valid();iter.Next()){
        accumulated_material_n+=volumes_n(iter.Data())*state_n(iter.Data());}

    TV_DIMENSION accumulated_material_np1=TV_DIMENSION();
    for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,int>::ITERATOR iter(indices_np1);iter.Valid();iter.Next()){
        accumulated_material_np1+=state_np1(iter.Data());}
    if((accumulated_material_np1-accumulated_material_n-hybrid_flux_sum).Magnitude_Squared() >= (T)10*std::numeric_limits<T>::epsilon()){
        std::stringstream ss;
        ss<<"ERROR: Conservation error; np1 - n - hybrid_flux = "<<accumulated_material_np1 - accumulated_material_n - hybrid_flux_sum<<std::endl
                 <<"\tTime N   = "<<accumulated_material_n<<std::endl<<"\tTime Np1 = "<<accumulated_material_np1<<std::endl<<"\tHybrid Boundary Flux = "<<hybrid_flux_sum<<std::endl;
        LOG::filecout(ss.str());}
#endif
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) template void Advect_Near_Interface_Data(const GRID<VECTOR<T,d> >&,const T,const T,const ARRAY<bool,VECTOR<int,d> >&,                                               \
                                                                           const ARRAY<bool,VECTOR<int,d> >&,const ARRAY<VECTOR<T,d+2>,FACE_INDEX<d> >&,const ARRAY<T,VECTOR<int,d> >&,               \
                                                                           HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >&,                                                      \
                                                                           const ARRAY<bool,VECTOR<int,d> >&,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >&,const ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >&, \
                                                                           const ARRAY<bool,VECTOR<int,d> >&,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >&,ARRAY<VECTOR<T,d+2>,VECTOR<int,d> >&);      \
  template void Compute_Hybrid_Boundary_Fluxes(const GRID<VECTOR<T,d> >&,const T dt,const ARRAY<bool,VECTOR<int,d> >&,const ARRAY<bool,VECTOR<int,d> >&,const ARRAY<VECTOR<T,d+2>,FACE_INDEX<d> >&,   \
                                                                           ARRAY<TRIPLE<VECTOR<int,d>,VECTOR<int,d>,VECTOR<T,d+2> > >&);
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
};
