//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION_CUT
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Images/EPS_FILE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/GENERALIZED_FLUID_MASS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_CUT.h>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
namespace PhysBAM{template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
FLUID_TO_SOLID_INTERPOLATION_CUT(const COLLISION_AWARE_INDEX_MAP<TV>& map,SEGMENTED_CURVE_2D<T>& curve_input,T density_input)
    :FLUID_TO_SOLID_INTERPOLATION_BASE<TV>(map),curve(curve_input),fluid_mass(0),gradient(0),density(density_input),outside_density(0),outside_fluid(0),psi_N(0),use_cut_volume(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
~FLUID_TO_SOLID_INTERPOLATION_CUT()
{
}
//#####################################################################
// Function Setup_Before_Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Setup_Before_Compute(ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_input)
{
    outside_fluid=&outside_fluid_input;
    psi_N=&psi_N_input;
    Compute_Dirichlet_Cells();
}
//#####################################################################
// Function Compute_Dirichlet_Cells
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Compute_Dirichlet_Cells()
{
    const ARRAY_VIEW<TV> X(curve.particles.X);

    CLIP_ENTRY ce;
    for(ce.i=1;ce.i<=curve.mesh.elements.m;ce.i++){
        SEGMENT_2D<T> segment(X(curve.mesh.elements(ce.i).x),X(curve.mesh.elements(ce.i).y));
        PHYSBAM_ASSERT(index_map.grid.domain.Lazy_Inside(segment.x1) && index_map.grid.domain.Lazy_Inside(segment.x2));
        RANGE<TV_INT> box(index_map.grid.Cell(segment.x1,3));
        box.Enlarge_To_Include_Point(index_map.grid.Cell(segment.x2,3));
        for(UNIFORM_ARRAY_ITERATOR<TV::m> it(box);it.Valid();it.Next())
            if(segment.Clip_To_Box(index_map.grid.Cell_Domain(it.index),ce.a,ce.b))
                cut_cells.Get_Or_Insert(it.index).clipped_segments.Append(ce);}

    Remove_Degeneracy();

    for(typename  HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        int number_cut=0;
        const ARRAY<CLIP_ENTRY>& cut_segments=it.Data().clipped_segments;
        for(int i=1;i<=cut_segments.m;i++){const CLIP_ENTRY& e=cut_segments(i);number_cut+=(e.a!=0)+(e.b!=1);}
        if(number_cut>2){
            EPS_FILE<T> eps("fail-config.eps");
            eps.Draw_Object(index_map.grid.Cell_Domain(it.Key()));
            for(int i=1;i<=cut_segments.m;i++){
                const CLIP_ENTRY& e=cut_segments(i);
                SEGMENT_2D<T> segment(X(curve.mesh.elements(e.i).x),X(curve.mesh.elements(e.i).y));
                TV e1=segment.Point_From_Barycentric_Coordinates(e.a);
                TV e2=segment.Point_From_Barycentric_Coordinates(e.b);
                eps.Line_Color(VECTOR<T,3>::Axis_Vector(i%3+1));
                eps.Draw_Line(e1,e2);}}
        PHYSBAM_ASSERT(number_cut<=2);}

    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next())
        (*outside_fluid)(it.Key())=false;

    const_cast<COLLISION_AWARE_INDEX_MAP<TV>&>(index_map).number_extra_cells=cut_cells.Size();
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Compute(const int ghost_cells)
{
    const ARRAY_VIEW<TV> X(curve.particles.X); //TV = VECTOR<T,2>, not VECTOR<T,1> or VECTOR<T,3>
    ARRAY<T> edge_length(curve.mesh.elements.m);
    ARRAY<T> particle_length(X.m);
    for(int i=1;i<=curve.mesh.elements.m;i++){
        T length=(X(curve.mesh.elements(i).x)-X(curve.mesh.elements(i).y)).Magnitude();
        edge_length(i)=length;
        particle_length(curve.mesh.elements(i).x)+=length/2;
        particle_length(curve.mesh.elements(i).y)+=length/2;}

    if(!index_map.two_phase){
        for(int i=1;i<=index_map.indexed_faces.m;i++){
            FACE_INDEX<TV::m> face=index_map.indexed_faces(i);
            if(!index_map.cell_indices(face.First_Cell_Index()) || !index_map.cell_indices(face.Second_Cell_Index())){
                unused_faces.Set(i);}}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("unused faces",0,1);}

    entries.Resize(X.m);
    int last_face=index_map.indexed_faces.m;
    int last_cell=index_map.indexed_cells.m;
    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        const TV_INT& cell=it.Key();
        CUT_CELL& cc=it.Data();
        cc.other_cell=index_map.two_phase?++last_cell:0;
        cc.face=++last_face;
        for(int i=1;i<=cc.clipped_segments.m;i++){const CLIP_ENTRY& e=cc.clipped_segments(i);
            SEGMENT_2D<T> segment(X(curve.mesh.elements(e.i).x),X(curve.mesh.elements(e.i).y));
            TV n=(segment.x2-segment.x1).Rotate_Clockwise_90()*(e.b-e.a);
            for(int a=1;a<=TV::m;a++)
                for(int k=0;k<=1;k++){
                    FACE_INDEX<TV::m> f(a,cell);
                    f.index(a)+=k;
                    Cut_Face(f,n,segment);}}
        TV N=(cc.segment.x2-cc.segment.x1).Normalized().Rotate_Clockwise_90();
        for(int i=1;i<=cc.clipped_segments.m;i++){const CLIP_ENTRY& e=cc.clipped_segments(i);
            T len=(e.b-e.a)*edge_length(e.i);
            T alpha=(e.a+e.b)/2;
            T weight_a=(1-alpha)*len/particle_length(curve.mesh.elements(e.i).x);
            T weight_b=alpha*len/particle_length(curve.mesh.elements(e.i).y);
            ENTRY ea={weight_a*N,last_face};
            ENTRY eb={weight_b*N,last_face};
            entries(curve.mesh.elements(e.i).x).Append(ea);
            entries(curve.mesh.elements(e.i).y).Append(eb);}}

    T dx=index_map.grid.dX(1);
    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        const TV_INT& cell=it.Key();
        CUT_CELL& cc=it.Data();
        VECTOR<TV,TV::m> edges;
        int count=0;
        for(int a=1;a<=TV::m;a++)
            for(int k=0;k<=1;k++){
                FACE_INDEX<TV::m> f(a,cell);
                f.index(a)+=k;
                switch(Face_Type(index_map.face_indices(f))){
                    case outside:
                    case unused:
                        edges(a)(k+1)=(T)0;
                        count++;
                        break;
                    case inside:
                        edges(a)(k+1)=dx;
                        break;
                    case cut:
                        edges(a)(k+1)=face_lengths.Get(f);
                        break;
                    default:
                        PHYSBAM_FATAL_ERROR("Should not be a coupled face!");
                        break;}}
        if(use_cut_volume)
        switch(count){
            case 0:
                cc.inside_volume=dx*dx-(T).5*abs(edges(1)(1)-edges(1)(2))*abs(edges(2)(1)-edges(2)(2));
                break;
            case 1:
            case 2:
                cc.inside_volume=(T).5*(edges(1)(1)+edges(1)(2))*(edges(2)(1)+edges(2)(2));
                break;
            default:
                PHYSBAM_FATAL_ERROR("Invalid cut cell!");
                break;}}

//    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next())
//        Add_Debug_Particle(index_map.grid.X(it.Key()),VECTOR<T,3>(0,1,0));

//    for(typename HASHTABLE<FACE_INDEX<TV::m>,T>::ITERATOR it(face_lengths);it.Valid();it.Next())
//        Add_Debug_Particle(index_map.grid.Axis_X_Face(it.Key()),VECTOR<T,3>(1,0,0));
//    PHYSBAM_DEBUG_WRITE_SUBSTEP("registered",0,1);

    Compute_Beta();
    Compute_Gradient();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("setup interp",0,1);
}
//#####################################################################
// Function Cut_Face
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Cut_Face(FACE_INDEX<TV::m>& f,const TV& normal,const SEGMENT_2D<T>& segment)
{
    int face_index=index_map.face_indices(f);
    if(!face_index || unused_faces.Contains(face_index)) return;
    if(!cut_cells.Contains(f.First_Cell_Index()) || !cut_cells.Contains(f.Second_Cell_Index())) return;

    T dx=index_map.grid.dX(3-f.axis);
    TV face_center(index_map.grid.Axis_X_Face(f)),dir=TV::Axis_Vector(3-f.axis);
    RAY<TV> ray(face_center-(T).5*dx*dir,dir);
    ray.semi_infinite=false;
    ray.t_max=dx;
    if(!INTERSECTION::Intersects(ray,segment,(T)1e-6*dx)) return;
    T length=clamp(ray.t_max,(T)1e-8*dx,dx);
    if(normal(3-f.axis)<0) length=dx-length;
    face_lengths.Set(f,length);
}
namespace{
template<class TV>
struct LOOP_ENTRY
{
    typename FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::CLIP_ENTRY e;
    VECTOR<int,TV::m> cell;
    bool operator<(const LOOP_ENTRY& h) const
    {
        if(e.i!=h.e.i) return e.i<h.e.i;
        return e.a<h.e.a;
    }
    LOOP_ENTRY * next;

    LOOP_ENTRY()
        :next(0)
    {}
};
}
//#####################################################################
// Function Remove_Degeneracy
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Remove_Degeneracy()
{
    const ARRAY_VIEW<TV> X(curve.particles.X);
    if(cut_cells.Size()<3) return;
    ARRAY<LOOP_ENTRY<TV> > list;
    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        const ARRAY<CLIP_ENTRY>& cut_segments=it.Data().clipped_segments;
        for(int i=1;i<=cut_segments.m;i++){
            const CLIP_ENTRY& e=cut_segments(i);
            LOOP_ENTRY<TV> le;
            le.cell=it.Key();
            le.e=e;
            list.Append(le);}}
    Sort(list);
    LOOP_ENTRY<TV> * pa = &list(list.m);
    for(int i=1;i<=list.m;i++) pa=pa->next=&list(i);
    for(int i=1;i<=list.m;i++,pa=pa->next){
        if(pa->cell==pa->next->cell) continue;
        LOOP_ENTRY<TV> *x=pa;
        while(x->next!=pa && x->next->cell==x->next->next->cell) x=x->next;
        if(x->next==pa) continue;
        LOOP_ENTRY<TV> *z=x->next->next;
        if(pa->cell!=z->cell) continue;
        LOG::filecout("DEGENERACY FIX: ");
        if(x->next->e.i==z->e.i){x->next=z;LOG::filecout("a");}
        if(pa->next->e.i==pa->e.i){pa->next=pa->next->next;LOG::filecout("b");}
        pa->e.b=1;
        z->e.a=0;
        for(x=pa->next;x!=z;x=x->next){x->cell=pa->cell;LOG::filecout("c");}
        LOG::filecout("\n");}

    for(LOOP_ENTRY<TV>* p=pa;;p=p->next){
        std::stringstream ss;
        while(p->next->e.b<1e-6){
            ss<<"CLIP "<<p->next->e.a<<"  "<<p->next->e.b<<std::endl;
            p->next=p->next->next;
            p->next->e.a=0;}
        if(p->next==pa) break;
        LOG::filecout(ss.str());}

    for(LOOP_ENTRY<TV>* p=pa;;p=p->next){
        std::stringstream ss;
        while(p->next->e.a>1-1e-6){
            ss<<"CLIP "<<p->next->e.a<<"  "<<p->next->e.b<<std::endl;
            p->next=p->next->next;
            p->e.b=1;}
        if(p->next==pa) break;
        LOG::filecout(ss.str());}

    for(LOOP_ENTRY<TV>* p=pa;;p=p->next){
        if(p->e.i!=p->next->e.i){
            PHYSBAM_ASSERT(p->e.b>1-1e-6 && p->next->e.a<1e-6);
            p->e.b=1;
            p->next->e.a=0;}
        else p->next->e.a=p->e.b;
        if(p->next==pa) break;}

    cut_cells.Remove_All();
    for(LOOP_ENTRY<TV>* p=pa;;p=p->next){
        cut_cells.Get_Or_Insert(p->cell).clipped_segments.Append(p->e);
        if(p->next==pa) break;}

    if(!cut_cells.Size()) return;
    PHYSBAM_ASSERT(cut_cells.Size()>=2);

    while(pa->cell==pa->next->cell) pa=pa->next;
    for(LOOP_ENTRY<TV>* p=pa->next;;p=p->next){
        LOOP_ENTRY<TV>* x=p;
        while(x->cell==x->next->cell) x=x->next;
        TV X1=X(curve.mesh.elements(p->e.i).x);
        TV X2=X(curve.mesh.elements(p->e.i).y);
        TV X3=X(curve.mesh.elements(x->e.i).x);
        TV X4=X(curve.mesh.elements(x->e.i).y);
        SEGMENT_2D<T> segment(X1+(X2-X1)*p->e.a,X3+(X4-X3)*x->e.b);
        cut_cells.Get(x->cell).segment=segment;
        p=x;
        if(p==pa) break;}
}
//#####################################################################
// Function Face_Type
//#####################################################################
template<class TV> typename FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::FACE_TYPE FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Face_Type(int f) const
{
    if(f>index_map.indexed_faces.m) return coupling;
    FACE_INDEX<TV::m> face_index=index_map.indexed_faces(f);
    if(face_lengths.Contains(face_index)) return cut;
    const CUT_CELL* a=cut_cells.Get_Pointer(face_index.First_Cell_Index());
    if(!a){
        if(!(*outside_fluid)(face_index.First_Cell_Index())) return inside;
        return index_map.two_phase?outside:unused;}
    const CUT_CELL* b=cut_cells.Get_Pointer(face_index.Second_Cell_Index());
    if(!b){
        if(!(*outside_fluid)(face_index.Second_Cell_Index())) return inside;
        return index_map.two_phase?outside:unused;}
    TV dx=index_map.grid.Axis_X_Face(face_index)-a->segment.x1;
    TV N=(a->segment.x2-a->segment.x1).Rotate_Clockwise_90();
    if(TV::Dot_Product(dx,N)<=0) return inside;
    return index_map.two_phase?outside:unused;
}
//#####################################################################
// Function Compute_Beta
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Compute_Beta()
{
    static VECTOR<T,3> colors[5]={VECTOR<T,3>(1,0,0),VECTOR<T,3>(1,1,0),VECTOR<T,3>(0,1,0),VECTOR<T,3>(0,0,1),VECTOR<T,3>(1,0,1)};
    T dx=index_map.grid.dX(1),full_volume=dx*dx,imi=Inverse(full_volume*density),imo=index_map.two_phase?Inverse(full_volume*outside_density):FLT_MAX;
    VECTOR_ND<T>& beta_inverse=(*fluid_mass)->one_over_fluid_mass_at_faces;
    T in_length=0,length=index_map.grid.Face_Size(1);
    for(int i=1;i<=index_map.indexed_faces.m;i++){
        FACE_INDEX<TV::m> face=index_map.indexed_faces(i);
        TV_INT cell1=face.index,cell2=cell1;cell1(face.axis)--;
        T volume1,volume2,inside_volume;
        if(CUT_CELL* cc=cut_cells.Get_Pointer(cell1)) volume1=cc->inside_volume;
        else if((*outside_fluid)(cell1)) volume1=(T)0;
        else volume1=full_volume;
        if(CUT_CELL* cc=cut_cells.Get_Pointer(cell2)) volume2=cc->inside_volume;
        else if((*outside_fluid)(cell2)) volume2=(T)0;
        else volume2=full_volume;
        inside_volume=(T).5*(volume1+volume2);
        //Add_Debug_Particle(index_map.grid.Axis_X_Face(index_map.indexed_faces(i)),colors[Face_Type(i)]);
        switch(Face_Type(i)){
            case coupling: break;
            case unused: beta_inverse(i)=0;break;
            case outside: 
                beta_inverse(i)=use_cut_volume?Inverse((full_volume-inside_volume)*outside_density):imo;
                break;
            case inside: 
                beta_inverse(i)=use_cut_volume?Inverse(inside_volume*density):imi;
                break;
            case cut:
                in_length=face_lengths.Get(index_map.indexed_faces(i));
                beta_inverse(i)=use_cut_volume?Inverse(inside_volume*density+(full_volume-inside_volume)*outside_density):Inverse((in_length*density+(length-in_length)*outside_density)*dx);
                break;}}

    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        //Add_Debug_Particle(it.Data().segment.Center(),VECTOR<T,3>(.5,.5,.5));
        beta_inverse(it.Data().face)=use_cut_volume?Inverse((T).5*(density*it.Data().inside_volume+outside_density*(full_volume-it.Data().inside_volume))):Inverse((T).5*((density+outside_density)*dx*it.Data().segment.Length()));}
    PHYSBAM_DEBUG_WRITE_SUBSTEP(__FUNCTION__,0,1);
}
//#####################################################################
// Function Cut_Face
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Add_Gradient_Entry(int fi,const FACE_INDEX<TV::m>& f,int side,bool outside)
{
    TV_INT cell=f.Cell_Index(side);
    int cell_index=index_map.cell_indices(cell);
    TV DX=index_map.grid.Axis_X_Face(f)+TV::Axis_Vector(f.axis)*((side==1?-1:1)*index_map.grid.dX.Min()/6);(void)DX;
    if(!cell_index){
        //Add_Debug_Particle(DX,VECTOR<T,3>(1,0,0));
        return;}
    if(outside)
        if(CUT_CELL* cc=cut_cells.Get_Pointer(cell))
            cell_index=cc->other_cell;

    T value=index_map.grid.Face_Sizes()(3-f.axis)*(side==1?-1:1);
    gradient->gradient.Append_Entry_To_Current_Row(cell_index,value);
    //Add_Debug_Particle(DX,outside?VECTOR<T,3>(0,1,0):VECTOR<T,3>(0,0,1));
}
//#####################################################################
// Function Cut_Face
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Add_Cut_Gradient_Entry(int fi,const FACE_INDEX<TV::m>& f,int side)
{
    TV_INT cell=f.Cell_Index(side);
    int cell_index=index_map.cell_indices(cell);
    TV DX=index_map.grid.Axis_X_Face(f)+TV::Axis_Vector(f.axis)*((side==1?-1:1)*index_map.grid.dX.Min()/6);(void)DX;
    if(!cell_index){
        //Add_Debug_Particle(DX,VECTOR<T,3>(1,1,0));
        return;}
    int sign=side==1?-1:1;
    T dx=index_map.grid.dX(1);
    const CUT_CELL& cc=cut_cells.Get(cell);
    T len=face_lengths.Get(f);
    gradient->gradient.Append_Entry_To_Current_Row(cell_index,sign*len);
    if(index_map.two_phase) gradient->gradient.Append_Entry_To_Current_Row(cc.other_cell,sign*(dx-len));
    //Add_Debug_Particle(DX,VECTOR<T,3>(1,0,1));
}
//#####################################################################
// Function Compute_Beta
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Compute_Gradient()
{
    gradient->gradient.Reset(index_map.Number_Cells());
    TV face_areas(index_map.grid.Face_Sizes());
    for(int i=1;i<=index_map.indexed_faces.m;i++){
        if((*psi_N)(index_map.indexed_faces(i))){
            gradient->gradient.Finish_Row();
            //Add_Debug_Particle(index_map.grid.Axis_X_Face(index_map.indexed_faces(i)),VECTOR<T,3>(1,.5,0));
            continue;}
        switch(Face_Type(i)){
            case unused: break;
            case coupling: break;
            case inside:
                Add_Gradient_Entry(i,index_map.indexed_faces(i),1,false);
                Add_Gradient_Entry(i,index_map.indexed_faces(i),2,false);
                break;
            case outside:
                Add_Gradient_Entry(i,index_map.indexed_faces(i),1,true);
                Add_Gradient_Entry(i,index_map.indexed_faces(i),2,true);
                break;
            case cut:
                Add_Cut_Gradient_Entry(i,index_map.indexed_faces(i),1);
                Add_Cut_Gradient_Entry(i,index_map.indexed_faces(i),2);
                break;}
        gradient->gradient.Finish_Row();}

    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        const CUT_CELL& cc=it.Data();
        T len=cc.segment.Length();
        //Add_Debug_Particle(cc.segment.Center(),VECTOR<T,3>(.5,.5,.5));
        gradient->gradient.Append_Entry_To_Current_Row(index_map.cell_indices(it.Key()),-len);
        if(index_map.two_phase) gradient->gradient.Append_Entry_To_Current_Row(cc.other_cell,len);
        gradient->gradient.Finish_Row();}
    gradient->gradient.Sort_Entries();
    PHYSBAM_DEBUG_WRITE_SUBSTEP(__FUNCTION__,0,1);
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Times_Add(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const
{
    for(int i=1;i<=entries.m;i++){const ARRAY<ENTRY>& array=entries(i);
        for(int j=1;j<=array.m;j++){const ENTRY& e=array(j);
            solid_velocity.V.array(i)+=e.w*fluid_velocity(e.i);}}
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,VECTOR_ND<T>& fluid_force) const
{
    for(int i=1;i<=entries.m;i++){const ARRAY<ENTRY>& array=entries(i);
        for(int j=1;j<=array.m;j++){const ENTRY& e=array(j);
            fluid_force(e.i)+=TV::Dot_Product(e.w,solid_force.V.array(i));}}
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("H-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("H",G.Raw_Size(),fluid_faces);
    ARRAY<int> reverse_map_deformable(G.V.array.Size());
    reverse_map_deformable.Subset(G.V.indices)=IDENTITY_ARRAY<>(G.V.Size());

    for(int i=1;i<=entries.m;i++){const ARRAY<ENTRY>& array=entries(i);
        for(int a=1;a<=TV::m;a++){
            HASHTABLE<int,T> face_weights;
            for(int j=1;j<=array.m;j++){const ENTRY& e=array(j);
                if(face_weights.Contains(e.i)) face_weights.Get(e.i)+=e.w(a);
                else face_weights.Insert(e.i,e.w(a));}
            for(typename HASHTABLE<int,T>::CONST_ITERATOR it(face_weights);it.Valid();it.Next())
                oo.Add_Sparse_Entry((reverse_map_deformable(i)-1)*TV::m+a,it.Key(),it.Data());}}
    oo.End_Sparse_Matrix();

    ARRAY<VECTOR<int,TV::m+1> > extra_faces;
    for(typename HASHTABLE<TV_INT,CUT_CELL>::CONST_ITERATOR it(cut_cells);it.Valid();it.Next())
        extra_faces.Append(it.Key().Append(it.Data().face));
    OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("extra-map-%i.txt",n).c_str()).Write("extra_map",extra_faces);
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    ARRAY<int> reverse_map_deformable(this->V_size);
    reverse_map_deformable.Subset(*this->V_indices)=IDENTITY_ARRAY<>(this->V_size);

    for(int i=1;i<=entries.m;i++){const ARRAY<ENTRY>& array=entries(i);
        for(int a=1;a<=TV::m;a++){
            HASHTABLE<int,T> face_weights;
            for(int j=1;j<=array.m;j++){const ENTRY& e=array(j);
                if(face_weights.Contains(e.i)) face_weights.Get(e.i)+=e.w(a);
                else face_weights.Insert(e.i,e.w(a));}
            for(typename HASHTABLE<int,T>::CONST_ITERATOR it(face_weights);it.Valid();it.Next())
                data.Append(TRIPLE<int,int,T>((reverse_map_deformable(i)-1)*TV::m+a,it.Key(),it.Data()));}}
}
//#####################################################################
// Function Fill_Extra_Velocities
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Fill_Extra_Velocities(VECTOR_ND<T>& fluid_velocity_vector) const
{
    VECTOR_ND<T> div(gradient->gradient.m);
    for(int i=index_map.indexed_faces.m+1;i<=gradient->gradient.m;i++) fluid_velocity_vector(i)=0;
    gradient->Transpose_Times(fluid_velocity_vector,div);
    for(int i=index_map.indexed_faces.m+1;i<=gradient->gradient.m;i++){
        int o=gradient->gradient.offsets(i);
        PHYSBAM_ASSERT(gradient->gradient.offsets(i+1)==o+1+index_map.two_phase);
        T v=-div(gradient->gradient.A(o).j)/gradient->gradient.A(o).a;
        if(index_map.two_phase) v=(v-div(gradient->gradient.A(o+1).j)/gradient->gradient.A(o+1).a)/2;
        fluid_velocity_vector(i)=v;}
}
//#####################################################################
// Function Dump_Extra_Velocities
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION_CUT<TV>::
Dump_Extra_Velocities(const VECTOR_ND<T>& fluid_velocity_vector)
{
/*
    T mn=FLT_MAX,mx=-mn;
    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        T x=fluid_velocity_vector(it.Data().face)/(it.Data().segment.x2-it.Data().segment.x1).Magnitude();
        if(x<mn) mn=x;
        if(x>mx) mx=x;}
    LOG::cout<<"RANGE: "<<mn<<" "<<mx<<std::endl;
    if(mx==mn){mx++;mn--;}

    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(mn,mx,false,false,false);

    for(typename HASHTABLE<TV_INT,CUT_CELL>::ITERATOR it(cut_cells);it.Valid();it.Next()){
        VECTOR<T,2> N=(it.Data().segment.x2-it.Data().segment.x1).Rotate_Clockwise_90();
        VECTOR<T,2> V=fluid_velocity_vector(it.Data().face)*N.Normalized();
        T x=fluid_velocity_vector(it.Data().face);
        Add_Debug_Particle(it.Data().segment.Center(),color_map(x));
        Debug_Particle_Set_Attribute<VECTOR<T,2>,VECTOR<T,2> >(ATTRIBUTE_ID_V,V);}*/
}
//#####################################################################
template class FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<float,2> >;
/*template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<float,1> >::Dump_Extra_Velocities(const VECTOR_ND<float>&);
template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<float,1> >::Fill_Extra_Velocities(VECTOR_ND<float>&) const;
template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<float,3> >::Dump_Extra_Velocities(const VECTOR_ND<float>&);
template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<float,3> >::Fill_Extra_Velocities(VECTOR_ND<float>&) const;*/
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<double,2> >;
/*template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<double,1> >::Dump_Extra_Velocities(const VECTOR_ND<double>&);
template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<double,1> >::Fill_Extra_Velocities(VECTOR_ND<double>&) const;
template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<double,3> >::Dump_Extra_Velocities(const VECTOR_ND<double>&);
template void FLUID_TO_SOLID_INTERPOLATION_CUT<VECTOR<double,3> >::Fill_Extra_Velocities(VECTOR_ND<double>&) const;*/
#endif
