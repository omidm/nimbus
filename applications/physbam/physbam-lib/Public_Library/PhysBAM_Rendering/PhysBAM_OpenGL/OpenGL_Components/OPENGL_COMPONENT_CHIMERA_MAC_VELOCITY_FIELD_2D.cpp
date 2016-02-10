//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D.h>
#include <PhysBAM_Tools/Math_Tools/Is_NaN.h>
using namespace PhysBAM;
    
template<class T,class RW> OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D(const ARRAY<GRID<TV> > &grid_input,const std::string &velocity_filename_set_input,const std::string filename_active_cells_input,const std::string filename_active_faces_input,bool is_moving_grid_input,std::string rigid_grid_frame_filename_set_input)
    : OPENGL_COMPONENT("MAC Velocity Field 2D"),draw_vorticity(false),
     velocity_filename_set(velocity_filename_set_input),filename_active_cells(filename_active_cells_input),filename_active_faces(filename_active_faces_input),current_grid(1),valid(false),draw_divergence(false),draw_all_grids(true),draw_streamlines(false),use_seed_for_streamlines(false),opengl_divergence_field(0),streamlines(*new SEGMENT_MESH(),*new GEOMETRY_PARTICLES<TV>()),opengl_streamlines(streamlines),psi_N_psi_D_basedir(""),min_vorticity(-1),max_vorticity(1),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input),vorticity_mode(0)
{
    Initialize(grid_input);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Initialize(const ARRAY<GRID<TV> > &grid_input)
{
    is_animation=true;
    frame_loaded=-1;

    int number_of_grids=grid_input.Size();

    opengl_mac_velocity_fields.Resize(number_of_grids);
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) 
        opengl_mac_velocity_fields(i)=new OPENGL_MAC_VELOCITY_FIELD_2D<T>(*(new GRID<TV>(grid_input(i))),*(new ARRAY<T,FACE_INDEX<2> >),0,0,is_moving_grid);
    opengl_mac_velocity_field=opengl_mac_velocity_fields(1);
    number_of_steps=2*opengl_mac_velocity_field->grid.counts.x;
    opengl_vorticity_magnitude=new OPENGL_SCALAR_FIELD_2D<T>(*(new T_GRID),*(new ARRAY<T,VECTOR<int,2> >),OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1),OPENGL_SCALAR_FIELD_2D<T>::DRAW_TEXTURE,is_moving_grid);
    opengl_vorticity_magnitude->Set_Scale_Range(min_vorticity,max_vorticity);

    OPENGL_COLOR_RAMP<T>* ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color((T)-1e+2,OPENGL_COLOR::Red());
    ramp->Add_Color(-1,OPENGL_COLOR::Yellow());
    ramp->Add_Color((T)-1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(0,OPENGL_COLOR::Black());
    ramp->Add_Color((T)1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(1,OPENGL_COLOR::Yellow());
    ramp->Add_Color((T)1e+2,OPENGL_COLOR::Red());
    opengl_divergence_field=new OPENGL_SCALAR_FIELD_2D<T>(*(new T_GRID),divergence,OPENGL_COLOR_RAMP<T>::Matlab_Jet(-1.5,1.5),OPENGL_SCALAR_FIELD_2D<T>::DRAW_TEXTURE,is_moving_grid);

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
~OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++){
        delete &opengl_mac_velocity_fields(i)->grid;
        delete &opengl_mac_velocity_fields(i)->u;
        delete &opengl_mac_velocity_fields(i)->v;}
    delete opengl_divergence_field;
}

template<class T,class RW> bool OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf(velocity_filename_set.c_str(),frame_input,1));
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_mac_velocity_field->Print_Selection_Info(stream,selection);
        if(draw_vorticity && ((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)->index;
            stream<<"vorticity (" << vorticity_mode << ")  magnitude = "<<opengl_vorticity_magnitude->values(index)<<std::endl;}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Display(const int in_color) const
{
    if(valid){
        if(draw){
            if(draw_all_grids) for(int i=1;i<=opengl_mac_velocity_fields.m;i++)
                opengl_mac_velocity_fields(i)->Display(in_color);
            else
                opengl_mac_velocity_fields(current_grid)->Display(in_color);
            if(draw_divergence) opengl_divergence_field->Display(in_color);
            if(draw_vorticity) opengl_vorticity_magnitude->Display(in_color);
            if(draw_streamlines) opengl_streamlines.Display(in_color);
        }
    }
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_mac_velocity_field->Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Reinitialize()
{
    if (draw || draw_divergence){
        if ((is_animation && (frame_loaded!=frame)) || (!is_animation && frame_loaded < 0)){
            valid = false;
            if(is_moving_grid){
                for(int i=1;i<=opengl_mac_velocity_fields.m;i++){
                    std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename_set.c_str(),frame,i);
                    if(FILE_UTILITIES::File_Exists(grid_filename)) FILE_UTILITIES::Read_From_File<T>(grid_filename,opengl_mac_velocity_fields(i)->rigid_grid_frame);}}
            for(int i=1;i<=opengl_mac_velocity_fields.m;i++){
                std::string tmp_filename=STRING_UTILITIES::string_sprintf(velocity_filename_set.c_str(),frame,i);
                if(FILE_UTILITIES::File_Exists(tmp_filename)){
                    FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_mac_velocity_fields(i)->face_velocities);}
                else return;
                opengl_mac_velocity_fields(i)->Update();}
            frame_loaded=frame;
            valid=true;
            Update_Divergence();
            Update_Streamlines();
            if(draw_vorticity){Update_Vorticity();opengl_vorticity_magnitude->Update();}}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Update_Divergence()
{
    if(draw_divergence && valid){
        opengl_divergence_field->rigid_grid_frame=opengl_mac_velocity_field->rigid_grid_frame;
        opengl_divergence_field->grid=opengl_mac_velocity_field->grid;
        GRID<TV>& grid=opengl_mac_velocity_field->grid;
        ARRAY_VIEW<T,VECTOR<int,2> > &u=opengl_mac_velocity_field->u,&v=opengl_mac_velocity_field->v;
        static T_FACE_ARRAYS_BOOL psi_N;
        static T_ARRAYS_BOOL psi_D;
        bool got_all_psi=true;
        if(!psi_N_psi_D_basedir.empty()){
            std::string psi_N_filename=STRING_UTILITIES::string_sprintf("%s/%d/psi_N_%d",psi_N_psi_D_basedir.c_str(),frame,current_grid);
            std::string psi_D_filename=STRING_UTILITIES::string_sprintf("%s/%d/psi_D_%d",psi_N_psi_D_basedir.c_str(),frame,current_grid);
            if(FILE_UTILITIES::File_Exists(psi_N_filename)) FILE_UTILITIES::Read_From_File<RW>(psi_N_filename,psi_N);
            else got_all_psi=false;
            if(FILE_UTILITIES::File_Exists(psi_D_filename)) FILE_UTILITIES::Read_From_File<RW>(psi_D_filename,psi_D);
            else got_all_psi=false;}
        else got_all_psi=false;
        if(!got_all_psi){psi_N.Clean_Memory();psi_D.Clean_Memory();}
        divergence.Resize(1,grid.counts.x,1,grid.counts.y);
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
            if(got_all_psi && (psi_D(i,j) || (psi_N.Component(1)(i,j) && psi_N.Component(1)(i+1,j) && psi_N.Component(2)(i,j) && psi_N.Component(2)(i,j+1)))) divergence(i,j)=0;
            else divergence(i,j)=grid.one_over_dX.x*(u(i+1,j)-u(i,j))+grid.one_over_dX.y*(v(i,j+1)-v(i,j));}
        opengl_divergence_field->Update();}
    else divergence.Clean_Memory();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Update_Streamlines()
{
    streamlines.Clean_Memory();streamlines.particles.array_collection->Clean_Memory();streamlines.mesh.Clean_Memory();
    if(!draw_streamlines || !valid) return;
    
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    int number_of_streamlines=100;
    T step_length=(T).5*grid.Minimum_Edge_Length();

    RANDOM_NUMBERS<T> random;
    if(use_seed_for_streamlines) random.Set_Seed(streamline_seed);
    T_LINEAR_INTERPOLATION_VECTOR linear_interpolation;
    T_FACE_ARRAYS_SCALAR mac_velocity_field(grid);
    mac_velocity_field.Component(1)=opengl_mac_velocity_field->u;
    mac_velocity_field.Component(2)=opengl_mac_velocity_field->v;
    FACE_LOOKUP_UNIFORM<GRID<TV> > V_lookup(mac_velocity_field);

    for(int i=1;i<=number_of_streamlines;i++){
        int p=streamlines.particles.array_collection->Add_Element();
        TV X=random.Get_Uniform_Vector(grid.domain);
        streamlines.particles.X(p)=opengl_mac_velocity_field->rigid_grid_frame*X;
        for(int step=1;step<=number_of_steps;step++){
            TV velocity=linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X);
            TV X_new=X+step_length*velocity;
            velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X_new));
            X_new=grid.Clamp(X+step_length*velocity);
            int new_particle=streamlines.particles.array_collection->Add_Element();
            streamlines.particles.X(new_particle)=opengl_mac_velocity_field->rigid_grid_frame*X_new;
            streamlines.mesh.elements.Append(VECTOR<int,2>(p,new_particle));
            p=new_particle;
            X=X_new;}}
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Velocity_Mode()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Toggle_Velocity_Mode();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Velocity_Mode_And_Draw()
{
    if (draw)
    {
        Toggle_Velocity_Mode();
        if ((int)opengl_mac_velocity_field->velocity_mode==0) Toggle_Draw();
    }
    else Toggle_Draw();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Increase_Vector_Size()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Scale_Vector_Size((T)1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Decrease_Vector_Size()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->Scale_Vector_Size(1/(T)1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Arrowhead()
{
    for(int i=1;i<=opengl_mac_velocity_fields.m;i++) opengl_mac_velocity_fields(i)->draw_arrowhead = !opengl_mac_velocity_fields(i)->draw_arrowhead;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_Divergence()
{
    draw_divergence=!draw_divergence;
    Update_Divergence();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_All_Grids()
{
    draw_all_grids=!draw_all_grids;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_Streamlines()
{
    draw_streamlines=!draw_streamlines;
    Update_Streamlines();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Use_Streamline_Seed()
{
    use_seed_for_streamlines=!use_seed_for_streamlines;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Streamline_Seed(const unsigned int seed)
{
    streamline_seed=seed;
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Lengthen_Streamlines()
{
    number_of_steps+=10;
    Update_Streamlines();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Shorten_Streamlines()
{
    number_of_steps=max(number_of_steps-10,0);
    Update_Streamlines();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Next_Grid(){
    current_grid=min(current_grid+1,opengl_mac_velocity_fields.m);
    opengl_mac_velocity_field=opengl_mac_velocity_fields(current_grid);
    //Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Previous_Grid(){
    current_grid=max(current_grid-1,1);
    opengl_mac_velocity_field=opengl_mac_velocity_fields(current_grid);
    //Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Draw_Vorticity()
{
    vorticity_mode=(vorticity_mode+1)%6;
    draw_vorticity=vorticity_mode!=0;
    if(draw_vorticity) valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Normalize_Vorticity_Color_Map()
{
    if(!draw_vorticity) return;
    LOG::cout << "vorticity range " << min_vorticity << " " << max_vorticity << std::endl;
    opengl_vorticity_magnitude->Set_Scale_Range(min_vorticity,max_vorticity);
    valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Increase_Vorticity_Range()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude->Increase_Scale_Range();
    valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Decrease_Vorticity_Range()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude->Decrease_Scale_Range();
    valid=false;
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<T,RW>::
Update_Vorticity()
{
    typedef VECTOR<int,2> TV_INT;

    int ghost_cells=2;
    opengl_vorticity_magnitude->rigid_grid_frame=opengl_mac_velocity_field->rigid_grid_frame;
    opengl_vorticity_magnitude->grid=opengl_mac_velocity_field->grid;
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    ARRAY_VIEW<T,VECTOR<int,2> > &u=opengl_mac_velocity_field->u,&v=opengl_mac_velocity_field->v;
    FACE_LOOKUP_UNIFORM<GRID<TV> > lookup(opengl_mac_velocity_field->face_velocities);
    TV dx=grid.DX();
    opengl_vorticity_magnitude->values.Resize(grid.Domain_Indices(ghost_cells));
    max_vorticity=0;
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,grid.Domain_Indices(ghost_cells));iterator.Valid();iterator.Next()){
        VECTOR<int,2> index=iterator.Cell_Index();
        T value=0;
        switch(vorticity_mode){
            case 1:
                value=VORTICITY_UNIFORM<TV>::Vorticity(grid,lookup,index)(1);
                break;
            case 2:
                value=(u(index+TV_INT(1,0))-u(index))/dx(1);
                break;
            case 3:
                value=(u(index+TV_INT(0,1))-u(index))/dx(2);
                break;
            case 4:
                value=(v(index+TV_INT(1,0))-v(index))/dx(1);
                break;
            case 5:
                value=(v(index+TV_INT(0,1))-v(index))/dx(2);
                break;
        }
        opengl_vorticity_magnitude->values(index)=value;
        if(!Is_NaN(value))
            max_vorticity=max(max_vorticity,(T)fabs(value));}
    min_vorticity=-max_vorticity;
}
template class OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CHIMERA_MAC_VELOCITY_FIELD_2D<double,double>;
#endif
