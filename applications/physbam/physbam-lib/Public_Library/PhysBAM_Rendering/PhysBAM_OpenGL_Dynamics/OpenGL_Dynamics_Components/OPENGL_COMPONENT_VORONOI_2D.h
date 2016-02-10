//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_VORONOI_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_VORONOI_2D__
#define __OPENGL_COMPONENT_VORONOI_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_2D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENT_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGLE_3D.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_VORONOI_2D: public OPENGL_COMPONENT
{
public:
    typedef VECTOR<T,2> TV;
    typedef GRID<TV> T_GRID;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    
    std::string filename;
    std::string filename_set;
    const std::string *color_map_filename;
    int frame_loaded;
    bool valid;
    VORONOI_GEOMETRY<TV> voronoi;
    ARRAY<VORONOI_GEOMETRY<TV> > voronois;
    const ARRAY<T_GRID>& grids;
    ARRAY<FRAME<TV> > frames;
    int current_grid;
    bool draw_all_grids;
    bool is_moving_grid;
    std::string rigid_grid_frame_filename_set;
    int number_of_voronois;
    bool draw_cell_status;
    bool draw_delaunay_simplices;
    int selected_simplex;

    T vector_size;

    OPENGL_COMPONENT_VORONOI_2D(const std::string &input_filename,const std::string& input_filename_set,const ARRAY<T_GRID>& input_grids,bool is_moving_grid_input=false,const std::string rigid_grid_frame_filename_set_input=""):
filename(input_filename),filename_set(input_filename_set),frame_loaded(-1),valid(false),grids(input_grids),current_grid(1),draw_all_grids(true),is_moving_grid(is_moving_grid_input),rigid_grid_frame_filename_set(rigid_grid_frame_filename_set_input),number_of_voronois(0),draw_cell_status(false),draw_delaunay_simplices(false),selected_simplex(0),vector_size(.025)
    {
        is_animation=true;
        while(filename_set!=""){
            std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,number_of_voronois+1);
            if(FILE_UTILITIES::File_Exists(filename)) number_of_voronois++;else break;}
        if(input_filename_set!="") voronois.Resize(number_of_voronois);
        Reinitialize();
    }
    ~OPENGL_COMPONENT_VORONOI_2D() {}
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE
    {
        return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
    }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE
    {
        OPENGL_COMPONENT::Set_Frame(frame_input);
        Reinitialize();
    }
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE
    {
        OPENGL_COMPONENT::Set_Draw(draw_input);
        Reinitialize();
    }

    void Toggle_Draw_All_Grids()
    {
        draw_all_grids=!draw_all_grids;
    }

    void Next_Grid()
    {
        current_grid=min(current_grid+1,number_of_voronois);
        voronoi=voronois(current_grid);
    }
    
    void Previous_Grid()
    {
        current_grid=max(current_grid-1,1);
        voronoi=voronois(current_grid);
    }

    void Display(const VORONOI_GEOMETRY<TV>& v) const{
        if(draw_cell_status){
            if(v.cell_indices_to_chimera_cell_status.Size()){
                for(int j=1;j<=grids.Size();j++){
                    if(v.cell_indices_to_chimera_cell_status(j).Size().Product())
                        for(CELL_ITERATOR iterator(grids(j),1);iterator.Valid();iterator.Next()){
                            int status=v.cell_indices_to_chimera_cell_status(j)(iterator.Cell_Index());
                            if(status){
                                glColor4f(0,1,0,1);
                                ARRAY<typename OPENGL_POLICY<T>::T_GL> valid_vertices;
                                if(status==1)
                                    glColor4f(0,1,0,1);
                                else
                                    glColor4f(0,0,1,1);
                                OpenGL_Vertex(frames(j)*iterator.Location(),valid_vertices);
                                OpenGL_Draw_Arrays(GL_POINTS,2,valid_vertices);}}}}}
            
        int n_faces=v.voronoi_faces_geometry.Size();

        for(int j=1;j<=n_faces;j++){
            glColor4f(1,0,0,1);
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            SEGMENT_2D<T> face=v.voronoi_faces_geometry(j);
            OpenGL_Line(face.x1,face.x2,vertices);
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);

            //PAIR<int,VECTOR<int,2> > i1=v.voronoi_faces(j)(1);
            //PAIR<int,VECTOR<int,2> > i2=v.voronoi_faces(j)(2);
            glColor4f(1,1,1,1);
            //OpenGL_String(face.x1*(T)0.55+face.x2*(T).45,STRING_UTILITIES::string_sprintf("%d(%d %d):%d(%d %d)",i1.x,i1.y(1),i1.y(2),i2.x,i2.y(1),i2.y(2)));

            if(v.voronoi_psi_N.Size() && v.voronoi_psi_N(j).x){
                glColor4f(0,1,0,1);            
                vertices.Remove_All();
                OpenGL_Vertex((face.x1+face.x2)*(T)0.5,vertices);
                OpenGL_Draw_Arrays(GL_POINTS,2,vertices);}}

        for(int j=1;j<=min(v.voronoi_face_velocities.Size(),n_faces);j++){
            SEGMENT_2D<T> face=v.voronoi_faces_geometry(j);
            TV location=(T).5*(face.x1+face.x2);
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            const VECTOR<GRID_CELL_INDEX,2>& indices=v.voronoi_faces(j).x;
            TV normal=(frames(indices(2).x)*grids(indices(2).x).X(indices(2).y)-frames(indices(1).x)*grids(indices(1).x).X(indices(1).y)).Normalized();
            TV end=location+vector_size*v.voronoi_face_velocities(j)*normal;
            glColor4f(0,1,0,1);
            OPENGL_SHAPES::Draw_Arrow(location,end,vertices);
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);
            glColor4f(1,1,1,1);
            OpenGL_String(end,STRING_UTILITIES::string_sprintf("%f",v.voronoi_face_velocities(j)));
        }

        for(int j=1;j<=v.vectors.Size();j++){
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            TV end=v.vectors(j)(1)+vector_size*v.vectors(j)(2);
            glColor4f(1,0,0,1);
            OPENGL_SHAPES::Draw_Arrow(v.vectors(j)(1),end,vertices);
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);
            glColor4f(1,1,1,1);
            OpenGL_String(end,STRING_UTILITIES::string_sprintf("(%f %f)",v.vectors(j)(2)(1),v.vectors(j)(2)(2)));}

        int n_grids=grids.Size();

        glColor4f(1,1,1,1);
        if(v.cell_indices_to_matrix_cell_indices.Size())
            for(int grid_index=1;grid_index<=n_grids;grid_index++){
                if(v.cell_indices_to_matrix_cell_indices(grid_index).Size().Product())
                for(CELL_ITERATOR iterator(grids(grid_index));iterator.Valid();iterator.Next()){
                    TV_INT cell_index=iterator.Cell_Index();
                    if(v.cell_indices_to_chimera_cell_status(grid_index)(cell_index)/*==2*/){
                        int mci=v.cell_indices_to_matrix_cell_indices(grid_index)(cell_index);
                        if(mci && (!selected_simplex || mci==selected_simplex))
                            OpenGL_String(frames(grid_index)*iterator.Location(),STRING_UTILITIES::string_sprintf("(%d %d):%d",cell_index(1),cell_index(2),mci));}}}

        /*if(v.face_indices_to_matrix_face_indices.Size())
          for(int grid_index=1;grid_index<=n_grids;grid_index++){
          //if(v.face_indices_to_matrix_face_indices(grid_index).Size().Product())
          for(FACE_ITERATOR iterator(grids(grid_index));iterator.Valid();iterator.Next()){
          if(v.cell_indices_to_chimera_cell_status(grid_index)(iterator.First_Cell_Index())==2 || v.cell_indices_to_chimera_cell_status(grid_index)(iterator.Second_Cell_Index())==2){
          FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
          int mfi=v.face_indices_to_matrix_face_indices(grid_index)(face_index);
          if(mfi)
          OpenGL_String(frames(grid_index)*iterator.Location(),STRING_UTILITIES::string_sprintf("(%d => %d %d):%d",face_index.axis,face_index.index(1),face_index.index(2),mfi));}}}*/

        if(v.voronoi_face_indices_to_matrix_face_indices.Size())
            for(int j=1;j<=v.voronoi_face_indices_to_matrix_face_indices.Size();j++){
                SEGMENT_2D<T> face=v.voronoi_faces_geometry(j);
                int mfi=v.voronoi_face_indices_to_matrix_face_indices(j);
                if(mfi && (!selected_simplex || mfi==selected_simplex))
                    OpenGL_String(face.x1*(T)0.45+face.x2*(T).55,STRING_UTILITIES::string_sprintf("%d:%d",j,mfi));}

        if(draw_delaunay_simplices){
            glDisable(GL_CULL_FACE);
            for(int j=1;j<=v.boundary_delaunay_simplex_vertices.Size();j++){
                if(!selected_simplex || j==selected_simplex){
                    glColor4f(0,1,1,0.2);
                    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                    VECTOR<TV,TV::dimension+1> triangle=v.boundary_delaunay_simplex_vertices(j);
                    OpenGL_Triangle(triangle(1),triangle(2),triangle(3),vertices);
                    OpenGL_Draw_Arrays(GL_TRIANGLES,2,vertices);}}
            for(int j=1;j<=v.boundary_delaunay_simplex_vertices.Size();j++){
                if(!selected_simplex || j==selected_simplex){
                    glColor4f(0,1,1,1);
                    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                    VECTOR<TV,TV::dimension+1> triangle=v.boundary_delaunay_simplex_vertices(j);
                    OpenGL_Line(triangle(1),triangle(2),vertices);
                    //OpenGL_Draw_Arrays(GL_LINES,2,vertices);
                    OpenGL_Line(triangle(3),triangle(2),vertices);
                    //OpenGL_Draw_Arrays(GL_LINES,2,vertices);
                    OpenGL_Line(triangle(1),triangle(3),vertices);
                    OpenGL_Draw_Arrays(GL_LINES,2,vertices);}}
            glEnable(GL_CULL_FACE);}
    }        

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE
    {
        if(!valid || !draw)
            return;

        glDisable(GL_LIGHTING);
        
        //glLineStipple(1,0x00FF);
        //glEnable(GL_LINE_STIPPLE);
        
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(2);
        glDisable(GL_DEPTH_TEST);
        //glPushName(0);

        if(filename_set=="" || !draw_all_grids)
            Display(voronoi);
        else for(int voronoi_index=1;voronoi_index<=number_of_voronois;voronoi_index++)
            Display(voronois(voronoi_index));
        
        //glPopName();
        glPopAttrib();
        //glDisable(GL_LINE_STIPPLE);
    }
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE
    {
        if (valid && draw)
        {
            const RANGE<TV>& domain=grids(1).Domain();
            return RANGE<VECTOR<float,3> >(VECTOR<float,3>(domain.min_corner.x,domain.min_corner.y,0),VECTOR<float,3>(domain.max_corner.x,domain.max_corner.y,0));
        }
        else
            return RANGE<VECTOR<float,3> >::Centered_Box();
    }

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Next_Grid, "Next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Previous_Grid, "Previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Toggle_Draw_All_Grids, "Toggle draw all grids");

    void Toggle_Draw_Cell_Status()
    {
        draw_cell_status=!draw_cell_status;
        selected_simplex=0;
    }
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Toggle_Draw_Cell_Status, "Toggle drawing of chimera cell status");

    void Toggle_Draw_Delaunay_Simplices()
    {
        draw_delaunay_simplices=!draw_delaunay_simplices;
        selected_simplex=0;
    }
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Toggle_Draw_Delaunay_Simplices, "Toggle drawing of delaunay triangles for interpolation");

    void Highlight_Simplex(){
        OPENGL_WORLD::Singleton()->Prompt_User("Enter Simplex Number: ",Highlight_Simplex_Response_CB());}
    void Highlight_Simplex_Response(){
        if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
            int index=0;std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);sstream>>index;
            selected_simplex=index;}}
    
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D,Highlight_Simplex,"Highlight Voronoi/Delaunay simplex");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D,Highlight_Simplex_Response,"");

    void Scale_Vector_Size(const T scale){vector_size*=scale;}
    void Increase_Vector_Size(){Scale_Vector_Size((T)1.1);}
    void Decrease_Vector_Size(){Scale_Vector_Size(1/(T)1.1);}
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_2D, Decrease_Vector_Size, "Decrease vector size");

private:
    void Reinitialize()
    {
        if (draw)
        {
            if ((is_animation && frame_loaded != frame) ||
                (!is_animation && frame_loaded < 0))
            {
                valid = false;
                if(filename_set==""){
                    std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
                    if (FILE_UTILITIES::File_Exists(tmp_filename))
                        FILE_UTILITIES::Read_From_File<RW>(tmp_filename,voronoi);
                    else
                        return;}
                else{
                    for(int i=1;i<=number_of_voronois;i++){
                        std::string tmp_filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,i);
                        if(FILE_UTILITIES::File_Exists(tmp_filename))
                            FILE_UTILITIES::Read_From_File<T>(tmp_filename,voronois(i));
                        else
                            return;}
                    voronoi=voronois(current_grid);
                }

                frames.Resize(grids.Size());
                for(int i=1;i<=grids.Size();i++)
                {
                    std::string grid_filename=STRING_UTILITIES::string_sprintf(rigid_grid_frame_filename_set.c_str(),frame,i);
                    if(FILE_UTILITIES::File_Exists(grid_filename))
                        FILE_UTILITIES::Read_From_File<T>(grid_filename,frames(i));
                }

                //voronoi.delaunay_mesh.Get_Segment_Mesh();
                frame_loaded = frame;
                valid = true;
            }
        }
    }
};

}

#endif
