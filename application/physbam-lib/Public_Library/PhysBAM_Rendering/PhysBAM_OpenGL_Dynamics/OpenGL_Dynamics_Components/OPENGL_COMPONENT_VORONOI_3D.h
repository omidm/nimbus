//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_VORONOI_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_VORONOI_3D__
#define __OPENGL_COMPONENT_VORONOI_3D__

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
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENT_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_POLYGON.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CHIMERA_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_VORONOI_3D: public OPENGL_COMPONENT
{
public:
    typedef VECTOR<T,3> TV;
    typedef GRID<TV> T_GRID;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    
    std::string filename;
    std::string filename_set;
    const std::string *color_map_filename;
    int frame_loaded;
    bool valid;
    bool draw_all_grids;
    int current_grid;
    VORONOI_GEOMETRY<TV> voronoi;
    ARRAY<VORONOI_GEOMETRY<TV> > voronois;
    bool is_moving_grid;
    int number_of_voronois;
    //const ARRAY<T_GRID>& grids;
    const ARRAY<OPENGL_GRID_3D<T>* >& opengl_grids;
    OPENGL_CHIMERA_SLICE<T>* chimera_slice;

    bool draw_voronoi_faces;
    bool draw_delaunay_simplices;
    bool draw_cell_status;
    int selected_simplex;
    TV_INT selected_cell_index;

    OPENGL_COMPONENT_VORONOI_3D(const std::string &input_filename,const std::string& input_filename_set,const ARRAY<OPENGL_GRID_3D<T>* >& input_opengl_grids,bool is_moving_grid_input=false):
        filename(input_filename),filename_set(input_filename_set),frame_loaded(-1),valid(false),draw_all_grids(true),current_grid(1),is_moving_grid(is_moving_grid_input),number_of_voronois(0),opengl_grids(input_opengl_grids),chimera_slice(NULL),draw_voronoi_faces(true),draw_delaunay_simplices(false),selected_simplex(0)
    {
        is_animation=true;
        while(filename_set!=""){
            std::string filename=STRING_UTILITIES::string_sprintf(filename_set.c_str(),frame,number_of_voronois+1);
            if(FILE_UTILITIES::File_Exists(filename)) number_of_voronois++;else break;}
        if(input_filename_set!="") voronois.Resize(number_of_voronois);
        Reinitialize();
    }
    ~OPENGL_COMPONENT_VORONOI_3D() {}
    
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
    void Toggle_Draw_All_Grids()
    {
        draw_all_grids=!draw_all_grids;
    }   
    void Set_Chimera_Slice(OPENGL_CHIMERA_SLICE<T>* slice_input)
    {
        chimera_slice=slice_input;
    }
    void Display(const VORONOI_GEOMETRY<TV>& v) const PHYSBAM_OVERRIDE
    {
        if(!valid || !draw)
            return;

        glEnable(GL_LIGHTING);
        
        //glLineStipple(1,0x00FF);
        //glEnable(GL_LINE_STIPPLE);
        
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glLineWidth(2);
        glEnable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
        glDisable(GL_COLOR_MATERIAL);
        //glPushName(0);
        
        //int n_simplices=0;
        //for(int j=1;j<=v.voronoi_face_simplices.Size();j++)
        //    n_simplices+=v.voronoi_face_simplices(j).y.Size();
        //LOG::cout << "voronoi display " << v.voronoi_face_simplices.Size() << " " << n_simplices << std::endl;
        
        typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX_FACE T_SIMPLEX_FACE;

        glDisable(GL_LIGHTING);
        if(draw_cell_status){
            if(voronoi.cell_indices_to_chimera_cell_status.Size()){
                for(int j=1;j<=opengl_grids.Size();j++){
                    if(voronoi.cell_indices_to_chimera_cell_status(j).Size().Product())
                        for(CELL_ITERATOR iterator(opengl_grids(j)->grid,1);iterator.Valid();iterator.Next()){
                            int status=voronoi.cell_indices_to_chimera_cell_status(j)(iterator.Cell_Index());
                            if(status){
                                glColor4f(0,1,0,1);
                                ARRAY<typename OPENGL_POLICY<T>::T_GL> valid_vertices;
                                if(status==1) glColor4f(0,1,0,1);
                                else glColor4f(0,0,1,1);
                                TV location=opengl_grids(j)->rigid_grid_frame*iterator.Location();
                                if(chimera_slice && chimera_slice->mode){
                                    bool draw_status=false;
                                    for(int i=1;i<=chimera_slice->plane1_array.Size();i++){
                                        if(!draw_all_grids && i!=j) continue;
                                        if(chimera_slice->plane1_array(i).Lazy_Inside(location) && chimera_slice->plane2_array(i).Lazy_Inside(location)){draw_status=true;break;}}
                                    if(!draw_status) continue;}
                                OpenGL_Vertex(location,valid_vertices);
                                OpenGL_Draw_Arrays(GL_POINTS,3,valid_vertices);}}}}}
        
        for(int mode=1;mode<=2;mode++)
        {
            if(mode==1)
            {
                OPENGL_MATERIAL::Matte(OPENGL_COLOR::Gray(0.5,1)).Send_To_GL_Pipeline(GL_FRONT_AND_BACK);
                glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
                glEnable(GL_LIGHTING);
                glEnable(GL_DEPTH_TEST);
            }
            else
            {
                glColor4f(1,1,1,1);
                glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
                glDisable(GL_LIGHTING);
                //glDisable(GL_DEPTH_TEST);               
            }

            if(draw_voronoi_faces)
                for(int j=1;j<=v.voronoi_faces.Size();j++)
                    if((!selected_simplex || 
                        ((v.voronoi_faces(j).x(1).x==selected_simplex && v.voronoi_faces(j).x(1).y==selected_cell_index) ||
                            (v.voronoi_faces(j).x(2).x==selected_simplex && v.voronoi_faces(j).x(2).y==selected_cell_index))) && v.voronoi_faces_geometry(j).Size()){
                        bool draw_simplex=true;
                    
                        if(!selected_simplex && chimera_slice && chimera_slice->mode){
                            TV center;
                            for(int k=1;k<=v.voronoi_faces_geometry(j).X.Size();k++)
                                center+=v.voronoi_faces_geometry(j).X(k);
                            center/=(T)v.voronoi_faces_geometry(j).X.Size();
                            draw_simplex=false;
                            for(int i=1;i<=chimera_slice->plane1_array.Size();++i){
                                if(!draw_all_grids && i!=current_grid) continue;
                                PLANE<T>& plane1=chimera_slice->plane1_array(i);
                                PLANE<T>& plane2=chimera_slice->plane2_array(i);
                                TV center_local=opengl_grids(i)->rigid_grid_frame.Inverse_Times(center);
                                if(plane1.Lazy_Inside(center_local) && plane2.Lazy_Inside(center_local)){draw_simplex=true;break;}}}
                        
                        if(draw_simplex){
                            glColor4f(0.5,0.5,0.5,1);
                            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                            ARRAY<GLfloat> normals;
                            TV normal=TRIANGLE_3D<T>::Normal(v.voronoi_faces_geometry(j).X(1),v.voronoi_faces_geometry(j).X(2),v.voronoi_faces_geometry(j).X(3));
                            for(int k=1;k<=v.voronoi_faces_geometry(j).X.Size();k++){
                                OpenGL_Vertex(v.voronoi_faces_geometry(j).X(k),vertices);
                                OpenGL_Normal(normal,normals);
                                OpenGL_Draw_Arrays_With_Normals(GL_POLYGON,3,vertices,normals);}}}
            
            if(draw_delaunay_simplices){
                for(int j=1;j<=v.boundary_delaunay_simplex_vertices.Size();j++)
                    if(!selected_simplex || j==selected_simplex){
                        VECTOR<TV,TV::dimension+1> verts=v.boundary_delaunay_simplex_vertices(j);
                        bool draw_this_tet=true;
                        if(chimera_slice && chimera_slice->mode){
                            draw_this_tet=false;
                            for(int i=1;i<=chimera_slice->plane1_array.Size();++i){
                                if(!draw_all_grids && i!=current_grid) continue;
                                PLANE<T>& plane1=chimera_slice->plane1_array(i);
                                PLANE<T>& plane2=chimera_slice->plane2_array(i);
                                TV centroid=opengl_grids(i)->rigid_grid_frame.Inverse_Times(T(0.25)*(verts(1)+verts(2)+verts(3)+verts(4)));
                                if(plane1.Lazy_Inside(centroid) && plane2.Lazy_Inside(centroid)){draw_this_tet=true;break;}
                            }
                        }
                        if(draw_this_tet){
                            glColor4f(0,1,1,1);
                            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                            ARRAY<GLfloat> normals;
                            VECTOR<TV,TV::dimension+1> norms;
                            OpenGL_Triangle(verts(1),verts(2),verts(3),vertices);
                            OpenGL_Triangle(verts(1),verts(3),verts(4),vertices);
                            OpenGL_Triangle(verts(1),verts(2),verts(4),vertices);
                            OpenGL_Triangle(verts(2),verts(3),verts(4),vertices);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(2),verts(3)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(2),verts(3)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(2),verts(3)),normals);
                        
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(3),verts(4)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(3),verts(4)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(3),verts(4)),normals);
                        
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(2),verts(4)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(2),verts(4)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(1),verts(2),verts(4)),normals);
                        
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(2),verts(3),verts(4)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(2),verts(3),verts(4)),normals);
                            OpenGL_Normal(TRIANGLE_3D<T>::Normal(verts(2),verts(3),verts(4)),normals);
                            OpenGL_Draw_Arrays_With_Normals(GL_TRIANGLES,3,vertices,normals);                    
                        }}}}

        /*glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        for(int j=1;j<=v.delaunay_mesh.elements.Size();j++)
        {
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            glColor4f(0,0,1,1);
            const VECTOR<int,4>& element=v.delaunay_mesh.elements(j);

            TETRAHEDRON<T> tetrahedron(v.cell_vertices(element(1)),v.cell_vertices(element(2)),v.cell_vertices(element(3)),v.cell_vertices(element(4)));

            OpenGL_Triangle(tetrahedron.triangle1.x1,tetrahedron.triangle1.x2,tetrahedron.triangle1.x3);
            OpenGL_Triangle(tetrahedron.triangle2.x1,tetrahedron.triangle2.x2,tetrahedron.triangle2.x3);
            OpenGL_Triangle(tetrahedron.triangle3.x1,tetrahedron.triangle3.x2,tetrahedron.triangle3.x3);
            OpenGL_Triangle(tetrahedron.triangle4.x1,tetrahedron.triangle4.x2,tetrahedron.triangle4.x3);
            
            OpenGL_Draw_Arrays(GL_TRIANGLES,3,vertices);
        }
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);*/
        
        /*glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        ARRAY<VECTOR<int,3> > fronts;
        v.tesselation_fronts.Get_Keys(fronts);
        for(int j=1;j<=fronts.Size();j++)
        {
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            glColor4f(1,0,0,1);
            //LOG::cout << "drawing front " << v.cell_vertices(fronts(j).x) << " " << v.cell_vertices(fronts(j).y) << " " << v.cell_vertices(fronts(j).z) << std::endl;
            OpenGL_Triangle(v.cell_vertices(fronts(j).x),v.cell_vertices(fronts(j).y),v.cell_vertices(fronts(j).z),vertices);
            OpenGL_Draw_Arrays(GL_TRIANGLES,3,vertices);
            vertices.Remove_All();
            TRIANGLE_3D<T> triangle(v.cell_vertices(fronts(j).x),v.cell_vertices(fronts(j).y),v.cell_vertices(fronts(j).z));
            TV center=triangle.Center();
            OpenGL_Line(center,center+triangle.Normal()*triangle.Maximum_Edge_Length()/5,vertices);
            OpenGL_Draw_Arrays(GL_LINES,3,vertices);
        }
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);*/
/*
        if(v.voronoi_vertices.Size())
        {
            glColor4f(0,1,0,1);
            ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
            OpenGL_Line(v.voronoi_vertices(1),v.voronoi_vertices(2),vertices);
            OpenGL_Draw_Arrays(GL_LINES,2,vertices);
        }

        for(int j=1;j<=v.voronoi_faces.Size();j++)
        {
          glColor4f(0,1,0,1);
          ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
          OpenGL_Line(v.voronoi_faces(j)(1),v.voronoi_faces(j)(2),vertices);
          OpenGL_Draw_Arrays(GL_LINES,2,vertices);
        }

        SEGMENT_MESH& segment_mesh=*v.delaunay_mesh.segment_mesh;
        if(v.voronoi_faces.Size() && v.mesh_face_indices_to_matrix_face_indices.Size())
        {
            for(int j=1;j<=segment_mesh.elements.Size();j++)
            {
                glColor4f(1,0,0,1);
                VECTOR<int,2> segment_indices=segment_mesh.elements(j);
                TV location=(v.cell_vertices(segment_indices(1))+v.cell_vertices(segment_indices(2)))*(T)0.5;
                TV normal=(v.cell_vertices(segment_indices(2))-v.cell_vertices(segment_indices(1))).Normalized();
                ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
                OpenGL_Line(location,location+(T)0.1*normal*v.voronoi_face_velocities(j),vertices);
                OpenGL_Draw_Arrays(GL_LINES,2,vertices);

                glColor4f(1,1,1,1);
                std::stringstream ss;
                ss << "vfi: " << j << "(" << v.mesh_face_indices_to_matrix_face_indices(j) << ")";
                OpenGL_String(location,ss.str());
            }
        }
        
        for(int j=1;j<=v.cell_vertices.Size();j++)
        {
            glColor4f(1,1,1,1);
            std::stringstream ss;
            ss << "cvi: " << j;
            if(v.mesh_cell_indices_to_matrix_cell_indices.Size())
                ss << "(" << v.mesh_cell_indices_to_matrix_cell_indices(j) << ")";;
            OpenGL_String(v.cell_vertices(j),ss.str());
            }*/
        
        //glPopName();
        glPopAttrib();
        //glDisable(GL_LINE_STIPPLE);
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
            const RANGE<TV>& domain=opengl_grids(1)->grid.Domain();
            return RANGE<VECTOR<float,3> >(VECTOR<float,3>((float)domain.min_corner.x,(float)domain.min_corner.y,0),VECTOR<float,3>((float)domain.max_corner.x,(float)domain.max_corner.y,0));
        }
        else
            return RANGE<VECTOR<float,3> >::Centered_Box();
    }
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Next_Grid,"Switch to next grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Previous_Grid,"Switch to previous grid");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Toggle_Draw_All_Grids,"Toggle to draw all the grids or not");
    
    void Toggle_Draw_Voronoi_Faces()
    {draw_voronoi_faces=!draw_voronoi_faces;selected_simplex=0;}
    void Toggle_Draw_Delaunay_Simplices()
    {draw_delaunay_simplices=!draw_delaunay_simplices;selected_simplex=0;}
    
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Toggle_Draw_Voronoi_Faces,"Toggle the drawing of Voronoi faces");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Toggle_Draw_Delaunay_Simplices,"Toggle the drawing of Delaunay simplices");
    
    void Highlight_Simplex(){
        OPENGL_WORLD::Singleton()->Prompt_User("Enter grid-cell index separated by spaces: ",Highlight_Simplex_Response_CB());}
    void Highlight_Simplex_Response(){
        if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
            std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);
            sstream>>selected_simplex>>selected_cell_index(1)>>selected_cell_index(2)>>selected_cell_index(3);}}
    
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Highlight_Simplex,"Highlight Voronoi/Delaunay simplex");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Highlight_Simplex_Response,"");

    void Toggle_Draw_Cell_Status()
    {
        draw_cell_status=!draw_cell_status;
    }
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_VORONOI_3D,Toggle_Draw_Cell_Status,"Toggle drawing of chimera cell status");

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

                //voronoi.delaunay_mesh.Get_Segment_Mesh();
                frame_loaded = frame;
                valid = true;
            }
        }
    }
};

}

#endif
