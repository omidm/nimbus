//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Level_Sets/LEVELSET_RED_GREEN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Level_Sets/READ_WRITE_LEVELSET_RED_GREEN.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_ADAPTIVE_NODE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D(const std::string& prefix,const int start_frame)
    :OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T,RW>(prefix,start_frame),
    deformable_body_collection(0,collision_body_list),has_embedded_objects(false)
{
    deformable_geometry_collection=&deformable_body_collection.deformable_geometry;
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D()
{
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Reinitialize(bool force)
{
    if(!(draw && (force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)))) return;
    static bool first_time=true;
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/deformable_object_structures",prefix.c_str(),frame);
    bool read_static_variables=!deformable_body_collection.deformable_geometry.structures.m;
    int static_frame=-1;
    if(FILE_UTILITIES::File_Exists(filename)){static_frame=frame;read_static_variables=true;}
    if(read_static_variables && first_time) LOG::filecout("Deformable bodies static variables will be read each frame\n");
    deformable_body_collection.Read(STREAM_TYPE(RW()),prefix,prefix,frame,static_frame,read_static_variables,true); // Currently this will exit if any of the files don't exist... we should
                                                                                                                     // change it to not do that

    BASE::Reinitialize(force);

    if(read_static_variables){
        int m=deformable_body_collection.deformable_geometry.structures.m;
        embedded_curve_objects.Delete_Pointers_And_Clean_Memory();embedded_curve_objects.Resize(m);
        segmented_curve_objects.Delete_Pointers_And_Clean_Memory();segmented_curve_objects.Resize(m);
        triangulated_area_objects.Delete_Pointers_And_Clean_Memory();triangulated_area_objects.Resize(m);
        triangles_of_material_objects.Delete_Pointers_And_Clean_Memory();triangles_of_material_objects.Resize(m);
        free_particles_objects.Delete_Pointers_And_Clean_Memory();free_particles_objects.Resize(m);
        free_particles_indirect_arrays.Delete_Pointers_And_Clean_Memory();free_particles_indirect_arrays.Resize(m);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        grid_list.Delete_Pointers_And_Clean_Memory();grid_list.Resize(m);
#endif
        phi_list.Delete_Pointers_And_Clean_Memory();phi_list.Resize(m);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        phi_objects.Delete_Pointers_And_Clean_Memory();phi_objects.Resize(m);
#endif
        int color_map_index=15;
        for(int i=1;i<=m;i++){
            STRUCTURE<TV>* structure=deformable_body_collection.deformable_geometry.structures(i);
            if(EMBEDDED_MATERIAL_SURFACE<TV,2>* embedding=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,2>*>(structure)){
                has_embedded_objects=true;
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": embedded triangulated area\n";LOG::filecout(ss.str());}
                triangles_of_material_objects(i)=new OPENGL_TRIANGULATED_AREA<T>(embedding->material_surface,false,OPENGL_COLOR::Red(),OPENGL_COLOR::Blue());
                embedding->embedded_object.simplicial_object.mesh.Initialize_Segment_Mesh();
                triangulated_area_objects(i)=new OPENGL_TRIANGULATED_AREA<T>(embedding->embedded_object.simplicial_object,true);
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": embedded segmented curve\n";LOG::filecout(ss.str());}
                embedded_curve_objects(i)=new OPENGL_SEGMENTED_CURVE_2D<T>(embedding->embedded_object.embedded_object);
                embedded_curve_objects(i)->draw_velocities=draw_velocities;
                embedded_curve_objects(i)->velocity_scale=velocity_scale;} // apply current parameters
            else if(SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(structure)){
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": segmented curve\n";LOG::filecout(ss.str());}
                segmented_curve_objects(i)=new OPENGL_SEGMENTED_CURVE_2D<T>(*segmented_curve,color_map->Lookup(color_map_index--));
                segmented_curve_objects(i)->draw_velocities=draw_velocities;
                segmented_curve_objects(i)->velocity_scale=velocity_scale;} // apply current parameters
            else if(TRIANGULATED_AREA<T>* triangulated_area=dynamic_cast<TRIANGULATED_AREA<T>*>(structure)){
                if(first_time) {std::stringstream ss;ss<<"object "<<i<<": triangulated area\n";LOG::filecout(ss.str());}
                triangulated_area->mesh.Initialize_Segment_Mesh(); // to enable segment selection
                triangulated_area_objects(i)=new OPENGL_TRIANGULATED_AREA<T>(*triangulated_area,true);}
            else PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Weird object %d",i));}}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    for(int i=1;i<=deformable_geometry_collection->structures.m;i++){
        std::string suffix=STRING_UTILITIES::string_sprintf("_%d",i);
        std::string frame_prefix=STRING_UTILITIES::string_sprintf("%s/%d",prefix.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(frame_prefix+"levelset_red_green"+suffix)){
            if(first_time) {std::stringstream ss;ss<<"adding red green levelset "<<i<<"\n";LOG::filecout(ss.str());}
            grid_list(i)=new RED_GREEN_GRID_2D<T>;phi_list(i)=new ARRAY<T>;
            LEVELSET_RED_GREEN<TV> levelset(*grid_list(i),*phi_list(i));
            FILE_UTILITIES::Read_From_File<RW>(frame_prefix+"levelset_red_green"+suffix,levelset);
            static OPENGL_LEVELSET_COLOR_MAP<T> color_map(OPENGL_COLOR::Red(),OPENGL_COLOR::Blue());
            phi_objects(i)=new OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<T> >(*grid_list(i),*phi_list(i),&color_map);}}
#endif
    first_time=false;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Display(const int in_color) const
{
    BASE::Display(in_color);

    if(!draw||!valid) return;
    bool draw_embedded_curves=false;
    switch(display_mode){
        case 0:draw_embedded_curves=true;break;}
    for(int i=1;i<=segmented_curve_objects.m;i++){
        glPushName(i);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
        if(phi_objects(i)){glPushName(4);phi_objects(i)->Display(in_color);glPopName();}
#endif
        if(draw_embedded_curves && embedded_curve_objects(i)){glPushName(5);embedded_curve_objects(i)->Display(in_color);glPopName();}
        glPopName();}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=BASE::Bounding_Box();
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(draw && valid && deformable_geometry_collection->structures.m>0){
        for(int i=1;i<=phi_objects.m;i++) if(phi_objects(i)) box.Enlarge_To_Include_Box(phi_objects(i)->Bounding_Box());}
#endif
    return box;
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Set_Vector_Size(const T vector_size)
{
    BASE::Set_Vector_Size(vector_size);
    for(int i=1;i<=embedded_curve_objects.m;i++) if(embedded_curve_objects(i)) embedded_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Increase_Vector_Size()
{
    BASE::Increase_Vector_Size();
    for(int i=1;i<=embedded_curve_objects.m;i++) if(embedded_curve_objects(i)) embedded_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Decrease_Vector_Size()
{
    BASE::Decrease_Vector_Size();
    for(int i=1;i<=embedded_curve_objects.m;i++) if(embedded_curve_objects(i)) embedded_curve_objects(i)->velocity_scale=velocity_scale; // apply to objects which already exist
}
//#####################################################################
// Function Toggle_Draw_Velocities
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Toggle_Draw_Velocities()
{
    BASE::Toggle_Draw_Velocities();
    for(int i=1;i<=embedded_curve_objects.m;i++) if(embedded_curve_objects(i)) embedded_curve_objects(i)->draw_velocities=draw_velocities;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION* OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_2D<T>* selection=0;
    if(buffer_size >= 1){
        // TODO: resolve the potential conflict of selection numbers with the base class
        if(buffer[1]!=4 && buffer[1]!=5) return BASE::Get_Selection(buffer,buffer_size);

        selection=new OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_2D<T>(this);
        selection->body_index=buffer[0];selection->subobject_type=buffer[1];
        switch(selection->subobject_type){
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
            case 4:selection->subobject=phi_objects(buffer[0]);break;
#endif
            case 5:selection->subobject=embedded_curve_objects(buffer[0]);break;
            default:delete selection;return 0;}
        selection->body_selection=selection->subobject->Get_Selection(&buffer[2],buffer_size-2);
        if(!selection->body_selection){delete selection;selection=0;}}
    return selection;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T,RW>::
Clear_Highlight()
{
    BASE::Clear_Highlight();
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    for(int i=1;i<=phi_objects.m;i++)if(phi_objects(i))phi_objects(i)->Clear_Highlight();
#endif
    for(int i=1;i<=embedded_curve_objects.m;i++)if(embedded_curve_objects(i))embedded_curve_objects(i)->Clear_Highlight();
}
//#####################################################################
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<double,double>;
#endif
