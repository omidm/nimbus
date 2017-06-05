//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/LOCAL_GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Parallel_Computation/READ_WRITE_LOCAL_GRID.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/BUBBLE_PARTICLES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Read_Write/Particles/READ_WRITE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Class MERGER
//#####################################################################
template<class T,class T_GRID,class RW>
class MERGER
{
public:
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS; 
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_ARRAYS_SCALAR::template REBIND<VECTOR<T,T_GRID::dimension+2> >::TYPE T_ARRAYS_DIMENSION_SCALAR;
    

    const int number_of_processes;
    int number_of_fluid_processes;
    int fluid_proc_offset;
    const std::string input_directory,output_directory;
    T_GRID grid;
    ARRAY<LOCAL_GRID<T_GRID>*> local_grids;

    ARRAY<int> needs_init,needs_destroy;
    SOLID_BODY_COLLECTION<TV>* solid_body_collection;
    RIGID_BODY_COLLECTION<TV> rigid_body_collection;
    ARRAY<RIGID_BODY_COLLECTION<TV>*> local_rigid_collections;
    ARRAY<int> colors;

    bool merge_levelset;
    bool merge_object_levelset;
    bool merge_debug_data;
    bool merge_particles;
    bool merge_removed_particles;
    bool merge_removed_particle_times;
    bool merge_velocities;
    bool merge_temperatures;
    bool merge_densities;
    bool merge_compressible;
    bool merge_solid_fluid;
    bool print_max_errors;
    bool merge_rigids;

    MERGER(const int number_of_processes_input,const std::string& input_directory_input,const std::string& output_directory_input,const bool print_max_errors,const bool merge_solid_fluid_input,const bool merge_rigid_input)
        :number_of_processes(number_of_processes_input),number_of_fluid_processes(number_of_processes_input),fluid_proc_offset(0),input_directory(input_directory_input),
        output_directory(output_directory_input),rigid_body_collection(0,0),merge_levelset(true),merge_object_levelset(true),merge_debug_data(true),merge_particles(true),merge_removed_particles(true),
        merge_removed_particle_times(false),merge_velocities(true),merge_temperatures(true),merge_densities(true),merge_compressible(false),merge_solid_fluid(merge_solid_fluid_input),
        print_max_errors(print_max_errors),merge_rigids(merge_rigid_input)
    {
        if(merge_solid_fluid){
            number_of_fluid_processes=number_of_processes-1;
            fluid_proc_offset=number_of_processes-number_of_fluid_processes;
            FILE_UTILITIES::Read_From_File<RW>(input_directory+"/2/common/global_grid",grid);}
        else FILE_UTILITIES::Read_From_File<RW>(input_directory+"/1/common/global_grid",grid);
        grid=grid.Get_MAC_Grid();
    }

    ~MERGER()
    {delete solid_body_collection;}

    bool Need_Merge(const std::string &filename)
    {std::string output_filename=output_directory+"/"+filename;
    if(!FILE_UTILITIES::File_Exists(output_filename)) return true;
    // check file times
    std::string prefix=input_directory+"/"+filename+".";
    for(int i=1;i<=number_of_processes;i++)if(FILE_UTILITIES::Compare_File_Times(input_directory+STRING_UTILITIES::string_sprintf("/%d/",i)+filename,output_filename)>0) return true;
    return false;}

    bool Source_Files_Exist(const int frame) const
    {for(int i=1;i<=number_of_processes;i++)if(!FILE_UTILITIES::File_Exists(input_directory+"/"+STRING_UTILITIES::string_sprintf("%d/%d",i,frame)+"/time")) return false;
    return true;}

    void Merge_All_Frames(const int first_frame,const int last_frame)
    {
        if(merge_solid_fluid) solid_body_collection=new SOLID_BODY_COLLECTION<TV>(0,0);
        for(int frame=first_frame;frame<=last_frame;frame++){
            //if(!Source_Files_Exist(frame)){LOG::cout<<"missing source files for frame "<<frame<<std::endl;break;}
        Merge(frame);}
    }

//#####################################################################
    bool Merge_Rigid_Data(const int frame);
    void Merge(const int frame);
    template<class T_ARRAYS> bool Merge_Cell_Data(const std::string& filename,const int verify_bandwidth,const bool scale=false);
    template<class T_ARRAYS> bool Merge_Levelset(const std::string& filename,const int verify_bandwidth);
    template<class T_FACE_ARRAYS_2> bool Merge_Face_Data(const std::string& filename,const int verify_bandwidth);
    template<class T_LIST_2> bool Merge_Lists(const std::string& filename);
    template<class T_PARTICLES> bool Merge_Particles(const std::string& filename);
    template<class T_PARTICLES> bool Merge_Cell_Particles(const std::string& filename);
    template<class T_ARRAYS> void Scale_Cell_Data(T_ARRAYS& array){PHYSBAM_NOT_IMPLEMENTED();}
    void Scale_Cell_Data(T_ARRAYS_SCALAR& array);
//#####################################################################
};
//#####################################################################
// Function Merge
//#####################################################################
template<class T,class T_GRID,class RW> void MERGER<T,T_GRID,RW>::
Merge(const int frame)
{
    LOG::SCOPE scope("FRAME","Frame %d",frame,1);
    std::string f=STRING_UTILITIES::string_sprintf("%d/",frame);

    local_grids.Resize(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){
        local_grids(p)=new LOCAL_GRID<T_GRID>(grid);
        FILE_UTILITIES::Read_From_File<RW>(input_directory+STRING_UTILITIES::string_sprintf("/%d/common/grid",p+fluid_proc_offset),*local_grids(p));}
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid",grid);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    if(merge_compressible){
        Merge_Cell_Data<T_ARRAYS_DIMENSION_SCALAR>(f+"euler_U",0);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"density",0);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"pressure",0);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"temperature",0);
        Merge_Cell_Data<T_ARRAYS_VECTOR>(f+"centered_velocities",0);
        if(FILE_UTILITIES::File_Exists(input_directory+"/2/"+f+"soot")) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"soot",0);
        if(FILE_UTILITIES::File_Exists(input_directory+"/2/"+f+"soot_fuel")) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"soot_fuel",0);
        if(merge_debug_data){
            Merge_Cell_Data<T_ARRAYS_BOOL>(f+"psi_D",1);
            Merge_Face_Data<T_FACE_ARRAYS_BOOL>(f+"psi_N",1);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"density_gradient",0,true);
            Merge_Cell_Data<T_ARRAYS_BOOL>(f+"euler_psi",1);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"energy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"enthalpy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"entropy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"internal_energy",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"machnumber",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"speedofsound",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"velocity_plus_c",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"velocity_minus_c",0);
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"pressure_gradient",0,true);}
        T time;FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/1/%d/time",input_directory.c_str(),frame),time);
        FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/%d/time",output_directory.c_str(),frame),time);
        std::string levelset_file=input_directory+"/1"+f+"levelset";
        if(FILE_UTILITIES::File_Exists(levelset_file)){
            Merge_Levelset<T_ARRAYS_SCALAR>(f+"levelset",3);
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/levelset_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Levelset<T_ARRAYS_SCALAR>(filename,3)) break;}}}
    else{
        if(merge_levelset){
            Merge_Levelset<T_ARRAYS_SCALAR>(f+"levelset",3);
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/levelset_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Levelset<T_ARRAYS_SCALAR>(filename,3)) break;}}
        if(merge_object_levelset) Merge_Levelset<T_ARRAYS_SCALAR>(f+"object_levelset",3);
        if(merge_debug_data){
            Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"pressure",1);
            Merge_Cell_Data<T_ARRAYS_INT>(f+"colors",0); // TODO: consider changing this back to 1
            Merge_Cell_Data<T_ARRAYS_BOOL>(f+"psi_D",1);
            Merge_Face_Data<T_FACE_ARRAYS_BOOL>(f+"psi_N",1);
            Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(f+"negative_particles");}
        if(merge_velocities) Merge_Face_Data<T_FACE_ARRAYS>(f+"mac_velocities",3);
        if(merge_temperatures) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"temperature",3);
        Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"reaction_speed",3);
        if(merge_densities) Merge_Cell_Data<T_ARRAYS_SCALAR>(f+"density",3);
        if(merge_particles){
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/bubble_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Cell_Particles<BUBBLE_PARTICLES<TV> >(filename)) break;}
            Merge_Cell_Particles<BUBBLE_PARTICLES<TV> >(f+"bubble_particles");
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/negative_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(filename)) break;}
            Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(f+"positive_particles");
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/positive_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(filename)) break;}
            Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(f+"removed_negative_particles");}
        if(merge_removed_particles){
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/removed_negative_particles_%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(filename)) break;}
            Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(f+"removed_positive_particles");
            for(int number_of_sets=1;;number_of_sets++){
                std::string filename=STRING_UTILITIES::string_sprintf("%d/removed_positive_particles_%d.%d",frame,number_of_sets);LOG::cout<<"merging "<<filename<<std::endl;
                if(!Merge_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(filename)) break;}}
        if(merge_removed_particle_times) Merge_Lists<ARRAY<PAIR<int,T> > >(f+"removed_particle_times");}
    if(merge_solid_fluid){
        solid_body_collection->rigid_body_collection.Read(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/1/",input_directory.c_str()),frame,&needs_init);
        solid_body_collection->rigid_body_collection.rigid_geometry_collection.structure_list.Fill_Needs_Write();
        solid_body_collection->rigid_body_collection.Write(STREAM_TYPE(RW()),output_directory,frame);
        bool include_static_frame=frame==0;
        solid_body_collection->deformable_body_collection.Read(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/1/",input_directory.c_str()),frame,-1,include_static_frame,false);
        solid_body_collection->deformable_body_collection.Write(STREAM_TYPE(RW()),output_directory+"/",frame,-1,include_static_frame,false);}
    if(merge_rigids) Merge_Rigid_Data(frame);
    local_grids.Delete_Pointers_And_Clean_Memory();
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+f+"/colors",colors);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Merge_Cell_Data
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_ARRAYS> bool MERGER<T,T_GRID,RW>::
Merge_Cell_Data(const std::string& filename,const int verify_bandwidth,const bool scale)
{
    // read
    ARRAY<T_ARRAYS> local_data(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
    // merge
    T_ARRAYS global_data(grid.Cell_Indices(3));
    for(int p=1;p<=number_of_fluid_processes;p++)local_grids(p)->Put(local_data(p),global_data);
    // verify
    for(int p=1;p<=number_of_fluid_processes;p++){TV_INT index;
        T max_error=local_grids(p)->Maximum_Error(local_data(p),global_data,verify_bandwidth,index);
        if(max_error>0 && print_max_errors){LOG::cout<<filename<<": max error on process "<<p<<" = "<<max_error<<" ("<<index<<" = "<<local_data(p)(index)<<")"<<std::endl;}}
    if(scale) Scale_Cell_Data(global_data);
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    return true;
}
//#####################################################################
// Function Scale_Cell_Data
//#####################################################################
template<class T,class T_GRID,class RW> void MERGER<T,T_GRID,RW>::
Scale_Cell_Data(T_ARRAYS_SCALAR& array)
{
    T max_val=array.Max();
    if(max_val)
        for(typename T_GRID::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
            array(iterator.Cell_Index())=exp(-array(iterator.Cell_Index())/max_val);
}
//#####################################################################
// Function Merge_Levelset
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_ARRAYS> bool MERGER<T,T_GRID,RW>::
Merge_Levelset(const std::string& filename,const int verify_bandwidth)
{
    // read
    ARRAY<T_ARRAYS> local_data(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        T_GRID skip;FILE_UTILITIES::Read_From_File<RW>(name,skip,local_data(p));}
    // merge
    T_ARRAYS global_data(grid.Cell_Indices(3));
    for(int p=1;p<=number_of_fluid_processes;p++)local_grids(p)->Put(local_data(p),global_data);
    // verify
    for(int p=1;p<=number_of_fluid_processes;p++){TV_INT index;
        T max_error=local_grids(p)->Maximum_Error(local_data(p),global_data,verify_bandwidth,index);
        if(max_error>0 && print_max_errors){LOG::cout<<filename<<": max error on process "<<p<<" = "<<max_error<<" ("<<index<<" = "<<local_data(p)(index)<<")"<<std::endl;}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,grid,global_data);
    return true;
}
//#####################################################################
// Function Merge_Face_Data
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_FACE_ARRAYS_2> bool MERGER<T,T_GRID,RW>::
Merge_Face_Data(const std::string& filename,const int verify_bandwidth)
{
    // read
    ARRAY<T_FACE_ARRAYS_2*> local_data(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        local_data(p)=new T_FACE_ARRAYS_2;
        FILE_UTILITIES::Read_From_File<RW>(name,*local_data(p));}
    // merge
    T_FACE_ARRAYS_2 global_data(grid,3);
    for(int p=1;p<=number_of_fluid_processes;p++)local_grids(p)->Put_Faces(*local_data(p),global_data);
    // verify
    for(int p=1;p<=number_of_fluid_processes;p++){
        if(print_max_errors){
            std::string prefix=filename+": max error on process "+STRING_UTILITIES::string_sprintf("%d",p);
            local_grids(p)->Maximum_Error(prefix,*local_data(p),global_data,verify_bandwidth);}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    local_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Merge_Rigid_Data
//#####################################################################
template<class T,class TV,class RW> bool MERGER<T,TV,RW>::
Merge_Rigid_Data(const int frame)
{    
    //int current_index=0;
    int total_number_of_threads=number_of_processes;
    local_rigid_collections.Resize(total_number_of_threads);
    for(int p=1;p<=total_number_of_threads;p++){
        local_rigid_collections(p)=new RIGID_BODY_COLLECTION<TV>(0,0);
        local_rigid_collections(p)->Read(STREAM_TYPE(RW()),STRING_UTILITIES::string_sprintf("%s/%d",input_directory.c_str(),p),frame,&needs_init);
        local_rigid_collections(p)->rigid_geometry_collection.structure_list.Fill_Needs_Write();
        rigid_body_collection.rigid_body_particle.array_collection->Add_Elements(local_rigid_collections(p)->rigid_body_particle.array_collection->Size());
        colors.Resize(local_rigid_collections(p)->rigid_body_particle.array_collection->Size());
        for(int i=1;i<=local_rigid_collections(p)->rigid_body_particle.array_collection->Size();i++){
            if(!local_rigid_collections(p)->rigid_body_particle.rigid_geometry(i)) continue;
            if(!local_rigid_collections(p)->rigid_geometry_collection.Is_Active(i)) continue;
            colors(i)=p;
            rigid_body_collection.rigid_body_particle.array_collection->Copy_Element(*local_rigid_collections(p)->rigid_body_particle.array_collection,i,i);
            rigid_body_collection.Add_Rigid_Body_And_Geometry(&local_rigid_collections(p)->Rigid_Body(i));}}
    rigid_body_collection.Write(STREAM_TYPE(RW()),output_directory,frame);
    return true;
}
//#####################################################################
// Function Merge_Lists
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_LIST_2> bool MERGER<T,T_GRID,RW>::
Merge_Lists(const std::string& filename)
{
    // read
    ARRAY<T_LIST_2*> local_data(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        local_data(p)=new T_LIST_2;
        FILE_UTILITIES::Read_From_File<RW>(name,*local_data(p));}
    // merge
    T_LIST_2 global_data;
    for(int p=1;p<=number_of_fluid_processes;p++) global_data.Append_Elements(*local_data(p));
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    local_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Merge_Particles
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_PARTICLES> bool MERGER<T,T_GRID,RW>::
Merge_Particles(const std::string& filename)
{
    // read
    int process_without_particles=0;
    ARRAY<typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE> local_data(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
    // check for zero size
    if(process_without_particles){
        for(int p=1;p<=number_of_fluid_processes;p++)if(local_data(p).counts.x){
            LOG::cerr<<filename<<": process "<<p<<" has particles but process "<<process_without_particles<<" does not."<<std::endl;exit(1);}
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE());return true;}
    // merge
    typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE global_data(grid.Node_Indices());
    for(int p=1;p<=number_of_fluid_processes;p++){
        RANGE<TV_INT> region=local_grids(p)->Interior_Region(RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector())),interior_region=region.Thickened(-1);
        // copy interior particles
        {NODE_ITERATOR local(local_grids(p)->grid,interior_region),global(grid,interior_region+local_grids(p)->offset);
        for(;local.Valid();local.Next(),global.Next())exchange(global_data(global.Node_Index()),local_data(p)(local.Node_Index()));}
        // merge boundary particles
        {NODE_ITERATOR local(local_grids(p)->grid,region),global(grid,region+local_grids(p)->offset);
        for(;local.Valid();local.Next(),global.Next()){
            T_PARTICLES *&from_particles=local_data(p)(local.Node_Index()),*&to_particles=global_data(global.Node_Index());
            if(!from_particles) continue;
            if(!to_particles){exchange(to_particles,from_particles);continue;}
            to_particles->array_collection->Append(*from_particles->array_collection);}}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    // free memory
    for(int p=1;p<=number_of_fluid_processes;p++)local_data(p).Delete_Pointers_And_Clean_Memory();
    global_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Merge_Cell_Particles
//#####################################################################
template<class T,class T_GRID,class RW> template<class T_PARTICLES> bool MERGER<T,T_GRID,RW>::
Merge_Cell_Particles(const std::string& filename)
{
    // read
    int process_without_particles=0;
    ARRAY<typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE> local_data(number_of_fluid_processes);
    for(int p=1;p<=number_of_fluid_processes;p++){std::string name=input_directory+STRING_UTILITIES::string_sprintf("/%d/",(p+fluid_proc_offset))+filename;
        if(!FILE_UTILITIES::File_Exists(name)){LOG::cout<<"Missing "<<name<<"; skipping merge"<<std::endl;return false;}
        FILE_UTILITIES::Read_From_File<RW>(name,local_data(p));}
    // check for zero size
    if(process_without_particles){
        for(int p=1;p<=number_of_fluid_processes;p++)if(local_data(p).counts.x){
            LOG::cerr<<filename<<": process "<<p<<" has particles but process "<<process_without_particles<<" does not."<<std::endl;exit(1);}
        FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE());return true;}
    // merge
    typename T_ARRAYS_SCALAR::template REBIND<T_PARTICLES*>::TYPE global_data(grid.Domain_Indices());
    for(int p=1;p<=number_of_fluid_processes;p++){
        RANGE<TV_INT> region=local_grids(p)->Interior_Region(RANGE<TV_INT>(TV_INT(),TV_INT::All_Ones_Vector())),interior_region=region.Thickened(-1);
        // copy interior particles
        {CELL_ITERATOR local(local_grids(p)->grid,interior_region),global(grid,interior_region+local_grids(p)->offset);
        for(;local.Valid();local.Next(),global.Next())exchange(global_data(global.Cell_Index()),local_data(p)(local.Cell_Index()));}
        // merge boundary particles
        {CELL_ITERATOR local(local_grids(p)->grid,region),global(grid,region+local_grids(p)->offset);
        for(;local.Valid();local.Next(),global.Next()){
            T_PARTICLES *&from_particles=local_data(p)(local.Cell_Index()),*&to_particles=global_data(global.Cell_Index());
            if(!from_particles) continue;
            if(!to_particles){exchange(to_particles,from_particles);continue;}
            to_particles->array_collection->Append(*from_particles->array_collection);}}}
    // write
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+filename,global_data);
    // free memory
    for(int p=1;p<=number_of_fluid_processes;p++)local_data(p).Delete_Pointers_And_Clean_Memory();
    global_data.Delete_Pointers_And_Clean_Memory();
    return true;
}
//#####################################################################
// Function Do_Merge
//#####################################################################
template<class T> void
Do_Merge(PARSE_ARGS& parse_args)
{
    std::string input_directory=parse_args.Extra_Arg(1),output_directory=input_directory;
    if(parse_args.Is_Value_Set("-o")) output_directory=parse_args.Get_String_Value("-o");

    int first_frame,last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/common/first_frame",first_frame);
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/1/common/last_frame",last_frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",first_frame);
    bool print_maxerrors=parse_args.Is_Value_Set("-print_maxerrors");
    bool is_solid_fluid=parse_args.Is_Value_Set("-solid_fluid");
    bool merge_rigids=parse_args.Is_Value_Set("-merge_rigid");
    bool is_1d=false,is_2d=false,is_3d=false;

    if(parse_args.Is_Value_Set("-start_frame")) first_frame=parse_args.Get_Integer_Value("-start_frame");
    if(parse_args.Is_Value_Set("-last_frame")) last_frame=parse_args.Get_Integer_Value("-last_frame");
    if(parse_args.Is_Value_Set("-1d")) is_1d=true;
    else if(parse_args.Is_Value_Set("-2d")) is_2d=true;
    else if(parse_args.Is_Value_Set("-3d")) is_3d=true;
    //bool force=parse_args.Get_Option_Value("-f");
    int number_of_processes=parse_args.Get_Integer_Value("-np");
    if(number_of_processes<0){LOG::cerr<<"Invalid np<0"<<std::endl;exit(1);}
    else if(!number_of_processes){ // autodetect number of processes
        FILE_UTILITIES::Find_First_Nonexistent_Directory_In_Sequence(STRING_UTILITIES::string_sprintf("%s/%%d",input_directory.c_str()),1,&number_of_processes);--number_of_processes;
        LOG::cout<<"Autodetected "<<number_of_processes<<" processes"<<std::endl;}

    if(!(is_1d || is_2d || is_3d)) {
        VECTOR<int,3> mnmn;
        if(!is_solid_fluid) FILE_UTILITIES::Read_From_File<T>(input_directory+"/1/common/grid",mnmn);
        else FILE_UTILITIES::Read_From_File<T>(input_directory+"/2/common/grid",mnmn);
        if(mnmn.z>10 && mnmn.z<20000) is_3d=true;
        else if(mnmn.y>10 && mnmn.y<20000) is_2d=true;
        else is_1d=true;}

    if(is_1d){
        MERGER<T,GRID<VECTOR<T,1> >,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors,is_solid_fluid,merge_rigids);
        if(parse_args.Is_Value_Set("-skip_levelset"))merger.merge_levelset=false;
        if(parse_args.Is_Value_Set("-skip_object_levelset"))merger.merge_object_levelset=false;
        if(parse_args.Is_Value_Set("-skip_debug_data"))merger.merge_debug_data=false;
        if(parse_args.Is_Value_Set("-skip_particles"))merger.merge_particles=false;
        if(parse_args.Is_Value_Set("-skip_removed_particles"))merger.merge_removed_particles=false;
        if(parse_args.Is_Value_Set("-removed_particle_times"))merger.merge_removed_particle_times=true;
        if(parse_args.Is_Value_Set("-skip_velocities"))merger.merge_velocities=false;
        if(parse_args.Is_Value_Set("-compressible"))merger.merge_compressible=true;
        if(parse_args.Is_Value_Set("-minimal")){
            merger.merge_object_levelset=false;merger.merge_debug_data=false;merger.merge_particles=false;merger.merge_removed_particles=false;merger.merge_removed_particle_times=false;merger.merge_velocities=false;}
        merger.Merge_All_Frames(first_frame,last_frame);}
    else if(is_2d){
        MERGER<T,GRID<VECTOR<T,2> >,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors,is_solid_fluid,merge_rigids);
        if(parse_args.Is_Value_Set("-skip_levelset"))merger.merge_levelset=false;
        if(parse_args.Is_Value_Set("-skip_object_levelset"))merger.merge_object_levelset=false;
        if(parse_args.Is_Value_Set("-skip_debug_data"))merger.merge_debug_data=false;
        if(parse_args.Is_Value_Set("-skip_particles"))merger.merge_particles=false;
        if(parse_args.Is_Value_Set("-skip_removed_particles"))merger.merge_removed_particles=false;
        if(parse_args.Is_Value_Set("-removed_particle_times"))merger.merge_removed_particle_times=true;
        if(parse_args.Is_Value_Set("-skip_velocities"))merger.merge_velocities=false;
        if(parse_args.Is_Value_Set("-compressible"))merger.merge_compressible=true;
        if(parse_args.Is_Value_Set("-minimal")){
            merger.merge_object_levelset=false;merger.merge_debug_data=false;merger.merge_particles=false;merger.merge_removed_particles=false;merger.merge_removed_particle_times=false;merger.merge_velocities=false;}
        merger.Merge_All_Frames(first_frame,last_frame);}
    else if(is_3d){
        MERGER<T,GRID<VECTOR<T,3> >,T> merger(number_of_processes,input_directory,output_directory,print_maxerrors,is_solid_fluid,merge_rigids);
        if(parse_args.Is_Value_Set("-skip_levelset"))merger.merge_levelset=false;
        if(parse_args.Is_Value_Set("-skip_object_levelset"))merger.merge_object_levelset=false;
        if(parse_args.Is_Value_Set("-skip_debug_data"))merger.merge_debug_data=false;
        if(parse_args.Is_Value_Set("-skip_particles"))merger.merge_particles=false;
        if(parse_args.Is_Value_Set("-skip_removed_particles"))merger.merge_removed_particles=false;
        if(parse_args.Is_Value_Set("-removed_particle_times"))merger.merge_removed_particle_times=true;
        if(parse_args.Is_Value_Set("-skip_velocities"))merger.merge_velocities=false;
        if(parse_args.Is_Value_Set("-compressible"))merger.merge_compressible=true;
        if(parse_args.Is_Value_Set("-minimal")){
            merger.merge_object_levelset=false;merger.merge_debug_data=false;merger.merge_particles=false;merger.merge_removed_particles=false;merger.merge_removed_particle_times=false;merger.merge_velocities=false;}
        merger.Merge_All_Frames(first_frame,last_frame);}
    LOG::cout<<std::endl;
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc,char* argv[])
{
    Initialize_Particles();
    Initialize_Read_Write_Structures();

    LOG::Initialize_Logging(false,false,1<<30,true);

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-f","force");
    parse_args.Add_Integer_Argument("-start_frame",0,"start frame number");
    parse_args.Add_Integer_Argument("-last_frame",0,"last frame number");
    parse_args.Add_Integer_Argument("-np",0,"number of mpi processes (use 0 to autodetect)");
    parse_args.Add_String_Argument("-o","","output directory");
    parse_args.Add_Option_Argument("-skip_levelset","skip_levelset");
    parse_args.Add_Option_Argument("-skip_object_levelset","skip_object_levelset");
    parse_args.Add_Option_Argument("-skip_debug_data","skip_debug_data");
    parse_args.Add_Option_Argument("-skip_particles","skip_particles");
    parse_args.Add_Option_Argument("-skip_removed_particles","skip_removed_particles");
    parse_args.Add_Option_Argument("-removed_particle_times","removed_particle_times");
    parse_args.Add_Option_Argument("-skip_velocities","skip_velocities");
    parse_args.Add_Option_Argument("-minimal","skip everything but the levelset");
    parse_args.Add_Option_Argument("-print_maxerrors","print max errors");
    parse_args.Add_Option_Argument("-compressible","input data is compressible output");
    parse_args.Add_Option_Argument("-solid_fluid","input data is solid fluid data");
    parse_args.Add_Option_Argument("-merge_rigid","input data is rigid body collection");
    parse_args.Add_Option_Argument("-double","input data is in doubles");
    parse_args.Add_Option_Argument("-1d","input data is 1-D");
    parse_args.Add_Option_Argument("-2d","input data is 2-D");
    parse_args.Add_Option_Argument("-3d","input data is 3-D");
    parse_args.Set_Extra_Arguments(1,"<input_directory>");
    parse_args.Parse(argc,argv);

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(parse_args.Is_Value_Set("-double")) Do_Merge<double>(parse_args); 
    else Do_Merge<float>(parse_args);
#else
    Do_Merge<float>(parse_args);
#endif
}
//#####################################################################
