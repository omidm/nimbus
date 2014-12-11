#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"
#include <cstdio>
#include <ctime>
#include <time.h>

using namespace PhysBAM;

template<class T>
void Add_Source(WATER_EXAMPLE<VECTOR<T,1> >* example)
{
    PHYSBAM_FATAL_ERROR();
}

template<class T>
void Add_Source(WATER_EXAMPLE<VECTOR<T,2> >* example)
{
    typedef VECTOR<T,2> TV;
    TV point1,point2;BOX<TV> source;
    //point1=TV::All_Ones_Vector()*(T).5;point1(1)=.4;point1(2)=.95;point2=TV::All_Ones_Vector()*(T).65;point2(1)=.55;point2(2)=1;
    point1=TV::All_Ones_Vector()*(T).5;point1(1)=.95;point1(2)=.6;point2=TV::All_Ones_Vector()*(T).65;point2(1)=1;point2(2)=.75;
    source.min_corner=point1;source.max_corner=point2;
    example->sources.Append(new ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >(source));
}

template<class T>
void Add_Source(WATER_EXAMPLE<VECTOR<T,3> >* example)
{
    typedef VECTOR<T,3> TV;
    TV point1,point2;CYLINDER<T> source;
    //point1=TV::All_Ones_Vector()*(T).6;point1(1)=.4;point1(2)=.95;point2=TV::All_Ones_Vector()*(T).6;point2(1)=.4;point2(2)=1;
    point1=TV::All_Ones_Vector()*(T).8;point1(1)=.4;point1(3)=.95;point2=TV::All_Ones_Vector()*(T).8;point2(1)=.4;point2(3)=1;
    source.Set_Endpoints(point1,point2);source.radius=.1;
    IMPLICIT_OBJECT<TV>* analytic=new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(source);
    example->sources.Append(analytic);
}

template<class T>
void Add_Sphere(WATER_EXAMPLE<VECTOR<T,1> >* example)
{
    PHYSBAM_FATAL_ERROR();
}

template<class T>
void Add_Sphere(WATER_EXAMPLE<VECTOR<T,2> >* example)
{
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    std::string model_file_name="../../Public_Data/Archives/Rigid_Bodies_2D/circle";
    T radius=0.1;
    int rigid_particle_id=example->rigid_geometry_collection.Add_Rigid_Geometry(stream_type,model_file_name,radius,true,true,true,true);
    example->rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->is_static=false;
}

template<class T>
void Add_Sphere(WATER_EXAMPLE<VECTOR<T,3> >* example)
{
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    std::string model_file_name="../../Public_Data/Archives/Rigid_Bodies/sphere";
    T radius=0.1;
    int rigid_particle_id=example->rigid_geometry_collection.Add_Rigid_Geometry(stream_type,model_file_name,radius,true,true,true,true);
    example->rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->is_static=false;
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    //typedef VECTOR<T,2> TV;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::dimension> TV_INT;


    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-scale",128,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-e",1,"last frame");
    parse_args.Add_Integer_Argument("-refine",1,"refine levels");
    parse_args.Add_Integer_Argument("-threads",1,"number of threads");
    parse_args.Add_Double_Argument("-cfl",1,"cfl number");

    LOG::Initialize_Logging(false,false,1<<30,true,parse_args.Get_Integer_Value("-threads"));
    MPI_WORLD mpi_world(argc,argv);

    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
    
    WATER_EXAMPLE<TV>* example=new WATER_EXAMPLE<TV>(stream_type,parse_args.Get_Integer_Value("-threads"),parse_args.Get_Integer_Value("-refine"));

    int scale=parse_args.Get_Integer_Value("-scale");
    example->Initialize_Grid(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()));
    example->restart=parse_args.Get_Integer_Value("-restart");
    example->last_frame=parse_args.Get_Integer_Value("-e");
    example->write_substeps_level=parse_args.Get_Integer_Value("-substep");
    example->cfl=parse_args.Get_Double_Value("-cfl");
    Add_Source(example);
    //Add_Sphere(example);

    // Custom Partition
    TV_INT ppd=TV_INT::All_Ones_Vector();
    //ppd(1)=4;

    if(mpi_world.initialized){
        char buffer[100];
        int n = sprintf(buffer, "mpi%d.log", mpi_world.rank+1);
        buffer[n] = '\0';
	char pack_buffer[100];
	n = sprintf(pack_buffer, "pack_mpi%d.log", mpi_world.rank+1);
        pack_buffer[n] = '\0';
        PhysBAM::MPI_UTILITIES::SetLogFile(fopen(buffer, "w"));
	PhysBAM::MPI_UTILITIES::SetPackLogFile(fopen(pack_buffer, "w"));
        PhysBAM::MPI_UTILITIES::PrintTimeStamp("RESTART");
        // Custom Partition
        // example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3,false,ppd);
        // Original Partition
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    WATER_DRIVER<TV> driver(*example);

    driver.Execute_Main_Program();

    return 0;
}
