#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "projection_driver.h"
#include "projection_example.h"
#include "init_main.h"

using namespace PhysBAM;

void InitMain(int argc,char* argv[])
{
  typedef float T;
  typedef float RW;
  STREAM_TYPE *stream_type = new STREAM_TYPE((RW()));

  PARSE_ARGS parse_args;
  parse_args.Add_Integer_Argument("-restart",0,"restart frame");
  parse_args.Add_Integer_Argument("-scale",64,"fine scale grid resolution");
  parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
  parse_args.Add_Integer_Argument("-threads",1,"number of threads");
  parse_args.Add_Option_Argument("-3d","run in 3 dimensions");

  LOG::Initialize_Logging(false,false,1<<30,true,parse_args.Get_Integer_Value("-threads"));
  mpi_world = new MPI_WORLD(argc,argv);

  parse_args.Parse(argc,argv);
  parse_args.Print_Arguments(argc,argv);
   
  typedef VECTOR<T,2> TV;
  typedef VECTOR<int,TV::dimension> TV_INT;
  
  PROJECTION_EXAMPLE<TV>* example=new PROJECTION_EXAMPLE<TV>(*stream_type,parse_args.Get_Integer_Value("-threads"));

  int scale=parse_args.Get_Integer_Value("-scale");
  RANGE<TV> range(TV(),TV::All_Ones_Vector()*0.5);range.max_corner(2)=1;TV_INT counts=TV_INT::All_Ones_Vector()*scale/2;counts(2)=scale;
  example->Initialize_Grid(counts,range);
  example->restart=parse_args.Get_Integer_Value("-restart");
  example->write_substeps_level=parse_args.Get_Integer_Value("-substep");

  if(mpi_world->initialized){
      example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
      if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->mpi_grid->Number_Of_Processors(),(mpi_world->rank+1));}

  FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
  LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
  
  driver = new PhysBAM::PROJECTION_DRIVER<TV>(*example);
}
void FinishMain() {
  delete driver;
  delete mpi_world;
}
