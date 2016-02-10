#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "PROJECTION_DRIVER.h"
#include "PROJECTION_EXAMPLE.h"

using namespace PhysBAM;

template<class TV> void Execute_Main_Program(STREAM_TYPE& stream_type,PARSE_ARGS& parse_args,MPI_WORLD& mpi_world)
{
    typedef VECTOR<int,TV::dimension> TV_INT;
    
    PROJECTION_EXAMPLE<TV>* example=new PROJECTION_EXAMPLE<TV>(stream_type,parse_args.Get_Integer_Value("-threads"));

    int scale=parse_args.Get_Integer_Value("-scale");
    RANGE<TV> range(TV(),TV::All_Ones_Vector()*0.5);range.max_corner(2)=1;TV_INT counts=TV_INT::All_Ones_Vector()*scale/2;counts(2)=scale;
    example->Initialize_Grid(counts,range);
    example->restart=parse_args.Get_Integer_Value("-restart");
    example->write_substeps_level=parse_args.Get_Integer_Value("-substep");

    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    PROJECTION_DRIVER<TV> driver(*example);
//    driver.Execute_Main_Program();

  // Build projection_data and projection_internal_data implicitly.
  driver.PrepareForProjection();
  driver.PrepareForOneRegion();

  driver.pcg_mpi->ExchangePressure(driver.projection_internal_data, driver.projection_data);
  driver.pcg_mpi->InitializeResidual(driver.projection_internal_data, driver.projection_data);
  driver.pcg_mpi->SpawnFirstIteration(driver.projection_internal_data, driver.projection_data);
  if (driver.projection_internal_data->move_on) {
    driver.projection_internal_data->iteration = 0;
    do {
      driver.projection_internal_data->iteration++;
      driver.pcg_mpi->DoPrecondition(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->CalculateBeta(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->UpdateSearchVector(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->ExchangeSearchVector(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->UpdateTempVector(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->CalculateAlpha(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->UpdateOtherVectors(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->CalculateResidual(driver.projection_internal_data, driver.projection_data);
      driver.pcg_mpi->DecideToSpawnNextIteration(driver.projection_internal_data, driver.projection_data);
    } while (driver.projection_internal_data->move_on);
    driver.pcg_mpi->ExchangePressure(driver.projection_internal_data, driver.projection_data);
  }
  driver.WindUpForOneRegion();
  driver.ApplyPressureAndFinish();
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));


    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-scale",100,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-threads",1,"number of threads");
    parse_args.Add_Option_Argument("-3d","run in 3 dimensions");

    LOG::Initialize_Logging(false,false,1<<30,true,parse_args.Get_Integer_Value("-threads"));
    MPI_WORLD mpi_world(argc,argv);

    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
     
    if(parse_args.Is_Value_Set("-3d")){
        Execute_Main_Program<VECTOR<T,3> >(stream_type,parse_args,mpi_world);}
    else{
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);}

    return 0;
}
