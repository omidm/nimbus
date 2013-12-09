/*
 * The application specification of PhysBAM projection.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "projection_driver.h"
#include "projection_example.h"

#include "data_impl.h"
#include "job_impl.h"

#include "app.h"

typedef float T;

void App::Load() {
	dbg_add_mode("proj");
	dbg(DBG_PROJ, "Starts to load projection application.\n");
	RegisterJob("main", new Main(this));
	RegisterJob("Init", new Init(this));
	RegisterJob("Project_Forloop_Condition", new Project_Forloop_Condition(this));
	RegisterJob("Project_Forloop_Part1", new Project_Forloop_Part1(this));
	RegisterJob("Project_Forloop_Part2", new Project_Forloop_Part2(this));
	RegisterJob("Project_Forloop_Part3", new Project_Forloop_Part3(this));
	RegisterJob("Project_Forloop_Part4", new Project_Forloop_Part4(this));
	RegisterJob("Global_Sum", new Global_Sum(this));
	RegisterJob("Global_Max_Abs", new Global_Max_Abs(this));
	RegisterJob("Finish", new Finish(this));

	// Simulates parameter passing mechanism.
	// [TODO] Revise or delete.
	char a[] = "multiple_projection";
	char* temp[2];
	temp[0] = a;
	temp[1] = NULL;
	InitMain(1, temp);
	app_driver->PrepareForProjection();
	app_driver->PrepareForOneRegion();
	dbg(DBG_PROJ, "Starts to load projection application.\n");
}

// Create driver and example. Expected to be changed in the future.
void App::InitMain(int argc, char* argv[]) {
	using namespace PhysBAM;
	typedef float T;
	typedef float RW;
	STREAM_TYPE *stream_type = new STREAM_TYPE((RW()));

	PARSE_ARGS parse_args;
	parse_args.Add_Integer_Argument("-restart", 0, "restart frame");
	parse_args.Add_Integer_Argument("-scale", 64, "fine scale grid resolution");
	parse_args.Add_Integer_Argument("-substep", -1, "output-substep level");
	parse_args.Add_Integer_Argument("-threads", 1, "number of threads");
	parse_args.Add_Option_Argument("-3d", "run in 3 dimensions");

	LOG::Initialize_Logging(false, false, 1<<30, true, parse_args.Get_Integer_Value("-threads"));
	new MPI_WORLD(argc, argv);

	parse_args.Parse(argc, argv);
	parse_args.Print_Arguments(argc, argv);

	typedef VECTOR<T,2> TV;
	typedef VECTOR<int,TV::dimension> TV_INT;

	PROJECTION_EXAMPLE<TV>* example=new PROJECTION_EXAMPLE<TV>(*stream_type,parse_args.Get_Integer_Value("-threads"));

	int scale=parse_args.Get_Integer_Value("-scale");
	RANGE<TV> range(TV(), TV::All_Ones_Vector()*0.5);
	range.max_corner(2)=1;
	TV_INT counts=TV_INT::All_Ones_Vector()*scale/2;
	counts(2)=scale;
	example->Initialize_Grid(counts, range);
	example->restart=parse_args.Get_Integer_Value("-restart");
	example->write_substeps_level=parse_args.Get_Integer_Value("-substep");


	example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);	
	example->output_directory+=STRING_UTILITIES::string_sprintf("/%d", _rankID);


	FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
	LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt", false);

	app_driver = new PhysBAM::PROJECTION_DRIVER<TV>(*example);
}

void App::FinishMain() {
	delete app_driver;	
}

/*
 // The workflow for projection.
 PROJECTION_DRIVER<TV> &driver;
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
 */

