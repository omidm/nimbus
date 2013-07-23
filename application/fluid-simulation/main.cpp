#include "public/fluid-simulation-app.h"

#include <cassert>

#include "global-repo.h"
#include "myinclude.h"
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"

using namespace PhysBAM;

struct GlobalRepo *g_global_repo = NULL;

// Add the fluid source to the example.
void Add_Source(WATER_EXAMPLE* example) {
  typedef float T;
  // typedef VECTOR<T, 2> TV;
  typedef VECTOR<T, 3> TV;
  TV point1, point2;
  CYLINDER<T> source;
 
  /* 
     // Initialization for 2D case.
     point1 = TV::All_Ones_Vector()*(T).6;
     point1(1) = .4;
     point1(2) = .95;
     point2 = TV::All_Ones_Vector()*(T).6;
     point2(1) = .4;
     point2(2) = 1;
   */

  // Initialization for 3D case.
  point1 = TV::All_Ones_Vector() * (T) .8;
  point1(1) = .4;
  point1(3) = .95;
  point2 = TV::All_Ones_Vector() * (T) .8;
  point2(1) = .4;
  point2(3) = 1;
  source.Set_Endpoints(point1, point2);
  source.radius = .1;
  IMPLICIT_OBJECT<TV> *analytic = new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(
      source);
  example->sources.Append(analytic);
}

WATER_DRIVER *g_water_driver = NULL;

int main_job(int argc, char **argv) {
  typedef float T;
  typedef float RW;
  STREAM_TYPE stream_type((RW()));
  // typedef VECTOR<T, 2> TV;
  typedef VECTOR<T, 3> TV;
  typedef VECTOR<int, TV::dimension> TV_INT;

  PARSE_ARGS parse_args;
  parse_args.Add_Integer_Argument("-restart", 0, "restart frame");
  parse_args.Add_Integer_Argument("-scale", 128, "fine scale grid resolution");
  parse_args.Add_Integer_Argument("-substep", -1, "output-substep level");
  parse_args.Add_Integer_Argument("-e", 100, "last frame");
  parse_args.Add_Integer_Argument("-refine", 1, "refine levels");
  parse_args.Add_Integer_Argument("-threads", 1, "number of threads");
  parse_args.Add_Integer_Argument("-worker", 7, "number of workers");
  parse_args.Add_Integer_Argument("-cut", 20, "cut length");
  parse_args.Add_Double_Argument("-cfl", 1, "cfl number");

  LOG::Initialize_Logging(false, false, 1 << 30, true,
      parse_args.Get_Integer_Value("-threads"));

  parse_args.Parse(argc, argv);
  parse_args.Print_Arguments(argc, argv);

  WATER_EXAMPLE* example = new WATER_EXAMPLE(stream_type,
      parse_args.Get_Integer_Value("-threads"),
      parse_args.Get_Integer_Value("-refine"));

  int scale = parse_args.Get_Integer_Value("-scale");
  example->Initialize_Grid(TV_INT::All_Ones_Vector() * scale,
      RANGE < TV > (TV(), TV::All_Ones_Vector()));
  example->restart = parse_args.Get_Integer_Value("-restart");
  example->last_frame = parse_args.Get_Integer_Value("-e");
  example->write_substeps_level = parse_args.Get_Integer_Value("-substep");
  example->cfl = parse_args.Get_Double_Value("-cfl");
  Add_Source(example);

  FILE_UTILITIES::Create_Directory(example->output_directory + "/common");
  LOG::Instance()->Copy_Log_To_File(
      example->output_directory + "/common/log.txt", false);

  WATER_DRIVER *driver = new WATER_DRIVER(*example);
  driver->ADVECT_VELOCITY_WORKER.segment_len = parse_args.Get_Integer_Value(
      "-cut");
  driver->ADVECT_VELOCITY_WORKER.worker_num = parse_args.Get_Integer_Value(
      "-worker");

  assert(g_global_repo == NULL);
  g_global_repo = new struct GlobalRepo;
  g_global_repo->water_example = example;
  g_global_repo->water_driver = driver;

  return 0;
}

void run_job() {
  g_global_repo->water_driver->Execute_Main_Program();
} 
