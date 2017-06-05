/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "WATER_DRIVER.h"  // NOLINT
#include "WATER_EXAMPLE.h"  // NOLINT

namespace PhysBAM {
template<class TV> void Execute_Main_Program(
    STREAM_TYPE& stream_type,  // NOLINT
    PARSE_ARGS& parse_args, MPI_WORLD& mpi_world) {
  typedef VECTOR<int, TV::dimension> TV_INT;

  WATER_EXAMPLE<TV>* example = new WATER_EXAMPLE<TV>(stream_type,
      parse_args.Get_Integer_Value("-threads"));

  int scale = parse_args.Get_Integer_Value("-scale");
  RANGE<TV> range(TV(), TV::All_Ones_Vector());
  TV_INT counts = TV_INT::All_Ones_Vector() * scale;
  example->Initialize_Grid(counts, range);
  example->restart = parse_args.Get_Integer_Value("-restart");
  example->last_frame = parse_args.Get_Integer_Value("-e");
  example->write_substeps_level = parse_args.Get_Integer_Value("-substep");
  example->write_debug_data = true;

  TV point1 = TV::All_Ones_Vector() * .65, point2 = TV::All_Ones_Vector() * .75;
  point1(1) = 0;
  point2(1) = .05;
  example->source.min_corner = point1;
  example->source.max_corner = point2;

  if (mpi_world.initialized) {
    example->mpi_grid = new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid, 3);
    if (example->mpi_grid->Number_Of_Processors() > 1)
      example->output_directory += STRING_UTILITIES::string_sprintf("/%d",
          (mpi_world.rank + 1));
  }

  FILE_UTILITIES::Create_Directory(example->output_directory + "/common");

  WATER_DRIVER<TV> driver(*example);
  driver.Execute_Main_Program();
}
}  // namespace PhysBAM

int main(int argc, char *argv[]) {
  using namespace PhysBAM;  // NOLINT
  typedef float T;
  typedef float RW;
  STREAM_TYPE stream_type((RW()));

  PARSE_ARGS parse_args;
  parse_args.Add_Integer_Argument("-restart", 0, "restart frame");
  parse_args.Add_Integer_Argument("-scale", 100, "fine scale grid resolution");
  parse_args.Add_Integer_Argument("-substep", -1, "output-substep level");
  parse_args.Add_Integer_Argument("-e", 100, "last frame");
  parse_args.Add_Integer_Argument("-threads", 1, "number of threads");
  parse_args.Add_Option_Argument("-3d", "run in 3 dimensions");

  parse_args.Parse(argc, argv);
  parse_args.Print_Arguments(argc, argv);

  LOG::Initialize_Logging(true, true, 1 << 30, true,
      parse_args.Get_Integer_Value("-threads"));

  MPI_WORLD mpi_world(argc, argv);

  printf("Start\n");

  if (parse_args.Is_Value_Set("-3d")) {
    Execute_Main_Program<VECTOR<T, 3> >(stream_type, parse_args, mpi_world);
  } else {
    Execute_Main_Program<VECTOR<T, 2> >(stream_type, parse_args, mpi_world);
  }

  return 0;
}
