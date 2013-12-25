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
#ifndef __WATER_DRIVER__  // NOLINT
#define __WATER_DRIVER__  // NOLINT

#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
namespace PhysBAM {

template<class TV> class WATER_EXAMPLE;

template<class TV>
class WATER_DRIVER {
 public:
  typedef typename TV::SCALAR T;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;
  typedef EXTRAPOLATION_UNIFORM<GRID<TV>, T> T_EXTRAPOLATION_SCALAR;
  typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_BASE T_ARRAYS_BASE;

  int current_frame;
  T time;
  int output_number;

  WATER_EXAMPLE<TV>& example;

  WATER_DRIVER(WATER_EXAMPLE<TV>& example) : example(example) {}  // NOLINT
  virtual ~WATER_DRIVER() {}

  void Scalar_Advance(const T dt, const T time);
  void Convect(const T dt, const T time);
  void Add_Forces(const T dt, const T time);
  void Project(const T dt, const T time);
  void Execute_Main_Program();
  void Initialize();
  void Advance_To_Target_Time(const T target_time);
  void Simulate_To_Frame(const int frame_input);
  void Write_Output_Files(const int frame);
};

}  // namespace PhysBAM
#endif  // NOLINT
