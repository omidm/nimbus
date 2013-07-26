//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_DRIVER__
#define __WATER_DRIVER__
#include "myinclude.h"
namespace PhysBAM {

class WATER_EXAMPLE;

// typedef VECTOR<float, 2> template_TV;
typedef VECTOR<float, 3> template_TV;

class WATER_DRIVER {
public:
  typedef template_TV TV;
  typedef typename TV::SCALAR T;
  typedef typename TV::REBIND<int>::TYPE TV_INT;
  typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET; //?
  typedef typename ADVECTION_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR; //?
  typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
  typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
  typedef typename T_ARRAYS_SCALAR::REBIND<bool>::TYPE T_ARRAYS_BOOL;
  typedef typename T_FACE_ARRAYS_SCALAR::REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;

  int current_frame;
  T time;
  int output_number;

  WATER_EXAMPLE& example;
  KINEMATIC_EVOLUTION<TV> kinematic_evolution;
public:
  THREAD_QUEUE* thread_queue;

  WATER_DRIVER(WATER_EXAMPLE& example);
  virtual ~WATER_DRIVER();

  void Execute_Main_Program();
  void Initialize();
  void Advance_To_Target_Time(const T target_time);
  void Simulate_To_Frame(const int frame_input, const int tid = 1);
  void Write_Output_Files(const int frame);
  void Write_Substep(const std::string& title, const int substep,
      const int level = 0);
  void Run(RANGE<TV_INT>& domain, const T dt, const T time);

  // Used by Hang Qu.
  friend void* advect_velocity_worker(void* arg);
  friend void* advect_velocity_fetcher(void* arg);

  struct ADVECT_VELOCITY_WORKER_T {
    typedef AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > > AVERAGING_TYPE;
    typedef LINEAR_INTERPOLATION_UNIFORM<GRID<TV>, float,
        FACE_LOOKUP_UNIFORM<GRID<TV> > > INTERPOLATION_TYPE;

    TV_INT range_all, range_re, range_x, range_y, range_z;
    int segment_len;
    T my_dt;
    T_FACE_ARRAYS_SCALAR *my_face_velocities_ghost;

    bool fetcher_stop;
  } ADVECT_VELOCITY_WORKER;
//#####################################################################
};
}
#endif
