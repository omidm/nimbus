#include "myinclude.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <sched.h>
#include <errno.h>
#include <string.h>
#include <netdb.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <string>
#include <sstream>

#include <arpa/inet.h>
#include "advection-velocity.h"
#include "WATER_EXAMPLE.h"
#include "WATER_DRIVER.h"

using namespace PhysBAM;

namespace ADVECT_VELOCITY_NS {

// Set the core affinity of a thread to core_id.
void set_core_affinity(int core_id) {
  cpu_set_t temp_set;
  CPU_ZERO(&temp_set);
  CPU_SET(core_id, &temp_set);
  sched_setaffinity(0, sizeof(temp_set), &temp_set);
}

void run(WATER_DRIVER &driver, WATER_EXAMPLE &example,
    WATER_DRIVER::ADVECT_VELOCITY_WORKER_T &INFO, AVERAGING_TYPE &averaging,
    INTERPOLATION_TYPE &interpolation, const TV_INT &segment_start) {
  TV_INT face;
  int axis;
  int x, y, z;
  face = segment_start;
  // [TODO] Add some neat trick to handle the range here! Too ugly.
  for (x = segment_start(1); x < segment_start(1) + INFO.segment_len; x++)
    if (x < INFO.range_re(1))
      for (y = segment_start(2); y < segment_start(2) + INFO.segment_len; y++)
        if (y < INFO.range_re(2))
          for (z = segment_start(3); z < segment_start(3) + INFO.segment_len;
              z++)
            if (z < INFO.range_re(3))
              for (axis = 1; axis <= 3; axis++) {
                face = TV_INT(x, y, z);
                if (((axis == 1)
                    && ((face(1) < INFO.range_x(1))
                        && (face(2) < INFO.range_x(2))
                        && (face(3) < INFO.range_x(3))))
                    || ((axis == 2)
                        && ((face(1) < INFO.range_y(1))
                            && (face(2) < INFO.range_y(2))
                            && (face(3) < INFO.range_y(3))))
                    || ((axis == 3)
                        && ((face(1) < INFO.range_z(1))
                            && (face(2) < INFO.range_z(2))
                            && (face(3) < INFO.range_z(3))))) {
                  if (example.particle_levelset_evolution.phi(face)
                      <= example.mac_grid.min_dX * 2.)
                    // Calculate the velocity for the cell, index == FACE, axis = AXIS.
                    //printf("%d %d %d %d;",face(1), face(2), face(3), axis);
                    // The coupling is limited to usage of GRID, FACE_ARRAY.
                    example.face_velocities.Component(axis)(face) = // Write to the output of face_velocities!
                        interpolation.Clamped_To_Array_Face_Component(
                            // interpolation calculation.
                            axis, // part of index.
                            example.mac_grid, // Context.
                            ((typename AVERAGING_TYPE::FACE_LOOKUP) (*INFO.my_face_velocities_ghost)).Starting_Point_Face(
                                axis, face), // Read the velocities.
                            example.mac_grid.Face(axis, face) // Context, function of index.
                                - INFO.my_dt
                                    * averaging.Face_To_Face_Vector(
                                        example.mac_grid, axis, face,
                                        (typename AVERAGING_TYPE::FACE_LOOKUP) (*INFO.my_face_velocities_ghost)));
                  // Read the velocities.
                } // Running task over.
              }
}

/*
void* advect_velocity_worker(void *arg) {
  // Import water example and water driver.
  typedef typename WATER_DRIVER::ADVECT_VELOCITY_WORKER_T::ThreadInfo ThreadInfo;
  ThreadInfo *tinfo = (ThreadInfo*) arg;
  WATER_DRIVER &driver = *(tinfo->driver);
  WATER_EXAMPLE &example = driver.example;
  typename WATER_DRIVER::ADVECT_VELOCITY_WORKER_T &INFO =
      driver.ADVECT_VELOCITY_WORKER;

  // Import main function.
  AVERAGING_TYPE averaging;
  INTERPOLATION_TYPE interpolation;

  set_core_affinity(tinfo->assigned_core_num);

  // Now this vector is equivalant to the task.
  TV_INT segment_start;

  // Set the core affinity of a thread.
  while (1) {
    // Fetch a task, which is described as the buffer_segment_start.
    pthread_mutex_lock(&INFO.mutex_buffer);
    while (INFO.task_exec_buffer->top == 0)
      pthread_cond_wait(&INFO.cond_buffer_any, &INFO.mutex_buffer);
    INFO.ongoing_worker_num++;
    segment_start = INFO.task_exec_buffer->task_content[--INFO.task_exec_buffer->top];
    if (INFO.task_exec_buffer->top == 0)
      pthread_cond_signal(&INFO.cond_buffer_clear);
    pthread_mutex_unlock(&INFO.mutex_buffer);

    // Run the task.
    run(driver, example, INFO, averaging, interpolation, segment_start);

    pthread_mutex_lock(&INFO.mutex_buffer);
    INFO.ongoing_worker_num--;
    pthread_cond_signal(&INFO.cond_finish);
    pthread_mutex_unlock(&INFO.mutex_buffer);
  }
  return NULL;
}

void* advect_velocity_fetcher(void *arg) {
  typedef typename WATER_DRIVER::ADVECT_VELOCITY_WORKER_T::ThreadInfo ThreadInfo;
  ThreadInfo *tinfo = (ThreadInfo*) arg;
  WATER_DRIVER &driver = *(tinfo->driver);
  WATER_EXAMPLE &example = driver.example;
  typename WATER_DRIVER::ADVECT_VELOCITY_WORKER_T &INFO =
      driver.ADVECT_VELOCITY_WORKER;
  set_core_affinity(tinfo->assigned_core_num);

  INFO.range_all = example.mac_grid.counts;
  INFO.range_x = INFO.range_y = INFO.range_z = INFO.range_re = INFO.range_all;
  INFO.range_re += TV_INT(2, 2, 2);
  INFO.range_x += TV_INT(2, 1, 1);
  INFO.range_y += TV_INT(1, 2, 1);
  INFO.range_z += TV_INT(1, 1, 2);
  TV_INT segment_start(1, 1, 1);
  pthread_mutex_lock(&INFO.mutex_fetcher);
  //int t = 0;
  while (true) {
    while (INFO.fetcher_stop || (INFO.task_recv_buffer->top == INFO.TASK_LIST_LENGTH)) {
      pthread_cond_wait(&INFO.cond_fetcher_go, &INFO.mutex_fetcher);
    }

    if (INFO.fetcher_refresh) {
      segment_start = TV_INT(1, 1, 1);
      INFO.fetcher_refresh = false;
      continue;
    }

    while (!INFO.fetcher_stop && (INFO.task_recv_buffer->top < INFO.TASK_LIST_LENGTH)) {
      // Calculate the next task of the new region.
      INFO.task_recv_buffer->task_content[INFO.task_recv_buffer->top++] = segment_start;
      // t = (t+1)%10;
      // if (t == 0) {
      //   pthread_cond_signal(&INFO.cond_fetcher_ready);
      //   pthread_mutex_unlock(&INFO.mutex_fetcher);
      //}
      segment_start(3) += INFO.segment_len;
      if (segment_start(3) >= INFO.range_re(3)) {
        segment_start(3) = 1;
        segment_start(2) += INFO.segment_len;
        if (segment_start(2) >= INFO.range_re(2)) {
          segment_start(2) = 1;
          segment_start(1) += INFO.segment_len;
          if (segment_start(1) >= INFO.range_re(1)) {
            segment_start(1) = 1;
            INFO.fetcher_stop = true;
          }
        }
      }
      //if (t == 0) pthread_mutex_lock(&INFO.mutex_fetcher);
    } // End while
    pthread_cond_signal(&INFO.cond_fetcher_ready);
  }
  return NULL;
}
*/
}
