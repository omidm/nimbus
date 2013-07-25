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
#include "water-example.h"
#include "water-driver.h"

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
}  // namespace ADVECT_VELOCITY_NS
