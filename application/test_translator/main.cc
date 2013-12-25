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
  *
  */

#include <algorithm>

// Temporery import for nimbus internal data structures.
#include "data/physbam/physbam_include.h"
#include "data/physbam/physbam_data.h"

#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "worker/physical_data_instance.h"
#include "worker/worker.h"

// Temporery import for nimbus internal data structures. End.

#include "data/physbam/translator_physbam.h"
#include "data/physbam/physbam_data.h"
#include "data/physbam/physbam_include.h"

#define R1_VALUE 3.1400
#define R2_VALUE 32.0

/*
void printLdo(nimbus::LogicalDataObject* obj) {
  printf("**Object - ID: %lu, Name: %s", obj->id(), obj->variable().c_str());
  printf(" region: [%lu+%lu, %lu+%lu, %lu+%lu]\n", obj->region()->x(), obj->region()->dx(), obj->region()->y(), obj->region()->dy(), obj->region()->z(), obj->region()->dz());  // NOLINT
}
*/

const int_dimension_t X = 1;
const int_dimension_t Y = 1;
const int_dimension_t Z = 1;
const int_dimension_t DX = 2;
const int_dimension_t DY = 2;
const int_dimension_t DZ = 2;
const int_dimension_t SIZE = DX * DY * DZ;
const int AVG_PARTICLES = 10;
const int TOTAL_PARTICLES = SIZE * AVG_PARTICLES * 5;
const int NUM_GHOST_CELLS = 2;
const int GRID_SCALE = 16;

float getX() {
  double val = drand48();
  val *= static_cast<double>(DX);
  val += static_cast<double>(X);
  return static_cast<float>(val);
}

float getY() {
  double val = drand48();
  val *= static_cast<double>(DY);
  val += static_cast<double>(Y);
  return static_cast<float>(val);
}

float getZ() {
  double val = drand48();
  val *= static_cast<double>(DZ);
  val += static_cast<double>(Z);
  return static_cast<float>(val);
}

typedef PhysBAM::VECTOR<float, 3> VECTOR_TYPE;
typedef VECTOR_TYPE TV;
typedef typename TV::SCALAR scalar_t;
typedef PhysBAM::VECTOR<int, 3> TV_INT;
typedef PhysBAM::GRID<VECTOR_TYPE> Grid;
typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleProfile;
typedef typename PhysBAM::PARTICLE_LEVELSET_UNIFORM<Grid> ParticleContainer;
typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> ParticlesUnit;
typedef typename PhysBAM::ARRAY<ParticlesUnit*, TV_INT> ParticlesArray;

template <typename GRID_T>
void PrintGrid(GRID_T& grid) {
  printf("Start\n");
  printf("Counts %d %d %d\n", grid.counts(1),
                        grid.counts(2),
                        grid.counts(3));
  printf("Min %1.0f %1.0f %1.0f\n", (float)grid.domain.min_corner(1),
                        (float)grid.domain.min_corner(2),
                        (float)grid.domain.min_corner(3));
  printf("Max %1.0f %1.0f %1.0f\n", grid.domain.max_corner(1),
                        (float)grid.domain.max_corner(2),
                        (float)grid.domain.max_corner(3));
}

bool ReadParticles(nimbus::GeometricRegion* region,
                   nimbus::CPdiVector* instances,
                   ParticleContainer& particle_container,
                   bool positive) {
  using namespace nimbus;
  ParticlesArray* particles;
  if (positive) {
    particles = &particle_container.positive_particles;
  } else {
    particles = &particle_container.negative_particles;
  }
  for (int z = region->z(); z < region->z()+region->dz(); z++)
    for (int y = region->y(); y < region->y()+region->dy(); y++)
      for (int x = region->x(); x < region->x()+region->dx(); x++) {
        printf("Access %d %d %d\n", x, y, z);
        printf(" %1.0f %1.0f %1.0f\n", (float)particles->domain.min_corner(1),
                              (float)particles->domain.min_corner(2),
                              (float)particles->domain.min_corner(3));
        printf(" %1.0f %1.0f %1.0f\n", (float)particles->domain.max_corner(1),
                              (float)particles->domain.max_corner(2),
                              (float)particles->domain.max_corner(3));
        TV_INT block_index(x, y, z);
        if (!(*particles)(block_index)) {
          (*particles)(block_index) = particle_container.Allocate_Particles(
              particle_container.template_particles);
        }
      }

  assert(instances != NULL);

  CPdiVector::iterator iter = instances->begin();
  for (; iter != instances->end(); ++iter) {
    const PhysicalDataInstance* instance = *iter;
    PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
    scalar_t* buffer = reinterpret_cast<scalar_t*>(data->buffer());

    printf("%d\n", (int) data->size());
    for (int i = 0; i < (int) data->size() / (int) sizeof(float); i+= 5) {
      // This instance may have particles inside the region. So
      // we have to iterate over them and insert them accordingly.
      // Particles are stored as (x,y,z) triples in data array
      scalar_t x = buffer[i];
      scalar_t y = buffer[i + 1];
      scalar_t z = buffer[i + 2];
      scalar_t radius = buffer[i + 3];

      uint16_t collision_distance = (uint16_t)buffer[i + 4];  // NOLINT

      VECTOR_TYPE position;
      position.x = x;
      position.y = y;
      position.z = z;
      int_dimension_t xi = (int_dimension_t)x;
      int_dimension_t yi = (int_dimension_t)y;
      int_dimension_t zi = (int_dimension_t)z;

      // If particle is within region, copy it to particles
      if (xi >= region->x() &&
          xi < (region->x() + region->dx()) &&
          yi >= region->y() &&
          yi < (region->y() + region->dy()) &&
          zi >= region->z() &&
          zi < (region->z() + region->dz())) {
        ParticlesUnit* cellParticles = (*particles)(TV_INT(xi, yi, zi));

        // Note that Add_Particle traverses a linked list of particle
        // buckets, so it's O(N^2) time. Blech.
        int index = particle_container.Add_Particle(cellParticles);
        cellParticles->quantized_collision_distance(index) =
          collision_distance;
        cellParticles->X(index) = position;
        cellParticles->radius(index) = radius;
        printf("Write %f %f %f to %d %d %d\n", x, y, z,
            (int)xi, (int)yi,(int) zi);
      }
    }
  }
  return true;
}

bool WriteParticles(nimbus::GeometricRegion* region,
                    nimbus::CPdiVector* instances,
                    ParticleContainer& particle_container,
                    bool positive) {
  using namespace nimbus;
  CPdiVector::iterator iter = instances->begin();
  for (; iter != instances->end(); ++iter) {
    const PhysicalDataInstance* instance = *iter;
    PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
    data->ClearTempBuffer();
  }

  for (int z = 1; z <= region->dz(); z++) {
    for (int y = 1; y <= region->dy(); y++) {
      for (int x = 1; x <= region->dx(); x++) {
        ParticlesArray* arrayPtr;
        if (positive) {
          arrayPtr = &particle_container.positive_particles;
        } else {
          arrayPtr = &particle_container.negative_particles;
        }
        // printf("Pointer to positive particles is %p\n", arrayPtr);
        ParticlesUnit* particles = (*arrayPtr)(TV_INT(x, y, z));
        while (particles) {
          for (int i = 1; i <= particles->array_collection->Size(); i++) {
            VECTOR_TYPE particle = particles->X(i);
            scalar_t x = particle.x;
            scalar_t y = particle.y;
            scalar_t z = particle.z;
            double xi = x;
            double yi = y;
            double zi = z;
            // If it's inside the region,
            if (xi >= region->x() &&
                xi <= (region->x() + region->dx()) &&
                yi >= region->y() &&
                yi <= (region->y() + region->dy()) &&
                zi >= region->z() &&
                zi <= (region->z() + region->dz())) {
              CPdiVector::iterator iter = instances->begin();
              // Iterate across instances, checking each one
              for (; iter != instances->end(); ++iter) {
                const PhysicalDataInstance* instance = *iter;
                GeometricRegion* instanceRegion = instance->region();
                PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());

                // If it's inside the region of the physical data instance
                if (xi >= instanceRegion->x() &&
                    xi <= (instanceRegion->x() + instanceRegion->dx()) &&
                    yi >= instanceRegion->y() &&
                    yi <= (instanceRegion->y() + instanceRegion->dy()) &&
                    zi >= instanceRegion->z() &&
                    zi <= (instanceRegion->z() + instanceRegion->dz())) {
                  printf("Read %f %f %f\n", x, y, z);
                  scalar_t particleBuffer[5];
                  particleBuffer[0] = particle.x;
                  particleBuffer[1] = particle.y;
                  particleBuffer[2] = particle.z;
                  particleBuffer[3] = particles->radius(i);
                  particleBuffer[4] = particles->quantized_collision_distance(i);
                  data->AddToTempBuffer(reinterpret_cast<char*>(particleBuffer),
                                        sizeof(scalar_t)*5);
                }
              }
            }
          }
          particles = particles->next;
        }
      }
    }
  }

  // Now that we've copied particles into the temporary data buffers,
  // commit the results.
  iter = instances->begin();
  for (; iter != instances->end(); ++iter) {
    const PhysicalDataInstance* instance = *iter;
    PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
    data->CommitTempBuffer();
  }
  return true;
}

int main(int argc, char *argv[]) {
  dbg_init();
  PhysBAM::LOG::Initialize_Logging();

  // 5 because a particle is a 5-tuple: x, y, z, radius, collision range

  int_dimension_t dimensions[] = {X, Y, Z, DX, DY, DZ};
  nimbus::GeometricRegion* region = new nimbus::GeometricRegion(dimensions);
  CPdiVector* vec = new CPdiVector();
  TranslatorPhysBAM<PhysBAM::VECTOR<float, 3> > translator;

  nimbus::GeometricRegion* r = new nimbus::GeometricRegion(dimensions);
  nimbus::LogicalDataObject* ldo = new LogicalDataObject(1, "velocity", r);
  float* floats = new float[TOTAL_PARTICLES];
  float* floatSource = new float[TOTAL_PARTICLES];
  for (int i = 0; i < TOTAL_PARTICLES; i+=5) {
    floatSource[i + 0] = floats[i + 0] = getX();
    floatSource[i + 1] = floats[i + 1] = getY();
    floatSource[i + 2] = floats[i + 2] = getZ();
    floatSource[i + 3] = floats[i + 3] = 1.0;
    floatSource[i + 4] = floats[i + 4] = 1.0;
  }
  PhysBAMData* pd = new PhysBAMData();
  pd->set_buffer((char*)floats, TOTAL_PARTICLES * sizeof(float));  // NOLINT
  PhysicalDataInstance* instance = new PhysicalDataInstance(1, ldo, pd, 0);

  vec->push_back(instance);

  PhysBAM::ARRAY<float, TV_INT> phi_array;
  phi_array.Resize(TV_INT(DX, DY, DZ));

  PhysBAM::RANGE<VECTOR_TYPE> range_input(X, X+DX-1, Y, Y+DY-1, Z, Z+DZ-1);

  Grid grid(TV_INT(DX, DY, DZ), range_input, true);

  /*
  if (grid.Is_MAC_Grid()) {
    printf("Is MAC grid.\n");
  } else {
    printf("Is not MAC grid.\n");
  }
  */
  // Critical Warning: Still cannot understand what is the memory fault.
  // "grid" will be changed in an unexpect way.  -quhang
  PrintGrid(grid);
  // Grid regular_grid = grid.Get_Regular_Grid();
  // PrintGrid(regular_grid);
  ParticleContainer particle_container(grid, phi_array, 0);
  particle_container.Initialize_Particle_Levelset_Grid_Values();
  // PrintGrid(regular_grid);
  // PrintGrid(particle_container.levelset.grid);
  float* result;
  bool pass = true;
  pass = ReadParticles(region, vec, particle_container, true);

  if (!pass) {
    printf("Failed to read particles.\n");
    goto error;
  }

  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    floats[i] = 1.0;
  }

  pass = WriteParticles(region, vec, particle_container, true);
  if (!pass) {
    printf("Failed to write particles.\n");
    goto error;
  }

  /*
  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    if (floats[i] != floatSource[i]) {
      dbg(DBG_ERROR, "Value in physical instance 1 should be %f, it's %f.\n", floatSource[i], floats[i]);  //  NOLINT
      pass = false;
    }
  }
  */
  // Cannot compare to float, which is deleted by PhysBAM data. --quhang
  // Assume particle values are different.
  result = (float*) ((PhysBAMData*)(instance->data()))->buffer();
  for (int i = 0; i < TOTAL_PARTICLES / 5; ++i) {
    bool any = false;
    printf("Search for %f\n", floatSource[i*5]);
    for (int j = 0; j < TOTAL_PARTICLES / 5; ++j) {
      bool flag = true;
      printf("Get %f\n", result[j*5]);
      for (int k = 0; k < 5; ++k)
        if (floatSource[i*5+k] != result[j*5+k]) {
          flag = false;
          break;
        }
      if (flag) {
        any = true;
        break;
      }
    }
    printf("coming to %d\n", i);
    if (!any) {
      pass = false;
      goto error;
    }
  }

 error:
  if (pass) {
    printf("Passed all tests successfully.\n");
  } else {
    printf("Tests failed. Use dbg=error to observe errors.\n");
  }
}
