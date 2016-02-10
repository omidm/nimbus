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
  * Tests simple read and write operations.
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <cmath>
#include <vector>

#include "data/physbam/physbam_data.h"
#include "data/physbam/physbam_include.h"
#include "../translator_physbam_test.h"  // NOLINT

typedef nimbus::TranslatorPhysBAMTest<PhysBAM::VECTOR<float, 3> > Translator;
typedef typename Translator::ParticleInternal ParticleType;

const int_dimension_t X = 1;
const int_dimension_t Y = 1;
const int_dimension_t Z = 1;
const int_dimension_t DX = 1;
const int_dimension_t DY = 1;
const int_dimension_t DZ = 1;
const int_dimension_t SIZE = DX * DY * DZ;
const int AVG_PARTICLES = 10;
const int TOTAL_PARTICLES = SIZE * AVG_PARTICLES;
const int GRID_SCALE = 64;
const int NUM_GHOST_CELL = 2;
const typename Translator::scalar_t MIN_R = 0.2;
const typename Translator::scalar_t MAX_R = 0.5;

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

float getR() {
  double val = drand48();
  val *= static_cast<double>(MAX_R - MIN_R);
  val += static_cast<double>(MIN_R);
  return static_cast<float>(val);
}

// The goto clause prevents me from declaring variables. So switch.
void Error() {
  printf("Tests failed. Use dbg=error to observe errors.\n");
  exit(1);
}

int main(int argc, char *argv[]) {
  dbg_init();
  PhysBAM::LOG::Initialize_Logging();

  // Sets up test data.
  ParticleType* particles = new ParticleType[TOTAL_PARTICLES];
  ParticleType* particlesSource = new ParticleType[TOTAL_PARTICLES];
  for (int i = 0; i < TOTAL_PARTICLES; ++i) {
    float temp_x = getX();
    float temp_y = getY();
    float temp_z = getZ();
    particlesSource[i].index[0] = floor(temp_x);
    particlesSource[i].index[1] = floor(temp_y);
    particlesSource[i].index[2] = floor(temp_z);
    particlesSource[i].delta[0] = temp_x - floor(temp_x);
    particlesSource[i].delta[1] = temp_y - floor(temp_y);
    particlesSource[i].delta[2] = temp_z - floor(temp_z);
    particlesSource[i].radius = getR();
    particlesSource[i].quantized_collision_distance = 0xEFEF;
    particlesSource[i].id = 0x0EDEDEDE;
    particles[i] = particlesSource[i];
  }

  Translator translator;

  int_dimension_t dimensions[] = {X, Y, Z, DX, DY, DZ};
  nimbus::GeometricRegion* r = new nimbus::GeometricRegion(dimensions);
  nimbus::LogicalDataObject* ldo = new LogicalDataObject(1, "velocity", r);

  PhysBAMData* pd = new PhysBAMData();
  pd->set_buffer(
      reinterpret_cast<char*>(particles),
      TOTAL_PARTICLES * sizeof(ParticleType));
  PhysicalDataInstance* instance = new PhysicalDataInstance(1, ldo, pd, 0);
  PdiVector vec;
  vec.push_back(instance);

  PhysBAM::ARRAY<float, Translator::TV_INT> phi_array;
  phi_array.Resize(
      PhysBAM::RANGE<Translator::TV_INT>(
          1, GRID_SCALE,
          1, GRID_SCALE,
          1, GRID_SCALE).Thickened(NUM_GHOST_CELL));
  Translator::Grid grid(
      Translator::TV_INT(GRID_SCALE, GRID_SCALE, GRID_SCALE),
      PhysBAM::RANGE<Translator::TV>::Unit_Box(),
      true);
  Translator::ParticleContainer particle_container(grid, phi_array, 0);
  particle_container.Initialize_Particle_Levelset_Grid_Values();

  nimbus::GeometricRegion region(dimensions);
  ParticleType* result;
  bool pass = true;
  pass = translator.ReadParticles(&region, &vec, particle_container, true);

  if (!pass) {
    Error();
  }

  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    memset(&particles[i], 0xEF, sizeof(particles[i]));
  }

  pass = translator.WriteParticles(&region, &vec, particle_container, true);
  if (!pass) {
    Error();
  }

  result = reinterpret_cast<ParticleType*>(
      reinterpret_cast<PhysBAMData*>(instance->data())->buffer());

  if (reinterpret_cast<PhysBAMData*>(instance->data())->size() !=
      TOTAL_PARTICLES * sizeof(ParticleType)) {
    printf("The number of particles is incorrect!\n");
    Error();
  }

  std::vector<bool> checked(TOTAL_PARTICLES, false);
  for (int i = 0; i < TOTAL_PARTICLES; ++i) {
    // Check whether the ith particle in particlesSource exists in result.
    bool findit = false;
    // Go throught each particle in result.
    for (int j = 0; j < TOTAL_PARTICLES; ++j)
      if (!checked[j]) {
        if (memcmp(&particlesSource[i], &result[j], sizeof(result[j]))) {
          checked[j] = true;
          findit = true;
          break;
        }
      }
    if (!findit) {
      printf("The %dth Particle set is not found in the output.\n", i);
      Error();
    }
  }

  printf("Passed all tests successfully.\n");
  return 0;
}
