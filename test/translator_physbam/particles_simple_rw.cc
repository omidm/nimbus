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
  * This file tests whether the PhysBAM translator correctly restores and
  * saves a face array consisting of two physical instances that perfectly
  * fit the region.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "data/physbam/translator_physbam.h"
#include "data/physbam/physbam_data.h"
#include "data/physbam/physbam_include.h"

#define R1_VALUE 3.1400
#define R2_VALUE 32.0

void printLdo(nimbus::LogicalDataObject* obj) {
  printf("**Object - ID: %lu, Name: %s", obj->id(), obj->variable().c_str());
  printf(" region: [%lu+%lu, %lu+%lu, %lu+%lu]\n", obj->region()->x(), obj->region()->dx(), obj->region()->y(), obj->region()->dy(), obj->region()->z(), obj->region()->dz());  // NOLINT
}

const int_dimension_t X = 10;
const int_dimension_t Y = 11;
const int_dimension_t Z = 12;
const int_dimension_t DX = 2;
const int_dimension_t DY = 2;
const int_dimension_t DZ = 21;
const int_dimension_t SIZE = DX * DY * DZ;
const int AVG_PARTICLES = 10;
const int TOTAL_PARTICLES = SIZE * AVG_PARTICLES * 5;
const int NUM_GHOST_CELLS = 3;

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
typedef PhysBAM::VECTOR<int, 3> TV_INT;
typedef PhysBAM::GRID<VECTOR_TYPE> Grid;
typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleLevelset;
typedef typename PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM<Grid> ParticleLevelsetEvolution;

int main(int argc, char *argv[]) {
  dbg_init();


  // 5 because a particle is a 5-tuple: x, y, z, radius, collision range

  int_dimension_t dimensions[] = {X, Y, Z, DX, DY, DZ};
  nimbus::GeometricRegion* region = new nimbus::GeometricRegion(dimensions);
  CPdiVector vec;
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

  vec.push_back(instance);

  Grid grid(TV_INT(), PhysBAM::RANGE<VECTOR_TYPE>::Unit_Box(), true);
  ParticleLevelsetEvolution evolution(grid, NUM_GHOST_CELLS);

  grid.Initialize(TV_INT(), PhysBAM::RANGE<VECTOR_TYPE>::Unit_Box(), true);
  // WARNING -- It seems broken that I have to set MAC_offset explicitly
  // for Initialize_Domain to work. I'm clearly missing an initialization step
  grid.MAC_offset = 0.5;
  evolution.Initialize_Domain(grid);
  evolution.particle_levelset.Set_Band_Width(6);
  evolution.Set_Time(0);
  evolution.Set_CFL_Number(0.9);
  evolution.Set_Number_Particles_Per_Cell(2 * AVG_PARTICLES);

  bool result;
  result = translator.ReadParticles(region, &vec, &evolution.particle_levelset, true);
  result = translator.ReadParticles(region, &vec, &evolution.particle_levelset, false);

  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    floats[i] = 1.0;
  }

  bool pass = true;
  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    if (floats[i] != floatSource[i]) {
      dbg(DBG_ERROR, "Value in physical instance 1 should be %f, it's %f.\n", floatSource[i], floats[i]);  //  NOLINT
      pass = false;
    }
  }

  if (pass) {
    printf("Passed all tests successfully.\n");
  } else {
    printf("Tests failed. Use dbg=error to observe errors.\n");
  }
}
