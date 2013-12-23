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
const int_dimension_t DZ = 4;
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
typedef PhysBAM::VECTOR<int, 3> TV_INT;
typedef PhysBAM::GRID<VECTOR_TYPE> Grid;
typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleLevelset;
typedef typename PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM<Grid> ParticleLevelsetEvolution;
typedef typename PhysBAM::GEOMETRY_BOUNDARY_POLICY<Grid>::BOUNDARY_PHI_WATER BoundaryPhiWater;
typedef typename TranslatorPhysBAM<VECTOR_TYPE>::FaceArray FaceArray;

template <class TV>
class CALLBACKS:public PhysBAM::LEVELSET_CALLBACKS<PhysBAM::GRID<TV> > {
  typedef typename TV::SCALAR T;
  void Adjust_Particle_For_Domain_Boundaries(PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV>& particles,
                                             const int index,TV& V,
                                             const PhysBAM::PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
                                             const T dt,
                                             const T time) {} // NOLINT
  void Get_Levelset_Velocity(const PhysBAM::GRID<TV>& grid,
                             PhysBAM::LEVELSET_MULTIPLE_UNIFORM<PhysBAM::GRID<TV> >& levelset_multiple,
                             PhysBAM::ARRAY<T, PhysBAM::FACE_INDEX<TV::dimension> >& V_levelset,
                             const T time) const PHYSBAM_OVERRIDE {} // NOLINT
};

typedef CALLBACKS<VECTOR_TYPE> Callbacks;
int main(int argc, char *argv[]) {
  dbg_init();
  
  PhysBAM::PARSE_ARGS parse_args;
  parse_args.Add_Integer_Argument("-restart",0,"restart frame");
  parse_args.Add_Integer_Argument("-scale",128,"fine scale grid resolution");
  parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
  parse_args.Add_Integer_Argument("-e",100,"last frame");
  parse_args.Add_Integer_Argument("-refine",1,"refine levels");
  parse_args.Add_Integer_Argument("-threads",1,"number of threads");
  parse_args.Add_Double_Argument("-cfl",1,"cfl number");
  PhysBAM::LOG::Initialize_Logging(false,false,1<<30,true,parse_args.Get_Integer_Value("-threads"));  // NOLINT

  Callbacks* callbacks = new Callbacks();
  PhysBAM::RANGE<PhysBAM::VECTOR<int, 3> > range(0, GRID_SCALE,
                                                 0, GRID_SCALE,
                                                 0, GRID_SCALE);

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

  Grid* grid = new Grid(TV_INT(),
                        PhysBAM::RANGE<VECTOR_TYPE>::Unit_Box(),
                        true);

  grid->Initialize(TV_INT::All_Ones_Vector() * GRID_SCALE,
                  PhysBAM::RANGE<VECTOR_TYPE>(VECTOR_TYPE(),
                                              VECTOR_TYPE::All_Ones_Vector()),
                  true);

  ParticleLevelsetEvolution* evolution = new
    ParticleLevelsetEvolution(*grid, NUM_GHOST_CELLS);

  printf("grid: %p\n", grid);
  printf("evolution: %p\n", evolution);
  printf("evolution.particle_levelset: %p\n", &evolution->particle_levelset);
  printf("evolution.particle_levelset.levelset: %p\n", &evolution->particle_levelset.levelset);
  printf("evolution.particle_levelset.levelset.grid: %p\n", &evolution->particle_levelset.levelset.grid);  // NOLINT

  evolution->Initialize_Domain(*grid);
  evolution->particle_levelset.Set_Band_Width(6);
  evolution->Set_Time(0);
  evolution->Set_CFL_Number(0.9);
  evolution->Set_Number_Particles_Per_Cell(2 * AVG_PARTICLES);
  evolution->Use_Particle_Levelset(true);
  evolution->Set_Levelset_Callbacks(*callbacks);
  
  FaceArray* faceVelocities = new FaceArray();
  faceVelocities->Resize(range);

  BoundaryPhiWater* phiWaterBoundary = new BoundaryPhiWater();
  phiWaterBoundary->Set_Velocity_Pointer(*faceVelocities);
  // phiWaterBoundary.Set_Constant_Extrapolation(domainOpenBoundaries);

  evolution->particle_levelset.levelset.Set_Custom_Boundary(*phiWaterBoundary);
  // evolution->particle_levelset.Use_Removed_Negative_Particles();

  // evolution->particle_levelset.Store_Unique_Particle_Id();
  evolution->Use_Particle_Levelset(true);
  // evolution.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);  // NOLINT
  // evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.incompressible.valid_mask); // NOLINT

  // evolution->particle_levelset.Set_Collision_Distance_Factors(.1,1);
  // evolution->particle_levelset.Use_Removed_Positive_Particles();

  evolution->Set_Seed(2606);
  evolution->Seed_Particles(0);

  bool pass = true;
  pass = translator.ReadParticles(region, vec, evolution->particle_levelset, true);
  if (!pass) {
    printf("Failed to read particles.\n");
    goto error;
  }

  //  result = translator.ReadParticles(region, &vec, evolution.particle_levelset, false);

  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    floats[i] = 1.0;
  }

  pass = translator.WriteParticles(region, vec, evolution->particle_levelset, true);
  if (!pass) {
    printf("Failed to write particles.\n");
    goto error;
  }

  for (int i = 0; i < TOTAL_PARTICLES; i++) {
    if (floats[i] != floatSource[i]) {
      dbg(DBG_ERROR, "Value in physical instance 1 should be %f, it's %f.\n", floatSource[i], floats[i]);  //  NOLINT
      pass = false;
    }
  }

 error:
  if (pass) {
    printf("Passed all tests successfully.\n");
  } else {
    printf("Tests failed. Use dbg=error to observe errors.\n");
  }
}
