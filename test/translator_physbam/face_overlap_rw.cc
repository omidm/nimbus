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

#define R1_VALUE 3.1400
#define R2_VALUE 32.0

void printLdo(nimbus::LogicalDataObject* obj) {
  printf("**Object - ID: %lu, Name: %s", obj->id(), obj->variable().c_str());
  printf(" region: [%lu+%lu, %lu+%lu, %lu+%lu]\n", obj->region()->x(), obj->region()->dx(), obj->region()->y(), obj->region()->dy(), obj->region()->z(), obj->region()->dz());  // NOLINT
}

int main(int argc, char *argv[]) {
  dbg_init();

  const int_dimension_t X = 10;
  const int_dimension_t Y = 10;
  const int_dimension_t Z = 10;
  const int_dimension_t DX = 10;
  const int_dimension_t DY = 1;
  const int_dimension_t DZ = 1;

  const int_dimension_t X1 = 8;
  const int_dimension_t Y1 = 9;
  const int_dimension_t Z1 = 10;
  const int_dimension_t DX1 = 6;
  const int_dimension_t DY1 = 2;
  const int_dimension_t DZ1 = 2;
  const int_dimension_t SIZE1 = DX1 * DY1 * DZ1;

  const int_dimension_t X2 = 14;
  const int_dimension_t Y2 = 10;
  const int_dimension_t Z2 = 10;
  const int_dimension_t DX2 = 6;
  const int_dimension_t DY2 = 2;
  const int_dimension_t DZ2 = 1;
  const int_dimension_t SIZE2 = DX2 * DY2 * DZ2;


  int_dimension_t dimensions[] = {X, Y, Z, DX, DY, DZ};
  nimbus::GeometricRegion* region = new nimbus::GeometricRegion(dimensions);
  CPdiVector vec;
  TranslatorPhysBAM<PhysBAM::VECTOR<float, 3> > translator;

  PhysBAM::ARRAY<float, PhysBAM::FACE_INDEX<3> >* result; // NOLINT

  int_dimension_t dimensions1[] = {X1, Y1, Z1, DX1, DY1, DZ1};
  nimbus::GeometricRegion* r1 = new nimbus::GeometricRegion(dimensions1);
  nimbus::LogicalDataObject* ldo1 = new LogicalDataObject(1, "velocity", r1);
  float* f1 = new float[3 * SIZE1];
  float* f1source = new float[3 * SIZE1];
  for (int i = 0; i < 3 * SIZE1; i++) {
    f1source[i] = R1_VALUE;
    f1[i] = f1source[i];
  }
  PhysBAMData* pd1 = new PhysBAMData();
  pd1->set_buffer((char*)f1, 3 * SIZE1 * sizeof(float));  // NOLINT
  PhysicalDataInstance* i1 = new PhysicalDataInstance(1, ldo1, pd1, 0);


  int_dimension_t dimensions2[] = {X2, Y2, Z2, DX2, DY2, DZ2};
  nimbus::GeometricRegion* r2 = new nimbus::GeometricRegion(dimensions2);
  nimbus::LogicalDataObject* ldo2 = new LogicalDataObject(2, "velocity", r2);
  float* f2 = new float[3 * SIZE2];
  float* f2source = new float[3 * SIZE2];
  for (int i = 0; i < 3 * SIZE2; i++) {
    f2source[i] = R2_VALUE;
    f2[i] = f2source[i];
  }
  PhysBAMData* pd2 = new PhysBAMData();
  pd2->set_buffer((char*)f2, 3 * SIZE2 * sizeof(float));  // NOLINT
  PhysicalDataInstance* i2 = new PhysicalDataInstance(2, ldo2, pd2, 0);

  vec.push_back(i1);
  vec.push_back(i2);

  result = translator.ReadFaceArray(region, &vec);

  for (int i = 0; i < 3 * SIZE1; i++) {
    f1[i] = 1.0;
  }
  for (int i = 0; i < 3 * SIZE2; i++) {
    f2[i] = 1.0;
  }

  bool pass = translator.WriteFaceArray(region, &vec, result);

  for (int z = Z; z < Z + DZ; z++) {
    for (int y = Y; y < Y + DY; y++) {
      for (int x = X; x < X + DX; x++) {
        if (z >= Z1 && z < (Z1 + DZ1) &&
            y >= Y1 && y < (Y1 + DY1) &&
            x >= X1 && z < (X1 + DZ1)) {
          int x_offset = x - X1;
          int y_offset = y - Y1;
          int z_offset = z - Z1;
          int offset = z_offset * (DX1 * DY1) + (y_offset * DX1) + x_offset;
          if (f1[offset] != f1source[offset]) {
            dbg(DBG_ERROR, "Value %i in PDI 1 should be %f, it's %f.\n",
                offset, f1source[offset], f1[offset]);
            pass = false;
          }
        }
        if (z >= Z2 && z < (Z2 + DZ2) &&
            y >= Y2 && y < (Y2 + DY2) &&
            x >= X2 && z < (X2 + DZ2)) {
          int x_offset = x - X2;
          int y_offset = y - Y2;
          int z_offset = z - Z2;
          int offset = z_offset * (DX2 * DY2) + (y_offset * DX2) + x_offset;
          if (f2[offset] != f2source[offset]) {
            dbg(DBG_ERROR, "Value %i in PDI 2 should be %f, it's %f.\n",
                offset, f2source[offset], f2[offset]);
            pass = false;
          }
        }
      }
    }
  }

  if (pass) {
    printf("Passed all tests successfully.\n");
  } else {
    printf("Tests failed. Use dbg=error to observe errors.\n");
  }
}
