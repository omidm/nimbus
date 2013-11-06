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
  * This file tests whether the Logical Data Object Index is working
  * properly.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#define DEBUG_MODE

#include <iostream> // NOLINT
#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"
#include "shared/ldo_index.h"

char relationship(nimbus::GeometricRegion* r1, nimbus::GeometricRegion* r2) {
  if (r1->Covers(r2)) {
    return 'C';
  } else if (r1->Intersects(r2)) {
    return 'I';
  } else if (r1->Adjacent(r2)) {
    return 'A';
  } else {
    return 'D';
  }
}

void printLdo(const nimbus::LogicalDataObject* obj) {
  printf("**Object - ID: %llu, Name: %s", obj->id(), obj->variable().c_str());
  printf(" region: [%llu+%llu, %llu+%llu, %llu+%llu]\n", obj->region()->x(), obj->region()->dx(), obj->region()->y(), obj->region()->dy(), obj->region()->z(), obj->region()->dz());  // NOLINT
}

void printLdoVector(nimbus::CLdoVector* v) {
  nimbus::CLdoVector::iterator iter = v->begin();
  while (iter != v->end()) {
    const nimbus::LogicalDataObject* obj = *iter;
    printLdo(obj);
    ++iter;
  }
}

int main(int argc, char *argv[]) {
  // Correct relationships:
  /*
       A  B  C  D  E
     A C  I  D  A  I
     B I  C  D  D  I
     C D  D  C  D  A
     D A  D  D  C  I
     E C  C  A  C  C

     Where C is covers, I is intersects, A is adjacent and
     D is disjoint.

     Row A, column E, value I means that A intersects E.
     Row E, column A, value C means that E covers A.
  */

  nimbus::GeometricRegion* regionA = new nimbus::GeometricRegion(10, 11, 12, 22, 29, 33);
  nimbus::LogicalDataObject ldoA(56, "velocity", regionA);

  nimbus::GeometricRegion* regionB = new nimbus::GeometricRegion(20, 11, 12, 22, 29, 33);
  nimbus::LogicalDataObject ldoB(12312, "velocity", regionB);

  nimbus::GeometricRegion* regionC = new nimbus::GeometricRegion(1, 11, 12, 1, 29, 33);
  nimbus::LogicalDataObject ldoC(37, "velocity", regionC);

  nimbus::GeometricRegion* regionD = new nimbus::GeometricRegion(9, 15, 12, 1, 29, 33);
  nimbus::LogicalDataObject ldoD(51, "velocity", regionD);

  nimbus::GeometricRegion* regionE = new nimbus::GeometricRegion(2, 9, 10, 50, 50, 66);
  nimbus::LogicalDataObject ldoE(43543821112342, "velocity", regionE);

  printf("Checking overlap matrix.\n");
  printf("Expected output: \n");
  printf("   A  B  C  D  E\n");
  printf("  +-------------\n");
  printf("A |C  I  D  A  I\n");
  printf("B |I  C  D  D  I\n");
  printf("C |D  D  C  D  A\n");
  printf("D |A  D  D  C  I\n");
  printf("E |C  C  A  C  C\n");


  nimbus::LogicalDataObject* objects[5];

  objects[0] = &ldoA;
  objects[1] = &ldoB;
  objects[2] = &ldoC;
  objects[3] = &ldoD;
  objects[4] = &ldoE;

  nimbus::LdoIndex index;

  for (int i = 0; i < 5; i++) {
    index.AddObject(objects[i]);
    printf("%c - ", 'A' + i);
    printLdo(objects[i]);
  }

  // Test operations
  nimbus::CLdoVector v1;
  printf("Objects intersecting A. Should be A, B, and E.\n");
  std::string velocity = "velocity";
  std::string pressure = "pressure";
  index.IntersectingObjects(velocity,
                            ldoA.region(),
                            &v1);
  printLdoVector(&v1);

  nimbus::CLdoVector v2;
  printf("Testing variable map. Should have no output.\n");
  index.IntersectingObjects(pressure,
                            ldoA.region(),
                            &v2);
  printLdoVector(&v2);

  nimbus::CLdoVector v3;
  printf("Objects E covers. Should be A, B, D, E.\n");
  index.CoveredObjects(velocity,
                       ldoE.region(),
                       &v3);
  printLdoVector(&v3);
}
