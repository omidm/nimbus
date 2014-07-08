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

#define DEBUG_MODE

#include <iostream> // NOLINT
#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"
#include "shared/ldo_index.h"
#include "shared/ldo_index_rtree.h"

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
  printf("**Object - ID: %llu, Name: %s", obj->id(), obj->variable().c_str());  // NOLINT
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

const int kPnum = 54;
void AddLdo(nimbus::LdoIndex* index) {
  nimbus::GeometricRegion kRegions[kPnum];
  kRegions[0] = nimbus::GeometricRegion(1, 1, 1, 3, 3, 3);
  kRegions[1] = nimbus::GeometricRegion(1, 1, 4, 3, 3, 34);
  kRegions[2] = nimbus::GeometricRegion(1, 1, 38, 3, 3, 3);
  kRegions[3] = nimbus::GeometricRegion(1, 4, 1, 3, 14, 3);
  kRegions[4] = nimbus::GeometricRegion(1, 4, 4, 3, 14, 34);
  kRegions[5] = nimbus::GeometricRegion(1, 4, 38, 3, 14, 3);
  kRegions[6] = nimbus::GeometricRegion(1, 18, 1, 3, 3, 3);
  kRegions[7] = nimbus::GeometricRegion(1, 18, 4, 3, 3, 34);
  kRegions[8] = nimbus::GeometricRegion(1, 18, 38, 3, 3, 3);
  kRegions[9] = nimbus::GeometricRegion(1, 21, 1, 3, 3, 3);
  kRegions[10] = nimbus::GeometricRegion(1, 21, 4, 3, 3, 34);
  kRegions[11] = nimbus::GeometricRegion(1, 21, 38, 3, 3, 3);
  kRegions[12] = nimbus::GeometricRegion(1, 24, 1, 3, 14, 3);
  kRegions[13] = nimbus::GeometricRegion(1, 24, 4, 3, 14, 34);
  kRegions[14] = nimbus::GeometricRegion(1, 24, 38, 3, 14, 3);
  kRegions[15] = nimbus::GeometricRegion(1, 38, 1, 3, 3, 3);
  kRegions[16] = nimbus::GeometricRegion(1, 38, 4, 3, 3, 34);
  kRegions[17] = nimbus::GeometricRegion(1, 38, 38, 3, 3, 3);
  kRegions[18] = nimbus::GeometricRegion(4, 1, 1, 34, 3, 3);
  kRegions[19] = nimbus::GeometricRegion(4, 1, 4, 34, 3, 34);
  kRegions[20] = nimbus::GeometricRegion(4, 1, 38, 34, 3, 3);
  kRegions[21] = nimbus::GeometricRegion(4, 4, 1, 34, 14, 3);
  kRegions[22] = nimbus::GeometricRegion(4, 4, 4, 34, 14, 34);
  kRegions[23] = nimbus::GeometricRegion(4, 4, 38, 34, 14, 3);
  kRegions[24] = nimbus::GeometricRegion(4, 18, 1, 34, 3, 3);
  kRegions[25] = nimbus::GeometricRegion(4, 18, 4, 34, 3, 34);
  kRegions[26] = nimbus::GeometricRegion(4, 18, 38, 34, 3, 3);
  kRegions[27] = nimbus::GeometricRegion(4, 21, 1, 34, 3, 3);
  kRegions[28] = nimbus::GeometricRegion(4, 21, 4, 34, 3, 34);
  kRegions[29] = nimbus::GeometricRegion(4, 21, 38, 34, 3, 3);
  kRegions[30] = nimbus::GeometricRegion(4, 24, 1, 34, 14, 3);
  kRegions[31] = nimbus::GeometricRegion(4, 24, 4, 34, 14, 34);
  kRegions[32] = nimbus::GeometricRegion(4, 24, 38, 34, 14, 3);
  kRegions[33] = nimbus::GeometricRegion(4, 38, 1, 34, 3, 3);
  kRegions[34] = nimbus::GeometricRegion(4, 38, 4, 34, 3, 34);
  kRegions[35] = nimbus::GeometricRegion(4, 38, 38, 34, 3, 3);
  kRegions[36] = nimbus::GeometricRegion(38, 1, 1, 3, 3, 3);
  kRegions[37] = nimbus::GeometricRegion(38, 1, 4, 3, 3, 34);
  kRegions[38] = nimbus::GeometricRegion(38, 1, 38, 3, 3, 3);
  kRegions[39] = nimbus::GeometricRegion(38, 4, 1, 3, 14, 3);
  kRegions[40] = nimbus::GeometricRegion(38, 4, 4, 3, 14, 34);
  kRegions[41] = nimbus::GeometricRegion(38, 4, 38, 3, 14, 3);
  kRegions[42] = nimbus::GeometricRegion(38, 18, 1, 3, 3, 3);
  kRegions[43] = nimbus::GeometricRegion(38, 18, 4, 3, 3, 34);
  kRegions[44] = nimbus::GeometricRegion(38, 18, 38, 3, 3, 3);
  kRegions[45] = nimbus::GeometricRegion(38, 21, 1, 3, 3, 3);
  kRegions[46] = nimbus::GeometricRegion(38, 21, 4, 3, 3, 34);
  kRegions[47] = nimbus::GeometricRegion(38, 21, 38, 3, 3, 3);
  kRegions[48] = nimbus::GeometricRegion(38, 24, 1, 3, 14, 3);
  kRegions[49] = nimbus::GeometricRegion(38, 24, 4, 3, 14, 34);
  kRegions[50] = nimbus::GeometricRegion(38, 24, 38, 3, 14, 3);
  kRegions[51] = nimbus::GeometricRegion(38, 38, 1, 3, 3, 3);
  kRegions[52] = nimbus::GeometricRegion(38, 38, 4, 3, 3, 34);
  kRegions[53] = nimbus::GeometricRegion(38, 38, 38, 3, 3, 3);
  for (int i = 0; i < kPnum; ++i) {
    index->AddObject(
        new nimbus::LogicalDataObject(i, "velocity", &kRegions[i]));
  }
}

int main(int argc, char *argv[]) {
  nimbus::LdoIndex index;
  nimbus::LdoIndexRtree index_rtree;
  AddLdo(&index);
  AddLdo(&index_rtree);

  nimbus::GeometricRegion kRegY2W3Outer[2];
  nimbus::GeometricRegion kRegY2W3Central[2];

  // kRegY2W3Outer
  kRegY2W3Outer[0].Rebuild(-2, -2, -2, 46, 26, 46);
  kRegY2W3Outer[1].Rebuild(-2, 18, -2, 46, 26, 46);
  // kRegY2W3Central
  kRegY2W3Central[0].Rebuild(1, 1, 1, 40, 20, 40);
  kRegY2W3Central[1].Rebuild(1, 21, 1, 40, 20, 40);

  // Test operations
  nimbus::CLdoVector v1, v2;
  int count_v1, count_v2;

  count_v1 = index.IntersectingObjects("velocity", &kRegY2W3Outer[0], &v1);
  count_v2 = index_rtree.IntersectingObjects("velocity", &kRegY2W3Outer[0], &v2);
  printf("%d for ldo_index; %d for ldo_index_rtree\n", count_v1, count_v2);
  assert(count_v1 == count_v2);

  count_v1 = index.IntersectingObjects("velocity", &kRegY2W3Outer[1], &v1);
  count_v2 = index_rtree.IntersectingObjects("velocity", &kRegY2W3Outer[1], &v2);
  printf("%d for ldo_index; %d for ldo_index_rtree\n", count_v1, count_v2);
  assert(count_v1 == count_v2);

  count_v1 = index.IntersectingObjects("velocity", &kRegY2W3Central[0], &v1);
  count_v2 = index_rtree.IntersectingObjects("velocity", &kRegY2W3Central[0], &v2);
  printf("Result for v1:\n");
  printLdoVector(&v1);
  printf("Result for v2:\n");
  printLdoVector(&v2);
  printf("%d for ldo_index; %d for ldo_index_rtree\n", count_v1, count_v2);
  assert(count_v1 == count_v2);

  count_v1 = index.IntersectingObjects("velocity", &kRegY2W3Central[1], &v1);
  count_v2 = index_rtree.IntersectingObjects("velocity", &kRegY2W3Central[1], &v2);
  printf("%d for ldo_index; %d for ldo_index_rtree\n", count_v1, count_v2);
  assert(count_v1 == count_v2);

  return 0;
}
