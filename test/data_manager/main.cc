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
  * This file tests the data manager and its memory use.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#define DEBUG_MODE
#define NOBJECTS 5
#define LOBJECTS 100
#include <iostream> // NOLINT
#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"
#include "scheduler/data_manager.h"
#include "shared/dbg.h"

using namespace nimbus;  // NOLINT

char relationship(GeometricRegion* r1, GeometricRegion* r2) {
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

void printTable(char* table) {
  printf("   A B C D E\n");
  printf("  +---------\n");
  for (int i = 0; i < NOBJECTS; i++) {
    printf("%c |", 'A' + i);
    for (int j = 0; j < NOBJECTS; j++) {
      printf("%c ", table[i*NOBJECTS + j]);
    }
    printf("\n");
  }
  printf("\n");
}

void printLdo(const LogicalDataObject* obj) {
  printf("**Object - ID: %llu, Name: %s", obj->id(), obj->variable().c_str());
  printf(" region: [%llu+%llu, %llu+%llu, %llu+%llu]\n", obj->region()->x(), obj->region()->dx(), obj->region()->y(), obj->region()->dy(), obj->region()->z(), obj->region()->dz());  // NOLINT
}

int main(int argc, char *argv[]) {
  dbg_init();
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

  GeometricRegion* regionA = new GeometricRegion(10, 11, 12, 22, 29, 33);
  dbg(DBG_MEMORY, "Allocated geo region 0x%x\n", regionA);
  LogicalDataObject ldoA(0, "velocity", regionA);

  GeometricRegion* regionB = new GeometricRegion(20, 11, 12, 22, 29, 33);
  dbg(DBG_MEMORY, "Allocated geo region 0x%x\n", regionB);
  LogicalDataObject ldoB(1, "velocity", regionB);

  GeometricRegion* regionC = new GeometricRegion(1, 11, 12, 1, 29, 33);
  dbg(DBG_MEMORY, "Allocated geo region 0x%x\n", regionC);
  LogicalDataObject ldoC(2, "velocity", regionC);

  GeometricRegion* regionD = new GeometricRegion(9, 15, 12, 1, 29, 33);
  dbg(DBG_MEMORY, "Allocated geo region 0x%x\n", regionD);
  LogicalDataObject ldoD(3, "velocity", regionD);

  GeometricRegion* regionE = new GeometricRegion(2, 9, 10, 50, 50, 66);
  dbg(DBG_MEMORY, "Allocated geo region 0x%x\n", regionE);
  LogicalDataObject ldoE(4, "velocity", regionE);

  LogicalDataObject* objects[NOBJECTS];
  objects[0] = &ldoA;
  objects[1] = &ldoB;
  objects[2] = &ldoC;
  objects[3] = &ldoD;
  objects[4] = &ldoE;

  char expected[] = {'C', 'I', 'D', 'A', 'I',
                     'I', 'C', 'D', 'D', 'I',
                     'D', 'D', 'C', 'D', 'A',
                     'A', 'D', 'D', 'C', 'I',
                     'C', 'C', 'A', 'C', 'C'};

  printf("\n");
  printf("** Checking overlap matrix **\n");
  printf("Expected output: \n");
  printTable(expected);

  char output[NOBJECTS * NOBJECTS];
  for (int i = 0; i < NOBJECTS; i++) {
    for (int j = 0; j < NOBJECTS; j++) {
      output[NOBJECTS * i + j] = relationship(objects[i]->region(), objects[j]->region());
    }
  }
  printf("Output:\n");
  printTable(output);

  DataManager mg;

  for (int i = 0; i < NOBJECTS; i++) {
    mg.AddLogicalObject(objects[i]->id(),
                        objects[i]->variable(),
                        *objects[i]->region());
  }

  CLdoVector results;

  int count = mg.FindLogicalObjects("velocity", &results);

  printf("** Printing %i velocity objects\n", count);
  CLdoVector::iterator it = results.begin();
  for (; it != results.end(); ++it) {
    const LogicalDataObject* obj = *it;
    printLdo(obj);
  }

  DataManager manager;
  //  LogicalDataObject* largeObjects[LOBJECTS];
  for (int i = 0; i < LOBJECTS; i++) {
    GeometricRegion region(lrand48(), lrand48(),
                           lrand48(), lrand48(),
                           lrand48(), lrand48());
    std::string str = "pressure";
    manager.AddLogicalObject(i, str, region);
    manager.RemoveLogicalObject(i);
  }
}
