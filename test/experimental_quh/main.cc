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
  * Experimental code to test the idea of job query. 
  * Not for testing or interface design.
  *
  * Author: Hang Qu <quhang@stanford.edu> 
  */

#include <string>
#include "shared/geometric_region.h"
#include "test/experimental_quh/query_test_engine.h"
#include "test/experimental_quh/variable_region_set.h"

// Basic test.
void test_1() {
  using nimbus::GeometricRegion;
  using nimbus::VariableRegionSet;
  GeometricRegion left(1, 1, 1, 64, 64, 1);
  GeometricRegion right(65, 1, 1, 64, 64, 1);
  GeometricRegion mid(2, 10, 1, 62, 10, 1);
  std::cout << left.toString() << std::endl;
  std::cout << right.toString() << std::endl;
  VariableRegionSet temp(
      VariableRegionSet("V", left) + VariableRegionSet("V", right));
  temp.DebugPrint();
  temp.IntersectsTest(VariableRegionSet("V", mid));
  temp.DebugPrint();
  temp.IntersectsAndDelete(VariableRegionSet("V", mid));
  temp.DebugPrint();
}

// The example shown in the group meeting.
void test2() {
  using nimbus::QueryTestEngine;
  using nimbus::QueryResult;
  using nimbus::GeometricRegion;
  using nimbus::VariableRegionSet;
  QueryTestEngine query_test_engine;
  QueryResult buffer[3];
  GeometricRegion left(1, 1, 1, 64, 64, 1);
  GeometricRegion mid(65, 1, 1, 64, 64, 1);
  GeometricRegion right(129, 1, 1, 64, 64, 1);
  GeometricRegion left_ghost(1, 1, 1, 65, 64, 1);
  GeometricRegion mid_ghost(64, 1, 1, 66, 64, 1);
  GeometricRegion right_ghost(128, 1, 1, 65, 64, 1);
  GeometricRegion core[3];
  core[0] = left;
  core[1] = mid;
  core[2] = right;
  GeometricRegion ghost[3];
  ghost[0] = left_ghost;
  ghost[1] = mid_ghost;
  ghost[2] = right_ghost;
  for (int iter = 0; iter < 5; iter++) {
    for (int i = 0; i < 3; ++i) {
      query_test_engine.Query(
          VariableRegionSet("V", ghost[i]),
          VariableRegionSet("P", core[i]),
          &buffer[i]);
    }
    for (int i = 0; i < 3; ++i) {
      query_test_engine.Add("VtoP", i, buffer[i]);
    }
    for (int i = 0; i < 3; ++i) {
      query_test_engine.Query(
          VariableRegionSet("P", ghost[i]),
          VariableRegionSet("V", core[i]),
          &buffer[i]);
    }
    for (int i = 0; i < 3; ++i) {
      query_test_engine.Add("PtoV", i, buffer[i]);
    }
  }
  query_test_engine.PruneResult();
  query_test_engine.WriteOut("output.dot");
}

// All the jobs in inner loop of projection.
void test3() {
  using nimbus::QueryTestEngine;
  using nimbus::QueryResult;
  using nimbus::GeometricRegion;
  using nimbus::VariableRegionSet;
  QueryTestEngine query_test_engine;
  QueryResult buffer[3];
  GeometricRegion left(1, 1, 1, 64, 64, 1);
  GeometricRegion mid(65, 1, 1, 64, 64, 1);
  GeometricRegion right(129, 1, 1, 64, 64, 1);
  GeometricRegion left_ghost(1, 1, 1, 65, 64, 1);
  GeometricRegion mid_ghost(64, 1, 1, 66, 64, 1);
  GeometricRegion right_ghost(128, 1, 1, 65, 64, 1);
  GeometricRegion whole(1, 1, 1, 192, 64, 1);
  GeometricRegion core[3];
  core[0] = left;
  core[1] = mid;
  core[2] = right;
  GeometricRegion ghost[3];
  ghost[0] = left_ghost;
  ghost[1] = mid_ghost;
  ghost[2] = right_ghost;
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Query(
        VariableRegionSet("r", core[i]),
        VariableRegionSet("z", core[i])
        + VariableRegionSet("pho_p", core[i])
        + VariableRegionSet("temp", core[i]),
        &buffer[i]);
  }
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Add("Calc rho:", i, buffer[i]);
  }
  query_test_engine.Query(
      VariableRegionSet("pho_p", whole),
      VariableRegionSet("pho", whole),
      &buffer[0]);
  query_test_engine.Add("Collect rho:", 9, buffer[0]);
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Query(
        VariableRegionSet("z", core[i])
        + VariableRegionSet("pho", whole)
        + VariableRegionSet("s", core[i]),
        VariableRegionSet("s", core[i]),
        &buffer[i]);
  }
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Add("Update vec:" , i, buffer[i]);
  }
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Query(
        VariableRegionSet("s", ghost[i]),
        VariableRegionSet("temp", core[i])
        + VariableRegionSet("tao_p", core[i]),
        &buffer[i]);
  }
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Add("Calc tao:" , i, buffer[i]);
  }
  query_test_engine.Query(
      VariableRegionSet("tao_p", whole),
      VariableRegionSet("tao", whole),
      &buffer[0]);
  query_test_engine.Add("Collect tao:", 9, buffer[0]);
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Query(
        VariableRegionSet("r", core[i])
        + VariableRegionSet("rho", whole)
        + VariableRegionSet("s", core[i])
        + VariableRegionSet("temp", core[i])
        + VariableRegionSet("x", core[i])
        + VariableRegionSet("tao", whole),
        VariableRegionSet("r", core[i])
        + VariableRegionSet("x", core[i])
        + VariableRegionSet("norm_p", core[i]),
        &buffer[i]);
  }
  for (int i = 0; i < 3; ++i) {
    query_test_engine.Add("Calc norm:" , i, buffer[i]);
  }
  query_test_engine.Query(
      VariableRegionSet("norm_p", whole),
      VariableRegionSet("norm", whole),
      &buffer[0]);
  query_test_engine.Add("Collect norm:", 9, buffer[0]);

  query_test_engine.PruneResult();
  query_test_engine.WriteOut("output.dot");
}

int main() {
  // test1();
  // test2();
  test3();
  return 0;
}
