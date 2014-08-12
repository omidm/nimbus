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
  * A test for building the data access pattern data structure to measure the
  * latency it adds.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <stdlib.h>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <iostream>  // NOLINT
#include <sstream>  // NOLINT
#include <string>
#include <map>
#include <vector>
#include <set>
#include "shared/nimbus.h"
#include "shared/log.h"
#include "shared/idset.h"
#include "shared/geometric_region.h"
#include "scheduler/region_map_entry.h"

using namespace nimbus; // NOLINT

int main(int argc, const char *argv[]) {
  nimbus_initialize();
  Log log;

  RegionMapEntry rme;

  GeometricRegion r1(0, 0, 0, 5, 5, 5);
  GeometricRegion r2(2, 2, 2, 5, 5, 5);

  rme.AddRegion(&r1);
  std::cout << rme.Print() << std::endl;

  rme.AddRegion(&r2);
  std::cout << rme.Print() << std::endl;

  RegionMapEntry::RegionList result;
  RegionMapEntry::RemoveIntersect(&r2, &r1, &result);
  RegionMapEntry::RegionListIter iter = result.begin();
  for (; iter != result.end(); ++iter) {
    std::cout << iter->ToNetworkData() << std::endl;
  }
  std::cout << std::endl;

  rme.RemoveRegion(&r1);
  std::cout << rme.Print() << std::endl;

  rme.RemoveRegion(&r2);
  std::cout << rme.Print() << std::endl;

  rme.AddRegion(&r1);
  std::cout << rme.Print() << std::endl;

  rme.RemoveRegion(&r2);
  std::cout << rme.Print() << std::endl;

  rme.AddRegion(&r2);
  std::cout << rme.Print() << std::endl;
}

