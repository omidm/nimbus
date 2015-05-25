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
  * A test for idset template class query. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <boost/unordered_set.hpp>
#include <iostream>  // NOLINT
#include <sstream>  // NOLINT
#include <string>
#include <map>
#include <set>
#include "shared/nimbus.h"
#include "application_utils/space_iterator.h"

using namespace nimbus; // NOLINT

int main(int argc, const char *argv[]) {
  if (argc < 4) {
    std::cout << "ERROR: provide three integers (SIZE, BW, PART_NUM)!" << std::endl;
    exit(-1);
  }

  int size;
  {
    std::string s(argv[1]);
    std::stringstream ss(s);
    ss >> size;
    if (ss.fail()) {
      std::cout << "ERROR: provide an integer for size!" << std::endl;
      exit(-1);
    }
  }

  int bw;
  {
    std::string s(argv[2]);
    std::stringstream ss(s);
    ss >> bw;
    if (ss.fail()) {
      std::cout << "ERROR: provide an integer for bw!" << std::endl;
      exit(-1);
    }
  }

  int part_num;
  {
    std::string s(argv[3]);
    std::stringstream ss(s);
    ss >> part_num;
    if (ss.fail()) {
      std::cout << "ERROR: provide an integer for part_num!" << std::endl;
      exit(-1);
    }
  }

  SpaceIterator::Cursor cursor;
  SpaceIterator iter(size, bw, part_num, true);
  iter.Initialize();
  do {
    iter.ReadCursor(&cursor);
    std::cout << " POINT: " << cursor.point_
              << " DELTA: " << cursor.delta_
              << " TYPE: " << cursor.type_
              << std::endl;
  } while (iter.Advance());
}

