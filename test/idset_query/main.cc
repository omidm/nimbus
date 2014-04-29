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
#include "shared/log.h"
#include "../../application/job_spawner/app.h"

#define START_ 100000

using namespace nimbus; // NOLINT

int main(int argc, const char *argv[]) {
  if (argc < 2) {
    std::cout << "ERROR: provide an integer!" << std::endl;
    exit(-1);
  }

  int elem_num;
  std::string s(argv[1]);
  std::stringstream ss(s);
  ss >> elem_num;
  if (ss.fail()) {
    std::cout << "ERROR: provide an integer!" << std::endl;
    exit(-1);
  }

  nimbus_initialize();
  Log log;

  log.StartTimer();
  // IDSet<job_id_t> set;
  // std::set<job_id_t> set;
  boost::unordered_set<job_id_t> set;
  for (int i = 0; i < elem_num; ++i) {
    set.insert(START_ + i);
  }
  log.StopTimer();
  std::cout << "Time elapsed idset insertion: " << log.timer() << std::endl;

  // IDSet<job_id_t> set_copy;
  // std::set<job_id_t> set_copy;
  boost::unordered_set<job_id_t> set_copy;
  log.StartTimer();
  set_copy = set;
  log.StopTimer();
  std::cout << "Time elapsed idset copy: " << log.timer() << std::endl;

//   log.StartTimer();
//   if (set.contains(START_ + 1)) {
//     log.StopTimer();
//     std::cout << "Time elapsed idset query: " << log.timer() << std::endl;
//   }

  log.StartTimer();
  std::map<int, int> map;
  for (int i = 0; i < elem_num; ++i) {
    map[START_ + i] = START_ + i;
  }
  log.StopTimer();
  std::cout << "Time elapsed map insertion: " << log.timer() << std::endl;

  std::map<int, int> map_copy;
  log.StartTimer();
  map_copy = map;
  log.StopTimer();
  std::cout << "Time elapsed map copy: " << log.timer() << std::endl;

  log.StartTimer();
  if (map.count(START_ + 1)) {
    log.StopTimer();
    std::cout << "Time elapsed map query: " << log.timer() << std::endl;
  }
}

