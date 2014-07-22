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
#include "scheduler/job_entry.h"

#define LDO_NUM 90000
#define JOB_NUM 946
#define LDO_NUM_PER_JOB 463
#define SEED 12

using namespace nimbus; // NOLINT

unsigned int seed = SEED;

void LoadRandomSet(IDSet<logical_data_id_t>* set) {
  // int max_num = LDO_NUM_PER_JOB * JOB_NUM;
  int max_num = LDO_NUM;
  if (max_num > RAND_MAX)
    max_num = RAND_MAX;
  while (set->size() < LDO_NUM_PER_JOB)
    set->insert(rand_r(&seed) % max_num + 1);
}


int main(int argc, const char *argv[]) {
  nimbus_initialize();
  Log log;


  typedef boost::unordered_set<JobEntry*> Pool;
  typedef boost::unordered_map<logical_data_id_t, Pool> Table;
  // typedef std::set<JobEntry*> Pool;
  // typedef std::map<logical_data_id_t, Pool> Table;

  std::vector<JobEntry*> jobs;
  Table table;

  log.StartTimer();

  for (int i = 0; i < JOB_NUM; ++i) {
    IDSet<logical_data_id_t> read_set;
    LoadRandomSet(&read_set);
    JobEntry* job = new JobEntry();
    job->set_read_set(read_set);
    jobs.push_back(job);
    // std::cout << job->read_set().toString() << std::endl;
  }
  log.StopTimer();
  std::cout << "Time elapsed making the read sets: " << log.timer() << std::endl;

  log.StartTimer();
  for (int i = 0; i < JOB_NUM; ++i) {
    JobEntry* job = jobs[i];
    IDSet<logical_data_id_t>::IDSetIter iter = job->read_set_p()->begin();
    for (; iter != job->read_set_p()->end(); ++iter) {
      table[*iter].insert(job);
    }
  }
  log.StopTimer();
  std::cout << "Time elapsed building the table: " << log.timer() << std::endl;
  std::cout << "Table size: " << table.size() << std::endl;
  std::cout << "RAND_MAX: " << RAND_MAX << std::endl;




/*
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


*/
}

