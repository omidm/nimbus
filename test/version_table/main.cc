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
  * A test for version table operations. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <iostream>  // NOLINT
#include <sstream>  // NOLINT
#include <string>
#include <map>
#include <vector>
#include "shared/nimbus.h"
#include "shared/log.h"
#include "scheduler/version_table.h"
#include "scheduler/version_operator.h"

#define VERSION_RANGE 20
#define SEED 123

using namespace nimbus; // NOLINT

int main(int argc, const char *argv[]) {
  nimbus_initialize();
  Log log;
  unsigned int seed = SEED;

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


  VersionOperator op;
  boost::shared_ptr<VersionTable> t_1(new VersionTable(op.GetNewVersionTableId()));
  boost::shared_ptr<VersionTable> t_2(new VersionTable(op.GetNewVersionTableId()));
  std::vector<boost::shared_ptr<VersionTable> > tables;
  std::vector<boost::shared_ptr<const VersionTable> > const_tables;

  log.StartTimer();
  for (int i = 0; i < elem_num; ++i) {
    t_1->set_entry(i, rand_r(&seed) % (VERSION_RANGE));
  }
  log.StopTimer();
  std::cout << "Time elapsed insertion: " << log.timer() << std::endl;

  log.StartTimer();
  VersionTable::Map content;
  for (int i = 0; i < elem_num; ++i) {
    content[i] = rand_r(&seed) % (VERSION_RANGE);
  }
  t_1->set_content(content);
  log.StopTimer();
  std::cout << "Time elapsed insertion: " << log.timer() << std::endl;

  std::cout << "Print before recomputing root:\n";
  // t_1->Print();

  std::cout << "Print after recomputing root:\n";
  tables.clear();
  tables.push_back(t_1);
  op.RecomputeRootForVersionTables(tables);
  // t_1->Print();


  data_version_t version;
  log.StartTimer();
  if (t_1->query_entry(elem_num * 2, &version)) {
    std::cout << "query result: " << version << std::endl;
    log.StopTimer();
    std::cout << "Time elapsed query: " << log.timer() << std::endl;
  }

  for (int i = 0; i < elem_num; ++i) {
    t_2->set_entry(i, rand_r(&seed) % (VERSION_RANGE));
  }

  tables.clear();
  tables.push_back(t_2);
  op.RecomputeRootForVersionTables(tables);
  // t_2->Print();

  const_tables.clear();
  const_tables.push_back(t_1);
  const_tables.push_back(t_2);
  boost::shared_ptr<VersionTable> merge;
  log.StartTimer();
  op.MergeVersionTables(const_tables, &merge);
  log.StopTimer();
  std::cout << "Time elapsed merging: " << log.timer() << std::endl;
  // merge->Print();


  op.FlushCache();

  std::cout << "Print after recomputing root for both tables:\n";
  tables.clear();
  tables.push_back(t_1);
  tables.push_back(t_2);
  op.RecomputeRootForVersionTables(tables);
  // t_1->Print();
  // t_2->Print();
  std::cout << "Root size: " << t_1->root()->size() << std::endl;
  std::cout << "content table 1 size: " << t_1->content_p()->size() << std::endl;
  std::cout << "content table 2 size: " << t_2->content_p()->size() << std::endl;

  const_tables.clear();
  const_tables.push_back(t_1);
  const_tables.push_back(t_2);
  log.StartTimer();
  op.MergeVersionTables(const_tables, &merge);
  log.StopTimer();
  std::cout << "Time elapsed merging: " << log.timer() << std::endl;
  // merge->Print();
}

