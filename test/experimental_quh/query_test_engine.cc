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

#include <list>
#include <fstream>  // NOLINT
#include <string>
#include "shared/nimbus_types.h"
#include "test/experimental_quh/job_index.h"
#include "test/experimental_quh/variable_region_set.h"
#include "test/experimental_quh/query_test_engine.h"

namespace nimbus {

QueryTestEngine::QueryTestEngine() {
  _fake_id_counter = START;
  _result_list.clear();
}

QueryTestEngine::~QueryTestEngine() {
  while (!_result_list.empty()) {
    delete _result_list.back();
    _result_list.pop_back();
  }
}

bool QueryTestEngine::Query(
    const VariableRegionSet& read_set,
    const VariableRegionSet& write_set,
    QueryResult* query_result) {
  assert(query_result != NULL);
  query_result->set_read_set_buffer(read_set);
  query_result->set_write_set_buffer(write_set);
  _job_index.QueryJobEntry(read_set, write_set,
      &query_result->get_before_set());
  return true;
}

bool QueryTestEngine::Add(const std::string& comment, const int label,
    const QueryResult& query_result) {
  ResultEntry* temp = new ResultEntry;
  temp->comment = comment;
  temp->label = label;
  temp->job_id = (_fake_id_counter++);
  temp->before_set = query_result.get_before_set();
  _job_index.AddJobEntry(
      temp->job_id,
      query_result.get_read_set_buffer(),
      query_result.get_write_set_buffer());
  _result_list.push_back(temp);
  return true;
}

bool QueryTestEngine::WriteOut(const std::string& filename) {
  std::ofstream fout(filename.c_str(), std::ofstream::out);
  fout << "digraph Workflow {" << std::endl;
  for (std::list<ResultEntry*>::iterator index = _result_list.begin();
       index != _result_list.end();
       index++) {
    fout << "\t" << "node" << (*index)->job_id
         << " [label=\"" << (*index)->comment << " " << (*index)->label
         << "\"];" << std::endl;
    for (JobIdSet::IDSetIter prec = (*index)->before_set.begin();
         prec != (*index)->before_set.end();
         prec++) {
      fout << "\tnode" << (*prec) << "->node" << (*index)->job_id << std::endl;
    }
  }
  fout << "}" << std::endl;
  fout.close();
  return true;
}

void QueryTestEngine::PruneResult() {
  typedef std::list<ResultEntry*>::iterator IterType;
  IterType index = _result_list.end();
  while (index != _result_list.begin()) {
    index--;
    for (IterType prec = _result_list.begin();
         prec != index;
         prec++)
      if ((*index)->before_set.contains((*prec)->job_id))
        for (JobIdSet::IDSetIter inner_index = (*prec)->before_set.begin();
             inner_index != (*prec)->before_set.end();
             inner_index++)
          if ((*inner_index) != (*prec)->job_id) {
            (*index)->before_set.remove(*inner_index);
          }
  }  // end while
}

}  // namespace nimbus
