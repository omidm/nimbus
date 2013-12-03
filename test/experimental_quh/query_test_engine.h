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
  * Helper class only for testing the idea of job query.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef QUERY_TEST_ENGINE_H_
#define QUERY_TEST_ENGINE_H_

#include <list>
#include <string>
#include "shared/nimbus_types.h"
#include "test/experimental_quh/job_index.h"
#include "test/experimental_quh/variable_region_set.h"

namespace nimbus {

class QueryResult {
 public:
  QueryResult() {}
  virtual ~QueryResult() {}
  VariableRegionSet& get_read_set_buffer() { return _read_set_buffer; }
  const VariableRegionSet& get_read_set_buffer() const { return _read_set_buffer; }
  void set_read_set_buffer(const VariableRegionSet& read_set_buffer) {
    _read_set_buffer.CopyFrom(read_set_buffer);
  }
  VariableRegionSet& get_write_set_buffer() { return _write_set_buffer; }
  const VariableRegionSet& get_write_set_buffer() const { return _write_set_buffer; }
  void set_write_set_buffer(const VariableRegionSet& write_set_buffer) {
    _write_set_buffer.CopyFrom(write_set_buffer);
  }
  JobIdSet& get_before_set() { return _before_set; }
  const JobIdSet& get_before_set() const { return _before_set; }
  void set_before_set(const JobIdSet& id_set) {
    _before_set = id_set;
  }
 private:
  VariableRegionSet _read_set_buffer;
  VariableRegionSet _write_set_buffer;
  JobIdSet _before_set;
};

class QueryTestEngine {
 public:
  QueryTestEngine();
  virtual ~QueryTestEngine();
  bool Query(
      const VariableRegionSet& read_set,
      const VariableRegionSet& write_set,
      QueryResult* query_result);
  bool Add(const std::string& comment, const int label,
           const QueryResult& query_result);
  void PruneResult();
  bool WriteOut(const std::string& filename);
 private:
  static const job_id_t START = 1000;
  job_id_t _fake_id_counter;
  JobIndex _job_index;
  struct ResultEntry {
    job_id_t job_id;
    JobIdSet before_set;
    std::string comment;
    int label;
  };
  std::list<ResultEntry*> _result_list;
};

}  // namespace nimbus

#endif  // QUERY_TEST_ENGINE_H_
