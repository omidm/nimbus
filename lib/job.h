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
  * A Nimbus job. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_LIB_JOB_H_
#define NIMBUS_LIB_JOB_H_

#include <vector>
#include <string>
#include <set>
#include <map>
#include "lib/data.h"
#include "lib/idset.h"

namespace nimbus {

class Application;
class Job;
typedef std::map<int, Job*> JobMap;
typedef std::map<std::string, Job*> JobTable;
typedef std::vector<Data*> DataArray;

enum JobType {JOB_COMP, JOB_SYNC};

class Job {
 public:
  Job(Application* app, JobType type);
  virtual ~Job() {}

  virtual void execute(std::string params, const DataArray& da) {}
  virtual Job* clone();
  virtual void sleep() {}
  virtual void cancel() {}

  uint64_t id();
  void set_id(uint64_t id);

 protected:
  Application* application_;
  JobType type_;

 private:
  uint32_t id_;
  IDSet read_set_;
  IDSet write_set_;
  IDSet before_set_;
  IDSet after_set_;
  std::string parameters_;
};

}  // namespace nimbus
#endif  // NIMBUS_LIB_JOB_H_


