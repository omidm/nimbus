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
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <boost/functional/hash.hpp>
#include <ctime>
#include <sstream>
#include <string>
#include "worker/worker.h"
#include "data/physbam/physbam_data.h"

#include "worker/util_dumping.h"

using boost::hash;

namespace nimbus {

void DumpVersionInformation(
    Job *job, const DataArray& da, Log *log, std::string tag) {
#ifndef MUTE_LOG
  std::string input = "";
  for (size_t i = 0; i < da.size(); ++i) {
    std::ostringstream ss_l;
    ss_l << da[i]->logical_id();
    input += ss_l.str();
    input += " : ";
    // std::ostringstream ss_p;
    // ss_p << da[i]->physical_id();
    // input += ss_p.str();
    // input += " : ";
    std::ostringstream ss_v;
    ss_v << da[i]->version();
    input += ss_v.str();
    input += " - ";
  }
  hash<std::string> hash_function;

  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "%s name: %s id: %llu  version_hash: %lu versions: %s",
           tag.c_str(), job->name().c_str(), job->id().elem(),
           hash_function(input), input.c_str());
  log->WriteToFile(std::string(buff), LOG_INFO);
#endif
}

void DumpDataHashInformation(Job *job, const DataArray& da, Log *log, std::string tag) {
#ifndef MUTE_LOG
  if ((dynamic_cast<CreateDataJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<LocalCopyJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<RemoteCopySendJob*>(job) == NULL) && // NOLINT
      (dynamic_cast<RemoteCopyReceiveJob*>(job) == NULL)) { // NOLINT
    std::string input = "";
    for (size_t i = 0; i < da.size(); ++i) {
      if (dynamic_cast<PhysBAMData*>(da[i]) != NULL) { // NOLINT
        std::ostringstream ss_l;
        ss_l << da[i]->logical_id();
        input += ss_l.str();
        input += " : ";
        // std::ostringstream ss_p;
        // ss_p << da[i]->physical_id();
        // input += ss_p.str();
        // input += " : ";
        std::ostringstream ss_v;
        ss_v << dynamic_cast<PhysBAMData*>(da[i])->HashCode(); // NOLINT
        input += ss_v.str();
        input += " - ";
      }
    }
    hash<std::string> hash_function;

    char buff[LOG_MAX_BUFF_SIZE];
    snprintf(buff, sizeof(buff),
        "%s name: %s id: %llu  aggregate_hash: %lu hashes: %s",
        tag.c_str(), job->name().c_str(), job->id().elem(),
        hash_function(input), input.c_str());
    log->WriteToFile(std::string(buff), LOG_INFO);
  }
#endif
}


void DumpDataOrderInformation(Job *job, const DataArray& da, Log *log,
                              std::string tag) {
#ifndef MUTE_LOG
  std::string input = "";
  for (size_t i = 0; i < da.size(); ++i) {
    std::ostringstream ss_l;
    ss_l << da[i]->logical_id();
    input += ss_l.str();
    input += " - ";
  }
  hash<std::string> hash_function;

  char buff[LOG_MAX_BUFF_SIZE];
  snprintf(buff, sizeof(buff),
      "%s name: %s id: %llu  order_hash: %lu logical ids: %s",
           tag.c_str(), job->name().c_str(), job->id().elem(),
           hash_function(input), input.c_str());
  log->WriteToFile(std::string(buff), LOG_INFO);
#endif
}

}  // namespace nimbus
