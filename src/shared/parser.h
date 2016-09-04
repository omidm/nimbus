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
  * Parser for Nimbus scheduler protocol. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_PARSER_H_
#define NIMBUS_SRC_SHARED_PARSER_H_

#include <stdint.h>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include "src/shared/id.h"
#include "src/shared/idset.h"
#include "src/shared/serialized_data.h"
#include "src/shared/escaper.h"
#include "src/shared/parameter.h"
#include "src/shared/nimbus_types.h"

namespace nimbus {

// TODO(omidm): remove obsolete functions, defs.

typedef std::set<std::string> CmSet;
void parseCommand(const std::string& string,
                  const CmSet& commandSet,
                  std::string& command,
                  std::vector<int>& arguments);

int parseCommandFile(const std::string& fname,
                     CmSet& cs);

bool isSet(const std::string& tag);

int countOccurence(std::string str, std::string substr);

// ********************************************************


bool ParseWorkerDataHeader(const std::string& input,
                           job_id_t& receive_job_id,
                           job_id_t& mega_rcr_job_id,
                           size_t& data_length,
                           data_version_t& version,
                           template_id_t& template_generation_id);


}  // namespace nimbus
#endif  // NIMBUS_SRC_SHARED_PARSER_H_
