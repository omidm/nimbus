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

#include "src/shared/parser.h"

using boost::tokenizer;
using boost::char_separator;

namespace nimbus {

// TODO(omidm): remove these obsolete parse functions.

void parseCommand(const std::string& str, const CmSet& cms,
                  std::string& cm, std::vector<int>& args) {
  cm.clear();
  args.clear();
  std::set<std::string>::iterator it;
  for (it = cms.begin(); it != cms.end(); ++it) {
    if (str.find(*it) == 0) {
      cm = *it;
      int arg;
      std::stringstream ss;
      ss << str.substr(it->length(), std::string::npos);
      while (true) {
        ss >> arg;
        if (ss.fail()) {
          break;
        }
        args.push_back(arg);
      }
      break;
    }
  }
  if (cm == "") {
    std::cout << "wrong command! try again." << std::endl;
  }
}


int parseCommandFile(const std::string& fname, CmSet& cs) {
  return 0;
}

int countOccurence(std::string str, std::string substr) {
  int count = 0;
  std::size_t pos = -1;

  while (true) {
    pos = str.find(substr, pos + 1);
    if (pos != std::string::npos)
      count++;
    else
      break;
  }
  return count;
}


bool isSet(const std::string& str) {
  if (str.find("{") == std::string::npos)
    return false;
  return true;
}

// ********************************************************

bool ParseWorkerDataHeader(const std::string& input,
                           job_id_t& receive_job_id,
                           job_id_t& mega_rcr_job_id,
                           size_t& data_length,
                           data_version_t& version,
                           template_id_t& template_generation_id) {
  int num = 5;
  job_id_t temp_j;
  job_id_t temp_m;
  size_t temp_d;
  data_version_t temp_v;
  template_id_t temp_t;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(input, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: Data header has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: Data Exchanger Header has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  std::stringstream ss_j(*iter);
  ss_j >> temp_j;
  if (ss_j.fail()) {
    std::cout << "ERROR: wrong element for receive job id." << std::endl;
    return false;
  } else {
    receive_job_id = temp_j;
  }

  iter++;
  std::stringstream ss_m(*iter);
  ss_m >> temp_m;
  if (ss_m.fail()) {
    std::cout << "ERROR: wrong element for mega rcr job id." << std::endl;
    return false;
  } else {
    mega_rcr_job_id = temp_m;
  }

  iter++;
  std::stringstream ss_d(*iter);
  ss_d >> temp_d;
  if (ss_d.fail()) {
    std::cout << "ERROR: wrong element for data length." << std::endl;
    return false;
  } else {
    data_length = temp_d;
  }

  iter++;
  std::stringstream ss_v(*iter);
  ss_v >> temp_v;
  if (ss_v.fail()) {
    std::cout << "ERROR: wrong element for data length." << std::endl;
    return false;
  } else {
    version = temp_v;
  }

  iter++;
  std::stringstream ss_t(*iter);
  ss_t >> temp_t;
  if (ss_t.fail()) {
    std::cout << "ERROR: wrong element for template generation id." << std::endl;
    return false;
  } else {
    template_generation_id = temp_t;
  }

  return true;
}

}  // namespace nimbus

