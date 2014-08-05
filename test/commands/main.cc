/*
 * Copyright 2014 Stanford University.
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
 * Program to test command serialization/deserialization.
 *
 * Author: Philip Levis <pal@cs.stanford.edu>
 */

#include <string>
#include <iostream>
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/spawn_compute_job_command.h"

using namespace nimbus; // NOLINT

int test_compute_job() {
  std::string name = "compute-test";
  ID<job_id_t> id(53);
  IDSet<logical_data_id_t> read;
  for (int i = 0; i < 125; i++) {
    read.insert(2000 + i);
  }
  IDSet<logical_data_id_t> write;
  for (int i = 0; i < 125; i++) {
    write.insert(1000 + i);
  }
  IDSet<job_id_t> before;
  before.insert(43);
  IDSet<job_id_t> after;

  ID<job_id_t> parent(23);
  ID<job_id_t> future(0);

  bool sterile = false;
  const Parameter params;

  SpawnComputeJobCommand cmd((const std::string)name,
                             (const ID<job_id_t>)id,
                             (const IDSet<logical_data_id_t>)read,
                             (const IDSet<logical_data_id_t>)write,
                             (const IDSet<job_id_t>)before,
                             (const IDSet<job_id_t>)after,
                             (const ID<job_id_t>)parent,
                             (const ID<job_id_t>)future,
                             (const bool)sterile,
                             params);

  std::string toStringResult = cmd.toString();
  SpawnComputeJobCommand cmd2;
  bool result = cmd2.Parse(toStringResult);
  
  if (result) {
    std::string tags = cmd2.toStringWTags();
    std::cout << "Successful parse of " << result << std::endl;
    std::cout << "With tags: " << tags.size() << ", protobuf: " << toStringResult.size() << std::endl;
    std::cout << "Tags: " << tags << std::endl;
    return 1;
  } else {
    std::cout << "Failed parse of buffer" << std::endl;
    std::cout << cmd2.toStringWTags() << std::endl;
    const unsigned char* buf = (const unsigned char*)toStringResult.c_str();
    for (unsigned int i = 0; i < toStringResult.size(); i++) {
      printf("%02x ", buf[i]);
    }
    return 0;
  }
}

int main(int argc, char** argv) {
  test_compute_job();
  
  return 0;
}
