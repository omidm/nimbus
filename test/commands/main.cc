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
#include "shared/spawn_copy_job_command.h"

using namespace nimbus; // NOLINT

IDSet<logical_data_id_t> ldo_set_a;
IDSet<logical_data_id_t> ldo_set_b;
IDSet<job_id_t> job_set_a;
IDSet<job_id_t> job_set_b;
ID<logical_data_id_t> ldo_a;
ID<logical_data_id_t> ldo_b;
ID<job_id_t> job_id;
ID<job_id_t> parent;
ID<job_id_t> future;
bool sterile;
Parameter params;
std::string name = "test-string";

void initalize_sets() {
  for (int i = 0; i < 125; i++) {
    ldo_set_a.insert(10000 + i);
  }
  for (int i = 0; i < 119; i++) {
    ldo_set_b.insert(11000 + i);
  }
  for (int i = 0; i < 21; i++) {
    job_set_a.insert(20000 + i);
  }
  for (int i = 0; i < 34; i++) {
    job_set_b.insert(21000 + i);
  }

  ldo_a.set_elem(10125);
  ldo_b.set_elem(11119);
  job_id.set_elem(2014);
  parent.set_elem(1940);
  future.set_elem(1977);
  sterile = false;

  SerializedData d("serialized data");
  params.set_ser_data(d);
}


int test_compute_job() {
  std::cout << "Testing SpawnComputeJob." << std::endl;
  SpawnComputeJobCommand cmd((const std::string)name,
                             (const ID<job_id_t>)job_id,
                             (const IDSet<logical_data_id_t>)ldo_set_a,
                             (const IDSet<logical_data_id_t>)ldo_set_b,
                             (const IDSet<job_id_t>)job_set_a,
                             (const IDSet<job_id_t>)job_set_b,
                             (const ID<job_id_t>)parent,
                             (const ID<job_id_t>)future,
                             (const bool)sterile,
                             (const Parameter)params);

  std::string ToNetworkDataResult = cmd.ToNetworkData();
  SchedulerPBuf buf;
  buf.ParseFromString(ToNetworkDataResult);
  SpawnComputeJobCommand cmd2;
  bool result = cmd2.Parse(buf);
  
  if (result) {
    std::string tags = cmd2.ToString();
    std::cout << "Successful parse of " << result << std::endl;
    std::cout << "With tags: " << tags.size() << ", protobuf: " << ToNetworkDataResult.size() << std::endl;
    std::cout << "Tags: " << tags << std::endl;
    return 1;
  } else {
    std::cout << "Failed parse of buffer" << std::endl;
    std::cout << cmd2.ToString() << std::endl;
    const unsigned char* buf = (const unsigned char*)ToNetworkDataResult.c_str();
    for (unsigned int i = 0; i < ToNetworkDataResult.size(); i++) {
      printf("%02x ", buf[i]);
    }
    return 0;
  }
}

int test_copy_job() {
  std::cout << "Testing SpawnCopyJob." << std::endl;
  
  SpawnCopyJobCommand cmd((const ID<job_id_t>)job_id,
                          (const ID<logical_data_id_t>)ldo_a,
                          (const ID<logical_data_id_t>)ldo_b,
                          (const IDSet<job_id_t>)job_set_a,
                          (const IDSet<job_id_t>)job_set_b,
                          (const ID<job_id_t>)parent);

  std::string ToNetworkDataResult = cmd.ToNetworkData();
  SchedulerPBuf buf;
  buf.ParseFromString(ToNetworkDataResult);
  SpawnCopyJobCommand cmd2;
  bool result = cmd2.Parse(buf);
  
  if (result) {
    std::string tags = cmd2.ToString();
    std::cout << "Successful parse of " << result << std::endl;
    std::cout << "With tags: " << tags.size() << ", protobuf: " << ToNetworkDataResult.size() << std::endl;
    std::cout << "Tags: " << tags << std::endl;
    return 1;
  } else {
    std::cout << "Failed parse of buffer" << std::endl;
    std::cout << cmd2.ToString() << std::endl;
    const unsigned char* buf = (const unsigned char*)ToNetworkDataResult.c_str();
    for (unsigned int i = 0; i < ToNetworkDataResult.size(); i++) {
      printf("%02x ", buf[i]);
    }
    return 0;
  }
  return 1;
}
int test_define_data() {
  return 1;
}
int test_handshake() {
  return 1;
}
int test_jobdone() {
  return 1;
}
int test_execute_job() {
  return 1;
}
int test_define_partition() {
  return 1;
}
int test_terminate() {
  return 1;
}
int test_profile() {
  return 1;
}


int main(int argc, char** argv) {
  initalize_sets();
  test_compute_job();
  test_copy_job();
  test_define_data();
  test_handshake();
  test_jobdone();
  test_execute_job();
  test_define_partition();
  test_terminate();
  test_profile();
  return 0;
}
