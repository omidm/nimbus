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
  * A Nimbus worker for job spawner application. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <iostream>  // NOLINT
#include <sstream> // NOLINT
#include "shared/nimbus.h"
#include "shared/distributed_db.h"

using namespace nimbus; // NOLINT

void PrintUsage() {
  std::cout << "ERROR: wrong arguments\n";
  std::cout << "Usage:\n";
  std::cout << "./worker\n";
  std::cout << "REQUIRED ARGUMENTS:\n";
  std::cout << "\t-ip [ip_address] -id [worker_id]\n";
  std::cout << "OPTIONIAL:\n";
  std::cout << "\t-? [?]\n";
}

int main(int argc, char *argv[]) {
  worker_id_t worker_id;
  std::string ip_address;
  bool ip_address_given = false;
  bool worker_id_given = false;

  std::cout << argc << std::endl;


  if (((argc - 1) % 2 != 0) || (argc != 5)) {
    PrintUsage();
    exit(-1);
  }

  for (int i = 1; i < argc; i = i + 2) {
    std::string tag = argv[i];
    std::string val = argv[i+1];
    if (tag == "-ip") {
      ip_address = val;
      ip_address_given = true;
    } else if (tag == "-id") {
      std::stringstream ss(val);
      ss >> worker_id;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
      worker_id_given = true;
    } else {
      PrintUsage();
      exit(-1);
    }
  }

  if (!ip_address_given || !worker_id_given) {
    PrintUsage();
    exit(-1);
  }

  DistributedDB ddb;
  ddb.Initialize("localhost", worker_id);

  std::cout << "Inserting key-value pair ...\n";

  std::string h1, h2;
  ddb.Put("key1", "value1", 3, &h1);
  ddb.Put("key2", "value2", 4, &h2);

  std::cout << "Looking up from same node ...\n";

  std::string v1, v2;
  ddb.Get(h1, &v1);
  std::cout << v1 << std::endl;
  ddb.Get(h2, &v2);
  std::cout << v2 << std::endl;

  std::cout << "Looking up from another node ...\n";

  DistributedDB ddb2;
  ddb2.Initialize("127.0.0.1", worker_id + 1);

  ddb2.Get(h1, &v1);
  std::cout << v1 << std::endl;
  ddb2.Get(h2, &v2);
  std::cout << v2 << std::endl;

  std::cout << "END\n";
}

