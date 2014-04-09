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
  * A test for hash function. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <boost/functional/hash.hpp>
#include <iostream>  // NOLINT
#include <sstream>  // NOLINT
#include <string>
#include <map>
#include "shared/nimbus.h"
#include "shared/log.h"

#define CAL_NUM 1000000

using namespace nimbus; // NOLINT
using boost::hash;

int main(int argc, const char *argv[]) {
  if (argc < 2) {
    std::cout << "ERROR: provide an input!" << std::endl;
    exit(-1);
  }

  std::string input(argv[1]);

  nimbus_initialize();
  Log log;
  hash<std::string> hash_function;
  std:: cout << "Hash value: " << hash_function(input) << std::endl;

  size_t result = 0;
  log.StartTimer();
  for (int i = 0; i < CAL_NUM; ++i) {
    result += hash_function(input);
  }
  log.StopTimer();
  std::cout << "result: " << result << std::endl;
  std::cout << "Time elapsed for " << CAL_NUM << " hash calculations: " << log.timer() << std::endl;
}

