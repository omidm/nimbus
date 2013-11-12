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
 * An Example application that is meant to run over multiple workers.
 * It is simply applying a stencil over a one dimensional array.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <string>
#include <iostream>
#include "protobufs/vector_msg.pb.h"

using vector_msg::VectorMsg;

bool Serialize() {
  VectorMsg vec_msg;
  int size_ = 1;

  for (int i = 0; i < size_; i++) {
    std::cout << "****" << i << std::endl;
    vec_msg.add_elem(i);
  }

  // Preliminary string test
  printf("Testing static allocation.\n");
  std::string staticStr;
  bool result = vec_msg.SerializeToString(&staticStr);
  
  printf("Testing dynamic allocation.\n");
  char* testArray = new char[128];
  //  strncpy(testArray, "test text", 10);
  testArray[0] = 0;
  std::string* prelim = new std::string(testArray);
  prelim->clear();
  prelim->resize(10);
  


  printf("Vector message string test.\n");
  //std::string str = "test";
  //  str = str + " test ";
  std::string* str = new std::string(testArray);
  result = vec_msg.SerializeToString(str);

  testArray[0] = 0;
  str = new std::string(testArray);
  result = vec_msg.SerializeToString(str);
  if (result) {
    printf("Serialzed to string correctly.\n");
    const char* ptr = str->c_str();
    printf("\t");
    for (uint32_t i = 0; i < str->length(); i++) {
      printf("%02hx ", ptr[i]);
    }
    printf("\n");
  } else {
    printf("Serialized to string incorrectly.\n");
  }
  return true;
}

int main(int argc, char** argv) {
  
  Serialize();
  
}
