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
 * Helper functions in the application.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <math.h>
#include "./utils.h"
#include "protobuf_compiled/parameter_msg.pb.h"

using parameter_msg::ParameterMsg;

bool LoadParameter(Parameter *parameter, size_t *value) {
  std::string params_str(parameter->ser_data().data_ptr_raw(),
                         parameter->ser_data().size());
  ParameterMsg data_msgs;
  data_msgs.ParseFromString(params_str);
  assert(data_msgs.elem_size() == 1);
  *value = data_msgs.elem(0);
  return true;
}

bool SerializeParameter(Parameter *parameter, size_t value) {
  std::string params_str;
  ParameterMsg data_msgs;
  data_msgs.add_elem(value);
  data_msgs.SerializeToString(&params_str);
  parameter->set_ser_data(SerializedData(params_str));
  return true;
}

double VectorDistance(const std::vector<double>* vec1,
                      const std::vector<double>* vec2) {
  assert(vec1->size() == vec2->size());
  double result = 0;
  for (size_t i = 0; i < vec1->size(); ++i) {
    result += pow(vec1->operator[](i) - vec2->operator[](i), 2);
  }
  return pow(result, 0.5);
}

void VectorScale(std::vector<double>* vec,
                 const double& scale) {
  for (size_t i = 0; i < vec->size(); ++i) {
    vec->operator[](i) *= vec->operator[](i) * scale;
  }
}

void VectorAddWithScale(std::vector<double>* acc,
                        const std::vector<double>* add,
                        const double& scale) {
  assert(acc->size() == add->size());
  for (size_t i = 0; i < acc->size(); ++i) {
    acc->operator[](i) += add->operator[](i) * scale;
  }
}

