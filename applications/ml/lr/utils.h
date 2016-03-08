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

#ifndef NIMBUS_APPLICATIONS_ML_LR_UTILS_H_
#define NIMBUS_APPLICATIONS_ML_LR_UTILS_H_

#include <iostream> // NOLINT
#include <vector>
#include <algorithm>
#include "src/worker/application.h"
#include "src/shared/nimbus_types.h"
#include "src/shared/parameter.h"
#include "./data.h"

using namespace nimbus; // NOLINT

bool LoadParameter(Parameter *parameter, size_t *value);

bool SerializeParameter(Parameter *parameter, size_t value);

double VectorDotProduct(const std::vector<double>* vec1,
                        const std::vector<double>* vec2);

void VrctorScale(std::vector<double>* vec,
                 const double& scale);

void VectorAddWithScale(std::vector<double>* acc,
                        const std::vector<double>* add,
                        const double& scale);

void PrintWeight(Weight* w,
                 size_t loop_counter,
                 size_t max_loop);

#endif  // NIMBUS_APPLICATIONS_ML_LR_UTILS_H_
