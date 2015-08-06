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
 * Distributed Logistic Regression application
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_LOGISTIC_REGRESSION_APP_H_
#define NIMBUS_APPLICATION_LOGISTIC_REGRESSION_APP_H_

#include <vector>
#include <iostream> // NOLINT
#include "worker/application.h"
#include "shared/nimbus_types.h"
// #include "./job.h"
#include "./data.h"

using nimbus::Application;

class LogisticRegression : public Application {
  public:
    LogisticRegression(const size_t& dimension_,
                       const size_t& iteration_num,
                       const size_t& partition_num,
                       const size_t& data_size_mb);
    ~LogisticRegression();
    virtual void Load();

    virtual size_t dimension();
    virtual size_t iteration_num();
    virtual size_t partition_num();
    virtual size_t data_size_mb();
    virtual size_t sample_num_per_partition();

  private:
    size_t dimension_;
    size_t iteration_num_;
    size_t partition_num_;
    size_t data_size_mb_;
    size_t sample_num_per_partition_;
};

#endif  // NIMBUS_APPLICATION_LOGISTIC_REGRESSION_APP_H_
