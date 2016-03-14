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

#ifndef NIMBUS_APPLICATIONS_ML_LOGISTIC_REGRESSION_APP_H_
#define NIMBUS_APPLICATIONS_ML_LOGISTIC_REGRESSION_APP_H_

#include <vector>
#include <iostream> // NOLINT
#include "src/worker/application.h"
#include "src/shared/nimbus_types.h"
#include "./data.h"

using nimbus::Application;

class LogisticRegression : public Application {
  public:
    LogisticRegression(const size_t& dimension_,
                       const size_t& iteration_num,
                       const size_t& partition_num,
                       const double& sample_num_m,
                       const size_t& reduction_partition_num);
    ~LogisticRegression();
    virtual void Load();

    virtual size_t dimension();
    virtual size_t iteration_num();
    virtual size_t partition_num();
    virtual double sample_num_m();
    virtual size_t reduction_partition_num();
    virtual size_t sample_num_per_partition();
    virtual bool automatic_reduction_active();
    virtual bool reduction_combiner_active();

    virtual void set_automatic_reduction_active(bool flag);
    virtual void set_reduction_combiner_active(bool flag);

  private:
    size_t dimension_;
    size_t iteration_num_;
    size_t partition_num_;
    double sample_num_m_;
    size_t reduction_partition_num_;
    size_t sample_num_per_partition_;
    bool automatic_reduction_active_;
    bool reduction_combiner_active_;
};

#endif  // NIMBUS_APPLICATIONS_ML_LOGISTIC_REGRESSION_APP_H_
