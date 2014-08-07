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

#ifndef NIMBUS_APPLICATION_JOB_SPAWNER_APP_H_
#define NIMBUS_APPLICATION_JOB_SPAWNER_APP_H_

#include <iostream> // NOLINT
#include "worker/application.h"
#include "shared/nimbus_types.h"

using nimbus::Application;

class JobSpawnerApp : public Application {
  public:
    JobSpawnerApp(size_t counter,
                  size_t part_num,
                  size_t chunk_per_part,
                  size_t chunk_size,
                  size_t bandwidth,
                  size_t stage_num,
                  size_t job_length_usec_);
    ~JobSpawnerApp();
    virtual void Load();

    size_t counter_;
    size_t part_num_;
    size_t chunk_per_part_;
    size_t chunk_size_;
    size_t bandwidth_;
    size_t stage_num_;
    size_t job_length_usec_;
};



#endif  // NIMBUS_APPLICATION_JOB_SPAWNER_APP_H_
