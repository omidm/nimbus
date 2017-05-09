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
 * Jobs are defined here.
 *
 */

#ifndef NIMBUS_APPLICATIONS_SCAFFOLD_JOB_H_
#define NIMBUS_APPLICATIONS_SCAFFOLD_JOB_H_

#include <unistd.h>
#include <iostream> // NOLINT
#include "src/worker/physical_data_instance.h"
#include "src/shared/nimbus.h"
#include "./app.h"
#include "src/shared/dbg.h"

#define LOOP_JOB_NAME "loop"
#define STENCIL_JOB_NAME "stencil"

using namespace nimbus; // NOLINT

class Main : public Job {
  public:
    explicit Main(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Loop : public Job {
  public:
    explicit Loop(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Stencil : public Job {
  public:
    explicit Stencil(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};


#endif  // NIMBUS_APPLICATIONS_SCAFFOLD_JOB_H_
