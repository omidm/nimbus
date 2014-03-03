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
 * Job classes for job spawner application.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_JOB_SPAWNER_JOB_H_
#define NIMBUS_APPLICATION_JOB_SPAWNER_JOB_H_

#include <unistd.h>
#include <iostream> // NOLINT
#include "worker/physical_data_instance.h"
#include "shared/nimbus.h"

#define INIT_JOB_NAME "init"
#define PRINT_JOB_NAME "print"
#define LOOP_JOB_NAME "for_loop"
#define STAGE_JOB_NAME "stage"
#define CONNECTOR_JOB_NAME "connector"

using namespace nimbus; // NOLINT

class Main : public Job {
  public:
    explicit Main(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Init : public Job {
  public:
    Init();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Print : public Job {
  public:
    Print();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class ForLoop : public Job {
  public:
    explicit ForLoop(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Stage : public Job {
  public:
    explicit Stage(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Connector : public Job {
  public:
    explicit Connector(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};


#endif  // NIMBUS_APPLICATION_JOB_SPAWNER_JOB_H_
