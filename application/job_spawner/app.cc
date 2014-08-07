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
 * An Example application that spawns a lot of jobs in the 3d space.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <vector>
#include "./app.h"
#include "./job.h"
#include "./data.h"

JobSpawnerApp::JobSpawnerApp(size_t counter,
                             size_t part_num,
                             size_t chunk_per_part,
                             size_t chunk_size,
                             size_t bandwidth,
                             size_t stage_num,
                             size_t job_length_usec) {
  counter_ = counter;
  part_num_ = part_num;
  chunk_per_part_ = chunk_per_part;
  chunk_size_ = chunk_size;
  bandwidth_ = bandwidth;
  stage_num_ = stage_num;
  job_length_usec_ = job_length_usec;
};

JobSpawnerApp::~JobSpawnerApp() {
};

void JobSpawnerApp::Load() {
  std::cout << "Start Creating Data and Job Tables" << std::endl;

  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(INIT_JOB_NAME, new Init());
  RegisterJob(LOOP_JOB_NAME, new ForLoop(this));
  RegisterJob(PRINT_JOB_NAME, new Print());
  RegisterJob(STAGE_JOB_NAME, new Stage(this));
  RegisterJob(CONNECTOR_JOB_NAME, new Connector(this));

  RegisterData(DATA_NAME, new Vec());

  std::cout << "Finished Creating Data and Job Tables" << std::endl;
};



