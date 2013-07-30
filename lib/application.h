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
  * Nimbus abstraction of an application. Programmers use this base class to
  * write various application served by Nimbus. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#ifndef NIMBUS_LIB_APPLICATION_H_
#define NIMBUS_LIB_APPLICATION_H_

#include <map>
#include <string>
#include <vector>
#include "lib/job.h"
#include "lib/data.h"
#include "lib/scheduler_client.h"

class Application;
typedef std::map<int, Application*> AppMap;

class Application {
  public:
    int id;
    int priority;

    JobTable jobTable;
    DataTable dataTable;
    SchedulerClient* client;

    Application();
    virtual void load();
    virtual void start(SchedulerClient* scheduler);

    void registerJob(std::string name, Job* j);

    void registerData(std::string name, Data* d);

    void spawnJob(std::string name, int id, IDSet bfore, IDSet after,
        IDSet read, IDSet write, std::string params);

    void defineData(std::string name, int id);

    void getNewJobID(int req_num, std::vector<int>* result);

    void getNewDataID(int req_num, std::vector<int>* result);
};


#endif  // NIMBUS_LIB_APPLICATION_H_
