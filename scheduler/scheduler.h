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
  * Nimbus scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_SCHEDULER_H_
#define NIMBUS_SCHEDULER_SCHEDULER_H_

#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include "shared/nimbus_types.h"
#include "shared/log.h"
#include "shared/cluster.h"
#include "shared/parser.h"
#include "shared/scheduler_server.h"
#include "scheduler/data_manager.h"
#include "scheduler/job_manager.h"
#include "scheduler/load_balancer.h"
#include "scheduler/job_assigner.h"
#include "scheduler/version_manager.h"
#include "shared/id_maker.h"

namespace nimbus {

class Scheduler {
  public:
    explicit Scheduler(port_t listening_port);
    virtual ~Scheduler();

    virtual void Run();

    virtual void set_max_job_to_assign(size_t num);
    virtual void set_min_worker_to_join(size_t num);
    virtual void set_job_assigner_thread_num(size_t num);
    virtual void set_max_command_process_num(size_t num);

    /* TODO(omidm): figure out what we want for these methods. */
    virtual void LoadClusterMap(std::string) {}
    virtual void DeleteWorker(SchedulerWorker * worker) {}
    virtual void AddWorker(SchedulerWorker * worker) {}
    virtual SchedulerWorker* GetWorker(int workerId) {return NULL;}

  protected:
    virtual void SchedulerCoreProcessor();

    virtual void RegisterInitialWorkers();

    virtual size_t RegisterPendingWorkers();

    virtual size_t AssignReadyJobs();

    virtual size_t RemoveObsoleteJobEntries();

    virtual void AddMainJob();

    virtual void TerminationProcedure();

    virtual void ProcessQueuedSchedulerCommands();

    virtual void ProcessSchedulerCommand(SchedulerCommand* cm);
    virtual void ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm);
    virtual void ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm);
    virtual void ProcessDefineDataCommand(DefineDataCommand* cm);
    virtual void ProcessHandshakeCommand(HandshakeCommand* cm);
    virtual void ProcessJobDoneCommand(JobDoneCommand* cm);
    virtual void ProcessDefinePartitionCommand(DefinePartitionCommand* cm);
    virtual void ProcessTerminateCommand(TerminateCommand* cm);

    virtual void CreateIDMaker();
    virtual void CreateSchedulerServer();
    virtual void CreateDataManager();
    virtual void CreateJobManager();
    virtual void CreateLoadBalancer();
    virtual void CreateJobAssigner();

    virtual void SetupIDMaker();
    virtual void SetupSchedulerServer();
    virtual void SetupDataManager();
    virtual void SetupJobManager();
    virtual void SetupLoadBalancer();
    virtual void SetupJobAssigner();
    virtual void SetupUserInterface();

    virtual void GetUserCommand();

    virtual void LoadUserCommands();
    virtual void LoadWorkerCommands();

    IDMaker *id_maker_;
    SchedulerServer* server_;
    JobManager* job_manager_;
    DataManager* data_manager_;
    JobAssigner *job_assigner_;
    LoadBalancer *load_balancer_;

    boost::thread* user_interface_thread_;
    boost::thread* scheduler_server_thread_;
    boost::thread* load_balancer_thread_;

    port_t listening_port_;
    size_t registered_worker_num_;
    bool terminate_application_flag_;
    exit_status_t terminate_application_status_;

    size_t max_job_to_assign_;
    size_t min_worker_to_join_;
    size_t job_assigner_thread_num_;
    size_t max_command_process_num_;

    CmSet user_command_set_;
    SchedulerCommand::PrototypeTable worker_command_table_;

    Log log_;
    Log log_receive_stamp_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_SCHEDULER_H_
