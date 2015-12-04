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
#include "scheduler/template_manager.h"
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

    virtual void set_cleaner_thread_active(bool flag);
    virtual void set_bouncer_thread_active(bool flag);

    virtual void set_controller_template_active(bool flag);
    virtual void set_complex_memoization_active(bool flag);
    virtual void set_binding_memoization_active(bool flag);
    virtual void set_cascaded_binding_active(bool flag);
    virtual void set_worker_template_active(bool flag);
    virtual void set_mega_rcr_job_active(bool flag);
    virtual void set_data_manager_query_cache_active(bool flag);

    virtual void set_load_balancing_active(bool flag);
    virtual void set_load_balancing_period(int64_t period);
    virtual void set_fault_tolerance_active(bool flag);
    virtual void set_checkpoint_creation_period(int64_t period);

    virtual void set_assign_batch_size(size_t num);
    virtual void set_remove_batch_size(size_t num);
    virtual void set_init_worker_num(size_t num);
    virtual void set_job_assigner_thread_num(size_t num);
    virtual void set_command_batch_size(size_t num);
    virtual void set_job_done_batch_size(size_t num);

    virtual void set_split_dimensions(const std::vector<size_t>& split);
    virtual void set_sub_split_dimensions(const std::vector<size_t>& split);
    virtual void set_global_region(const GeometricRegion& region);

    /* TODO(omidm): figure out what we want for these methods. */
    virtual void LoadClusterMap(std::string) {}
    virtual void DeleteWorker(SchedulerWorker * worker) {}
    virtual void AddWorker(SchedulerWorker * worker) {}
    virtual SchedulerWorker* GetWorker(int workerId) {return NULL;}

    virtual void PrintStats();

  protected:
    virtual void SchedulerCoreProcessor();

    virtual void RegisterInitialWorkers();

    virtual size_t RegisterPendingWorkers();

    virtual void QueryWorkerStats();

    virtual size_t AssignReadyJobs();

    virtual size_t RemoveObsoleteJobEntries();

    virtual void AddMainJob();

    virtual void TerminationProcedure();

    virtual size_t ProcessQueuedSchedulerCommands();

    virtual void ProcessSchedulerCommand(SchedulerCommand* cm);
    virtual void ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm);
    virtual void ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm);
    virtual void ProcessDefineDataCommand(DefineDataCommand* cm);
    virtual void ProcessHandshakeCommand(HandshakeCommand* cm);
    virtual void ProcessJobDoneCommand(JobDoneCommand* cm);
    virtual void ProcessDefinePartitionCommand(DefinePartitionCommand* cm);
    virtual void ProcessTerminateCommand(TerminateCommand* cm);
    virtual void ProcessSpawnTemplateCommand(SpawnTemplateCommand* cm);
    virtual void ProcessStartTemplateCommand(StartTemplateCommand* cm);
    virtual void ProcessEndTemplateCommand(EndTemplateCommand* cm);
    virtual void ProcessSaveDataJobDoneCommand(SaveDataJobDoneCommand* cm);
    virtual void ProcessWorkerDownCommand(WorkerDownCommand* cm);
    virtual void ProcessPrepareRewindCommand(PrepareRewindCommand* cm);
    virtual void ProcessRespondStatCommand(RespondStatCommand* cm);

    virtual void WaitForAllPrepareRewindResponses();

    virtual void CreateIDMaker();
    virtual void CreateAfterMap();
    virtual void CreateSchedulerServer();
    virtual void CreateDataManager();
    virtual void CreateJobManager();
    virtual void CreateTemplateManager();
    virtual void CreateLoadBalancer();
    virtual void CreateJobAssigner();

    virtual void SetupIDMaker();
    virtual void SetupSchedulerServer();
    virtual void SetupDataManager();
    virtual void SetupJobManager();
    virtual void SetupTemplateManager();
    virtual void SetupLoadBalancer();
    virtual void SetupJobAssigner();
    virtual void SetupCleaner();
    virtual void SetupJobDoneBouncer();
    virtual void SetupUserInterface();

    virtual void CleanerThread();
    virtual void JobDoneBouncerThread();
    virtual void GetUserCommandThread();

    virtual void LoadUserCommands();
    virtual void LoadWorkerCommands();

    IDMaker *id_maker_;
    AfterMap *after_map_;
    SchedulerServer* server_;
    JobManager* job_manager_;
    TemplateManager* template_manager_;
    DataManager* data_manager_;
    JobAssigner *job_assigner_;
    LoadBalancer *load_balancer_;

    std::map<job_id_t, std::string> template_spawner_map_;

    boost::thread* cleaner_thread_;
    boost::thread* user_interface_thread_;
    boost::thread* job_done_bouncer_thread_;
    boost::thread* scheduler_server_thread_;
    boost::thread* load_balancer_thread_;

    bool cleaner_thread_active_;
    bool bouncer_thread_active_;

    bool controller_template_active_;
    bool complex_memoization_active_;
    bool binding_memoization_active_;
    bool cascaded_binding_active_;
    bool worker_template_active_;
    bool mega_rcr_job_active_;
    bool data_manager_query_cache_active_;

    bool load_balancing_active_;
    int64_t load_balancing_period_;
    bool fault_tolerance_active_;
    int64_t checkpoint_creation_period_;

    std::vector<size_t> split_;
    std::vector<size_t> sub_split_;
    GeometricRegion global_region_;

    port_t listening_port_;
    size_t init_worker_num_;
    size_t command_batch_size_;
    size_t assign_batch_size_;
    size_t remove_batch_size_;
    size_t job_assigner_thread_num_;
    size_t job_done_batch_size_;

    counter_t query_stat_id_;
    int64_t last_query_stat_time_;
    size_t responded_worker_num_;
    size_t registered_worker_num_;
    bool terminate_application_flag_;
    exit_status_t terminate_application_status_;

    CmSet user_command_set_;
    SchedulerCommand::PrototypeTable worker_command_table_;

    Log log_;
    Log log_process_;
    Log log_overhead_;
    size_t processed_command_num_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_SCHEDULER_H_
