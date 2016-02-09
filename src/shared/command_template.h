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
  * This is CommandTemplate module to hold and instantiate a set of commands
  * sent from controller to worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_COMMAND_TEMPLATE_H_
#define NIMBUS_SRC_SHARED_COMMAND_TEMPLATE_H_

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <set>
#include "src/shared/nimbus_types.h"
#include "src/shared/scheduler_command_include.h"
#include "src/shared/dbg.h"
#include "src/shared/log.h"

namespace nimbus {

class SchedulerClient;

class CommandTemplate {
  public:
    CommandTemplate(const std::string& command_template_name,
                    const std::vector<job_id_t>& inner_job_ids,
                    const std::vector<job_id_t>& outer_job_ids,
                    const std::vector<physical_data_id_t>& phy_ids);

    ~CommandTemplate();


    bool finalized();
    size_t copy_job_num();
    size_t compute_job_num();
    std::string command_template_name();


    bool Finalize();

    bool Instantiate(const std::vector<job_id_t>& inner_job_ids,
                     const std::vector<job_id_t>& outer_job_ids,
                     const std::vector<Parameter>& parameters,
                     const std::vector<physical_data_id_t>& physical_ids,
                     SchedulerClient *client);

    bool AddComputeJobCommand(ComputeJobCommand* command);

    bool AddLocalCopyCommand(LocalCopyCommand* command);

    bool AddRemoteCopySendCommand(RemoteCopySendCommand* command);

    bool AddRemoteCopyReceiveCommand(RemoteCopyReceiveCommand* command);

    bool AddMegaRCRCommand(MegaRCRCommand* command);

  private:
    typedef boost::shared_ptr<job_id_t> JobIdPtr;
    typedef std::vector<JobIdPtr> JobIdPtrList;
    typedef boost::unordered_set<JobIdPtr> JobIdPtrSet;
    typedef boost::unordered_map<job_id_t, JobIdPtr> JobIdPtrMap;

    typedef boost::shared_ptr<physical_data_id_t> PhyIdPtr;
    typedef std::vector<PhyIdPtr> PhyIdPtrList;
    typedef boost::unordered_set<PhyIdPtr> PhyIdPtrSet;
    typedef boost::unordered_map<job_id_t, PhyIdPtr> PhyIdPtrMap;


    enum CommandTemplateType {
      BASE,
      COMPUTE,
      LC,
      RCS,
      RCR,
      MEGA_RCR
    };

    class BaseCommandTemplate {
      public:
        BaseCommandTemplate() {type_ = BASE;}
        ~BaseCommandTemplate() {}

        CommandTemplateType type_;
    };

    typedef std::vector<BaseCommandTemplate*> CommandTemplateVector;

    class ComputeJobCommandTemplate : public BaseCommandTemplate {
      public:
        ComputeJobCommandTemplate(const std::string& job_name,
                                  JobIdPtr job_id_ptr,
                                  const PhyIdPtrSet& read_set_ptr,
                                  const PhyIdPtrSet& write_set_ptr,
                                  const JobIdPtrSet& before_set_ptr,
                                  const JobIdPtrSet& after_set_ptr,
                                  JobIdPtr future_job_id_ptr,
                                  const bool& sterile,
                                  const GeometricRegion& region,
                                  const worker_id_t& worker_id)
          : job_name_(job_name),
            job_id_ptr_(job_id_ptr),
            read_set_ptr_(read_set_ptr),
            write_set_ptr_(write_set_ptr),
            before_set_ptr_(before_set_ptr),
            after_set_ptr_(after_set_ptr),
            future_job_id_ptr_(future_job_id_ptr),
            sterile_(sterile),
            region_(region),
            worker_id_(worker_id) {type_ = COMPUTE;}

        ~ComputeJobCommandTemplate() {}

        std::string job_name_;
        JobIdPtr job_id_ptr_;
        PhyIdPtrSet read_set_ptr_;
        PhyIdPtrSet write_set_ptr_;
        JobIdPtrSet before_set_ptr_;
        JobIdPtrSet after_set_ptr_;
        JobIdPtr future_job_id_ptr_;
        bool sterile_;
        GeometricRegion region_;
        worker_id_t worker_id_;
        size_t param_index_;
    };

    class LocalCopyCommandTemplate : public BaseCommandTemplate {
      public:
        LocalCopyCommandTemplate(JobIdPtr job_id_ptr,
                                 PhyIdPtr from_physical_data_id_ptr,
                                 PhyIdPtr to_physical_data_id_ptr,
                                 const JobIdPtrSet& before_set_ptr,
                                 const worker_id_t& worker_id)
          : job_id_ptr_(job_id_ptr),
            from_physical_data_id_ptr_(from_physical_data_id_ptr),
            to_physical_data_id_ptr_(to_physical_data_id_ptr),
            before_set_ptr_(before_set_ptr),
            worker_id_(worker_id) {type_ = LC;}

        ~LocalCopyCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        PhyIdPtr from_physical_data_id_ptr_;
        PhyIdPtr to_physical_data_id_ptr_;
        JobIdPtrSet before_set_ptr_;
        worker_id_t worker_id_;
    };


    class RemoteCopySendCommandTemplate : public BaseCommandTemplate {
      public:
        RemoteCopySendCommandTemplate(JobIdPtr job_id_ptr,
                                      JobIdPtr receive_job_id_ptr,
                                      JobIdPtr mega_rcr_job_id_ptr,
                                      PhyIdPtr from_physical_data_id_ptr,
                                      const ID<worker_id_t>& to_worker_id,
                                      const std::string to_ip,
                                      const ID<port_t>& to_port,
                                      const JobIdPtrSet& before_set_ptr,
                                      const worker_id_t& worker_id)
          : job_id_ptr_(job_id_ptr),
            receive_job_id_ptr_(receive_job_id_ptr),
            mega_rcr_job_id_ptr_(mega_rcr_job_id_ptr),
            from_physical_data_id_ptr_(from_physical_data_id_ptr),
            to_worker_id_(to_worker_id),
            to_ip_(to_ip),
            to_port_(to_port),
            before_set_ptr_(before_set_ptr),
            worker_id_(worker_id) {type_ = RCS;}

        ~RemoteCopySendCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        JobIdPtr receive_job_id_ptr_;
        JobIdPtr mega_rcr_job_id_ptr_;
        PhyIdPtr from_physical_data_id_ptr_;
        ID<worker_id_t> to_worker_id_;
        std::string to_ip_;
        ID<port_t> to_port_;
        JobIdPtrSet before_set_ptr_;
        worker_id_t worker_id_;
    };

    class RemoteCopyReceiveCommandTemplate : public BaseCommandTemplate {
      public:
        RemoteCopyReceiveCommandTemplate(JobIdPtr job_id_ptr,
                                         PhyIdPtr to_physical_data_id_ptr,
                                         const JobIdPtrSet& before_set_ptr,
                                         const worker_id_t& worker_id)
          : job_id_ptr_(job_id_ptr),
            to_physical_data_id_ptr_(to_physical_data_id_ptr),
            before_set_ptr_(before_set_ptr),
            worker_id_(worker_id) {type_ = RCR;}

        ~RemoteCopyReceiveCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        PhyIdPtr to_physical_data_id_ptr_;
        JobIdPtrSet before_set_ptr_;
        worker_id_t worker_id_;
    };

    class MegaRCRCommandTemplate : public BaseCommandTemplate {
      public:
        MegaRCRCommandTemplate(JobIdPtr job_id_ptr,
                               const JobIdPtrList& receive_job_id_ptrs,
                               const PhyIdPtrList& to_phy_id_ptrs,
                               const worker_id_t& worker_id)
          : job_id_ptr_(job_id_ptr),
            receive_job_id_ptrs_(receive_job_id_ptrs),
            to_phy_id_ptrs_(to_phy_id_ptrs),
            worker_id_(worker_id) {type_ = MEGA_RCR;}

        ~MegaRCRCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        JobIdPtrList receive_job_id_ptrs_;
        PhyIdPtrList to_phy_id_ptrs_;
        worker_id_t worker_id_;
    };

    bool finalized_;
    size_t compute_job_num_;
    size_t copy_job_num_;
    std::string command_template_name_;
    // Currently we do not support future job - omidm
    JobIdPtr future_job_id_ptr_;

    PhyIdPtrMap phy_id_map_;
    PhyIdPtrList phy_id_list_;

    JobIdPtrMap inner_job_id_map_;
    JobIdPtrList inner_job_id_list_;

    JobIdPtrMap outer_job_id_map_;
    JobIdPtrList outer_job_id_list_;

    CommandTemplateVector command_templates_;

    mutable boost::mutex mutex_;


    JobIdPtr GetExistingInnerJobIdPtr(job_id_t job_id);

    PhyIdPtr GetExistingPhyIdPtr(physical_data_id_t pdid);


    void PushComputeJobCommand(ComputeJobCommandTemplate* command,
                               const Parameter& parameter,
                               SchedulerClient *client);

    void PushLocalCopyCommand(LocalCopyCommandTemplate* command,
                              SchedulerClient *client);

    void PushRemoteCopySendCommand(RemoteCopySendCommandTemplate* command,
                                   SchedulerClient *client);

    void PushRemoteCopyReceiveCommand(RemoteCopyReceiveCommandTemplate* command,
                                      SchedulerClient *client);

    void PushMegaRCRCommand(MegaRCRCommandTemplate* command,
                            SchedulerClient *client);
};
}  // namespace nimbus
#endif  // NIMBUS_SRC_SHARED_COMMAND_TEMPLATE_H_
