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
  * This is BindingTemplate module to hold and instantiate copy and compute
  * commands sent to workers to bind physical data to a batch of jobs in a
  * template complex job.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_BINDING_TEMPLATE_H_
#define NIMBUS_SCHEDULER_BINDING_TEMPLATE_H_

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
#include "shared/nimbus_types.h"
#include "shared/scheduler_command_include.h"
#include "shared/dbg.h"
#include "shared/log.h"

namespace nimbus {

class BindingTemplate {
  public:
    BindingTemplate();
    ~BindingTemplate();

    bool Finalize() {
      return false;
    }

    bool Instantiate() {
      return false;
    }

    enum VERSION_TYPE {
      REGULAR,
      WILD_CARD
    };

    bool TrackDataObject(const logical_data_id_t& ldid,
                         const physical_data_id_t& pdid,
                         VERSION_TYPE version_type,
                         data_version_t version_diff_from_base = 0);

    bool AddComputeJobCommand(const ComputeJobCommand& command,
                              worker_id_t w_id);

    bool AddLocalCopyCommand(const LocalCopyCommand& command,
                             worker_id_t w_id);

    bool AddRemoteCopySendCommand(const RemoteCopySendCommand& command,
                                  worker_id_t w_id);

    bool AddRemoteCopyReceiveCommand(const RemoteCopyReceiveCommand& command,
                                     worker_id_t w_id);

  private:
    typedef boost::shared_ptr<job_id_t> JobIdPtr;
    typedef std::vector<JobIdPtr> JobIdPtrList;
    typedef boost::unordered_set<JobIdPtr> JobIdPtrSet;
    typedef boost::unordered_map<job_id_t, JobIdPtr> JobIdPtrMap;

    typedef boost::shared_ptr<physical_data_id_t> PhyIdPtr;
    typedef std::vector<PhyIdPtr> PhyIdPtrList;
    typedef boost::unordered_set<PhyIdPtr> PhyIdPtrSet;
    typedef boost::unordered_map<job_id_t, PhyIdPtr> PhyIdPtrMap;

    class ComputeJobCommandTemplate {
      public:
        ComputeJobCommandTemplate(const std::string& job_name,
                                  JobIdPtr job_id_ptr,
                                  const PhyIdPtrSet& read_set_ptr,
                                  const PhyIdPtrSet& write_set_ptr,
                                  const JobIdPtrSet& before_set_ptr,
                                  const JobIdPtrSet& after_set_ptr,
                                  JobIdPtr future_job_id_ptr,
                                  const bool& sterile,
                                  const GeometricRegion& region)
          : job_name_(job_name),
            job_id_ptr_(job_id_ptr),
            read_set_ptr_(read_set_ptr),
            write_set_ptr_(write_set_ptr),
            before_set_ptr_(before_set_ptr),
            after_set_ptr_(after_set_ptr),
            future_job_id_ptr_(future_job_id_ptr),
            sterile_(sterile),
            region_(region) {}

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
        Parameter params_;
    };

    class LocalCopyCommandTemplate {
      public:
        LocalCopyCommandTemplate(JobIdPtr job_id_ptr,
                                 PhyIdPtr from_physical_data_id_ptr,
                                 PhyIdPtr to_physical_data_id_ptr,
                                 const JobIdPtrSet& before_set_ptr)
          : job_id_ptr_(job_id_ptr),
            from_physical_data_id_ptr_(from_physical_data_id_ptr),
            to_physical_data_id_ptr_(to_physical_data_id_ptr),
            before_set_ptr_(before_set_ptr) {}

        ~LocalCopyCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        PhyIdPtr from_physical_data_id_ptr_;
        PhyIdPtr to_physical_data_id_ptr_;
        JobIdPtrSet before_set_ptr_;
    };


    class RemoteCopySendCommandTemplate {
      public:
        RemoteCopySendCommandTemplate(JobIdPtr job_id_ptr,
                                      JobIdPtr receive_job_id_ptr,
                                      PhyIdPtr from_physical_data_id_ptr,
                                      const ID<worker_id_t>& to_worker_id,
                                      const std::string to_ip,
                                      const ID<port_t>& to_port,
                                      const JobIdPtrSet& before_set_ptr)
          : job_id_ptr_(job_id_ptr),
            receive_job_id_ptr_(receive_job_id_ptr),
            from_physical_data_id_ptr_(from_physical_data_id_ptr),
            to_worker_id_(to_worker_id),
            to_ip_(to_ip),
            to_port_(to_port),
            before_set_ptr_(before_set_ptr) {}

        ~RemoteCopySendCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        JobIdPtr receive_job_id_ptr_;
        PhyIdPtr from_physical_data_id_ptr_;
        ID<worker_id_t> to_worker_id_;
        std::string to_ip_;
        ID<port_t> to_port_;
        JobIdPtrSet before_set_ptr_;
    };

    class RemoteCopyReceiveCommandTemplate {
      public:
        RemoteCopyReceiveCommandTemplate(JobIdPtr job_id_ptr,
                                         PhyIdPtr to_physical_data_id_ptr,
                                         const JobIdPtrSet& before_set_ptr)
          : job_id_ptr_(job_id_ptr),
            to_physical_data_id_ptr_(to_physical_data_id_ptr),
            before_set_ptr_(before_set_ptr) {}

        ~RemoteCopyReceiveCommandTemplate() {}

        JobIdPtr job_id_ptr_;
        PhyIdPtr to_physical_data_id_ptr_;
        JobIdPtrSet before_set_ptr_;
    };












    bool finalized_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_BINDING_TEMPLATE_H_
