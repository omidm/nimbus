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
#include "shared/dbg.h"
#include "shared/log.h"

namespace nimbus {

class JobManager;

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

    bool AddComputeJob() {
      return false;
    }

    bool AddLocalCopyJob() {
      return false;
    }

    bool AddRemoteCopySendJob() {
      return false;
    }

    bool AddRemoteCopyReceiveJob() {
      return false;
    }

  private:
    typedef std::vector<boost::shared_ptr<job_id_t> > PtrList;
    typedef boost::unordered_set<boost::shared_ptr<job_id_t> > PtrSet;
    typedef boost::unordered_map<job_id_t, boost::shared_ptr<job_id_t> > PtrMap;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_BINDING_TEMPLATE_H_
