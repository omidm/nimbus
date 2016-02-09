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
 * Memory profiler.
 *      
 * Author: Andrew Lim <alim16@stanford.edu>
 */

#ifndef NIMBUS_SRC_SHARED_PROFILER_H_
#define NIMBUS_SRC_SHARED_PROFILER_H_

#ifndef __MACH__
#include <sys/sysinfo.h>
#endif

#include <boost/asio.hpp>
#include <inttypes.h>
#include <string>
#include "src/shared/scheduler_client.h"

namespace nimbus {

class Profiler {
 public:
  void Run(SchedulerClient* client_, worker_id_t worker_id);
  void Profile();
  std::string ToNetworkData();
  std::string ToString();

 private:
#ifndef __MACH__
  struct sysinfo mem_info_;
#endif
  uint64_t totalVirtualMem;
  uint64_t usedVirtualMem;
  uint64_t procVirtualMem;
  uint64_t totalPhysMem;
  uint64_t usedPhysMem;
  uint64_t procPhysMem;

 private:
  void UpdateTotalVirtualMemory();
  void UpdateUsedVirtualMemory();
  void UpdateProcVirtualMemory();
  void UpdateTotalPhysicalMemory();
  void UpdateUsedPhysicalMemory();
  void UpdateProcPhysicalMemory();
  uint64_t ParseLine(std::string line);
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_SHARED_PROFILER_H_
