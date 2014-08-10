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
 * Author: Andrew Lim <alim16@stanford.edu>
 */

#include <shared/profiler.h>

#ifdef __MACH__
#include <sys/param.h>
#include <sys/mount.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#endif

#include <boost/asio.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <cassert>
#include <fstream>  // NOLINT
#include "boost/lexical_cast.hpp"
#include "shared/profile_command.h"
#include "shared/profiler_malloc.h"

#define PERIOD 200

namespace nimbus {

  void Profiler::Run(SchedulerClient* client_, worker_id_t worker_id_) {
    ID<worker_id_t> worker_id = ID<worker_id_t>(worker_id_);
    boost::asio::io_service io;
    while (true) {
      boost::asio::deadline_timer t(io, boost::posix_time::milliseconds(PERIOD));
      Profile();
      ProfileCommand cm(worker_id, totalVirtualMem, usedVirtualMem,
          procVirtualMem, totalPhysMem, usedPhysMem, procPhysMem);
      client_->SendCommand(&cm);
      boost::this_thread::interruption_point();
      t.wait();
    }
  }

  void Profiler::Profile() {
#ifndef __MACH__
    sysinfo(&mem_info_);
#endif
    UpdateTotalVirtualMemory();
    UpdateUsedVirtualMemory();
    UpdateProcVirtualMemory();
    UpdateTotalPhysicalMemory();
    UpdateUsedPhysicalMemory();
    UpdateProcPhysicalMemory();
  }

  std::string Profiler::ToNetworkData() {
    std::string str;
    str += (boost::lexical_cast<std::string>(totalVirtualMem) + " ");
    str += (boost::lexical_cast<std::string>(usedVirtualMem) + " ");
    str += (boost::lexical_cast<std::string>(procVirtualMem) + " ");
    str += (boost::lexical_cast<std::string>(totalPhysMem) + " ");
    str += (boost::lexical_cast<std::string>(usedPhysMem) + " ");
    str += (boost::lexical_cast<std::string>(procPhysMem));
    return str;
  }

  std::string Profiler::ToString() {
    std::string str;
    str += ("total_vm: " + boost::lexical_cast<std::string>(totalVirtualMem) + " ");
    str += ("used_vm: " + boost::lexical_cast<std::string>(usedVirtualMem) + " ");
    str += ("proc_vm: " + boost::lexical_cast<std::string>(procVirtualMem) + " ");
    str += ("total_pm: " + boost::lexical_cast<std::string>(totalPhysMem) + " ");
    str += ("used_pm: " + boost::lexical_cast<std::string>(usedPhysMem) + " ");
    str += ("proc_pm: " + boost::lexical_cast<std::string>(procPhysMem));
    return str;
  }

  void Profiler::UpdateTotalVirtualMemory() {
#ifdef __MACH__
    // totalVirtualMem = 0;
#else
    totalVirtualMem = mem_info_.totalram;
    totalVirtualMem += mem_info_.totalswap;
    totalVirtualMem *= mem_info_.mem_unit;
#endif
  }

  void Profiler::UpdateUsedVirtualMemory() {
#ifdef __MACH__
    // usedVirtualMem = 0;
#else
    usedVirtualMem = mem_info_.totalram - mem_info_.freeram;
    usedVirtualMem += mem_info_.totalswap - mem_info_.freeswap;
    usedVirtualMem *= mem_info_.mem_unit;
#endif
  }

  void Profiler::UpdateProcVirtualMemory() {
#ifdef __MACH__
    task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    if (KERN_SUCCESS == task_info(mach_task_self(),
                                  TASK_BASIC_INFO, (task_info_t)&t_info,
                                  &t_info_count)) {
      procVirtualMem = t_info.virtual_size;
    }
#else
    std::ifstream ifs;
    ifs.open("/proc/self/status");
    std::string line;
    while (getline(ifs, line)) {
      if (line.compare(0, 7, "VmSize:") == 0) {
        procVirtualMem = ParseLine(line) * 1024;
        break;
      }
    }
    ifs.close();
#endif
  }

  void Profiler::UpdateTotalPhysicalMemory() {
#ifdef __MACH__
    int mib[2];
    size_t length;

    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    length = sizeof(uint64_t);
    sysctl(mib, 2, &totalPhysMem, &length, NULL, 0);
#else
    totalPhysMem = mem_info_.totalram;
    totalPhysMem *= mem_info_.mem_unit;
#endif
  }

  void Profiler::UpdateUsedPhysicalMemory() {
#ifdef __MACH__
    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics_data_t vm_stats;

    mach_port = mach_host_self();
    count = sizeof(vm_stats) / sizeof(natural_t);
    if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
        KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO,
                                       (host_info_t)&vm_stats, &count)) {
      usedPhysMem = ((uint64_t)vm_stats.active_count +
                    (uint64_t)vm_stats.inactive_count +
                    (uint64_t)vm_stats.wire_count) * (uint64_t)page_size;
    }
#else
    usedPhysMem = mem_info_.totalram - mem_info_.freeram;
    usedPhysMem *= mem_info_.mem_unit;
#endif
  }

  void Profiler::UpdateProcPhysicalMemory() {
#ifdef __MACH__
    task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    if (KERN_SUCCESS == task_info(mach_task_self(),
                                  TASK_BASIC_INFO, (task_info_t)&t_info,
                                  &t_info_count)) {
      procPhysMem = t_info.resident_size;
    }
#else
    std::ifstream ifs;
    ifs.open("/proc/self/status");
    std::string line;
    while (getline(ifs, line)) {
      if (line.compare(0, 6, "VmRSS:") == 0) {
        procPhysMem = ParseLine(line) * 1024;
        break;
      }
    }
    ifs.close();
#endif
  }

  uint64_t Profiler::ParseLine(std::string line) {
    char *str = const_cast<char *>(line.c_str());
    int i = strlen(str);
    while (*str < '0' || *str > '9') str++;
    str[i-3] = '\0';
    i = atoi(str);
    return i;
  }

}  // namespace nimbus
