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
  * This file has the main function that launches Nimbus scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#define DEBUG_MODE

#include <iostream> // NOLINT
#include "./scheduler_multi.h"
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "shared/scheduler_command_include.h"
#include "shared/parser.h"
#include "shared/parameter.h"

int main(int argc, char *argv[]) {
  nimbus::nimbus_initialize();

//  std::string str = "createjob name:main id:{0} read:{1,2} write:{1,2} ";
//  str += " before:{} after:{1,2,3} type:operation param:t=20,g=6";
//  SchedulerCommand cm(str);
//  std::cout << cm.ToNetworkData() << std::endl;
//
//  str = "main {0} {1,2} {4,5} ";
//  str += " {6,7,8} {10,20,30} COMP t=20,g=6";
//  std::string job_name;
//  IDSet<job_id_t> job_id;
//  IDSet<data_id_t> read;
//  IDSet<data_id_t> write;
//  IDSet<job_id_t> before;
//  IDSet<job_id_t> after;
//  JobType job_type;
//  std::string params;
//
//  std::cout << "SpawnJob parameters string: " << str << std::endl;
//  bool cond = ParseSpawnJobCommand(str, job_name, job_id, read, write,
//      before, after, job_type, params);
//  if (cond) {
//    SpawnJobCommand sjc(job_name, job_id, read, write,
//      before, after, job_type, params);
//    std::cout << "Spawn Job correctly parsed as: " <<
//      sjc.ToString() << std::endl;
//  }
//
//  while (!cond) {}

//  std::string str("X\0;:,", 5);
//  SerializedData ser_data(str);
//  IDSet<param_id_t> idset;
//  idset.insert(13);
//  Parameter param(ser_data, idset);
//  std::string param_str = param.ToNetworkData();
//  std::cout << "param to string before: " << param_str << std::endl;
//
//  bool parsed = false;
//  SerializedData temp_ser_data;
//  IDSet<param_id_t> temp_idset;
//  parsed = ParseParameter(param_str, temp_ser_data, temp_idset);
//  Parameter temp_param(temp_ser_data, temp_idset);
//  std::cout << "parsed ser_data: " << temp_ser_data.ToNetworkData() << std::endl;
//  std::cout << "parsed idset: " << temp_idset.ToNetworkData() << std::endl;
//  std::cout << "param to string after: " << temp_param.ToNetworkData() << std::endl;
//  while (!parsed) {}


  Log log;
  log.writeToBuffer("**Start of the log file.");
  log.dbg_writeToBuffer("Some DEBUG information in the buffer!", LOG_DEBUG);
  log.dbg_writeToBuffer("Some more DEBUG information in the buffer!",
                        LOG_DEBUG);
  log.writeBufferToFile();
  Log::printLine("Nimbus is up!", LOG_INFO);
  Log::dbg_printLine("DEBUG information will be printed!", LOG_DEBUG);

  SimpleScheduler * s = new SimpleScheduler(NIMBUS_SCHEDULER_PORT);
  s->Run();
}

