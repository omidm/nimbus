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
  * Testing command parsing and communication between a scheduler and
  * a worker.
  * 
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_TEST_MESSAGING_CREATE_COMMANDS_H_
#define NIMBUS_TEST_MESSAGING_CREATE_COMMANDS_H_
using namespace nimbus;  // NOLINT

#include <string>
#include "shared/scheduler_command.h"

#define NUM_COMMANDS 50

void load_command_prototypes();
void create_commands();

nimbus::SchedulerCommand::PrototypeTable prototypes;
nimbus::SchedulerCommand* commands[NUM_COMMANDS*2];

void load_command_prototypes() {
  prototypes[SchedulerPBuf_Type_SPAWN_COMPUTE] = new SpawnComputeJobCommand();
  prototypes[SchedulerPBuf_Type_SPAWN_COPY] = new SpawnCopyJobCommand();
  prototypes[SchedulerPBuf_Type_DEFINE_DATA] = new DefineDataCommand();
  prototypes[SchedulerPBuf_Type_HANDSHAKE] = new HandshakeCommand();
  prototypes[SchedulerPBuf_Type_JOB_DONE] = new JobDoneCommand();
  prototypes[SchedulerPBuf_Type_EXECUTE_COMPUTE] = new ComputeJobCommand();
  prototypes[SchedulerPBuf_Type_CREATE_DATA] = new CreateDataCommand();
  prototypes[SchedulerPBuf_Type_REMOTE_SEND] = new RemoteCopySendCommand();
  prototypes[SchedulerPBuf_Type_REMOTE_RECEIVE] = new RemoteCopyReceiveCommand(); // NOLINT
  prototypes[SchedulerPBuf_Type_LOCAL_COPY] = new LocalCopyCommand();
  prototypes[SchedulerPBuf_Type_DEFINE_PARTITION] = new DefinePartitionCommand();
  prototypes[SchedulerPBuf_Type_LDO_ADD] = new LdoAddCommand();
  prototypes[SchedulerPBuf_Type_LDO_REMOVE] = new LdoRemoveCommand();
  prototypes[SchedulerPBuf_Type_PARTITION_ADD] = new PartitionAddCommand();
  prototypes[SchedulerPBuf_Type_PARTITION_REMOVE] = new PartitionRemoveCommand();
  prototypes[SchedulerPBuf_Type_TERMINATE] = new TerminateCommand();
  prototypes[SchedulerPBuf_Type_PROFILE] = new ProfileCommand();
}

IDSet<logical_data_id_t> ldo_set_a;
IDSet<logical_data_id_t> ldo_set_b;
IDSet<job_id_t> job_set_a;
IDSet<job_id_t> job_set_b;
ID<logical_data_id_t> ldo_a;
ID<logical_data_id_t> ldo_b;
ID<physical_data_id_t> phy_a;
ID<physical_data_id_t> phy_b;
IDSet<physical_data_id_t> phy_set_a;
IDSet<physical_data_id_t> phy_set_b;
ID<partition_id_t> partition;
IDSet<partition_id_t> partition_set;
ID<job_id_t> job_id;
ID<job_id_t> parent;
ID<job_id_t> future;
ID<worker_id_t> worker;
ID<port_t> port;
ID<exit_status_t> exit_val;
bool sterile;
Parameter params;
std::string name = "test-string";
std::string ip = "127.0.0.4";
GeometricRegion region(1, 2, 3, 101, 102, 103);
LogicalDataObject* logical_obj;

void initalize_sets() {
  for (int i = 0; i < 12; i++) {
    ldo_set_a.insert(10000 + i);
    phy_set_a.insert(40000 + i);
  }
  for (int i = 0; i < 11; i++) {
    ldo_set_b.insert(11000 + i);
    phy_set_b.insert(41000 + i);
  }
  for (int i = 0; i < 5; i++) {
    job_set_a.insert(20000 + i);
  }
  for (int i = 0; i < 7; i++) {
    job_set_b.insert(21000 + i);
  }
  for (int i = 0; i < 9; i++) {
    partition_set.insert(31000 + i);
  }

  ldo_a.set_elem(10125);
  ldo_b.set_elem(11119);
  phy_a.set_elem(40125);
  phy_b.set_elem(41119);
  job_id.set_elem(2014);
  parent.set_elem(1940);
  future.set_elem(1977);
  partition.set_elem(1947);
  worker.set_elem(22577);
  port.set_elem(39);
  exit_val.set_elem(666);
  sterile = false;

  SerializedData d("serialized data");
  params.set_ser_data(d);

  logical_obj = new LogicalDataObject(ldo_a.elem(), "velocity", &region);
}

void create_commands() {
  initalize_sets();

  const IDSet<logical_data_id_t> cldo_set_a = (const IDSet<logical_data_id_t>)ldo_set_a;
  const IDSet<logical_data_id_t> cldo_set_b = (const IDSet<logical_data_id_t>)ldo_set_b;
  const IDSet<job_id_t> cjob_set_a = (const IDSet<job_id_t>)job_set_a;
  const IDSet<job_id_t> cjob_set_b = (const IDSet<job_id_t>)job_set_b;
  const ID<logical_data_id_t> cldo_a = (const ID<logical_data_id_t>)ldo_a;
  const ID<logical_data_id_t> cldo_b = (const ID<logical_data_id_t>)ldo_b;
  const IDSet<physical_data_id_t> cphy_set_a = (const IDSet<physical_data_id_t>)phy_set_a;
  const IDSet<physical_data_id_t> cphy_set_b = (const IDSet<physical_data_id_t>)phy_set_b;
  const ID<physical_data_id_t> cphy_a = (const ID<physical_data_id_t>)phy_a;
  const ID<physical_data_id_t> cphy_b = (const ID<physical_data_id_t>)phy_b;
  const ID<partition_id_t> cpartition = (const ID<partition_id_t>)partition;
  const IDSet<partition_id_t> cpartition_set = (const IDSet<partition_id_t>)partition_set;
  const ID<job_id_t> cjob_id = (const ID<job_id_t>)job_id;
  const ID<job_id_t> cparent = (const ID<job_id_t>)parent;
  const ID<job_id_t> cfuture = (const ID<job_id_t>)future;
  const ID<worker_id_t> cworker = (const  ID<worker_id_t>)worker;
  const ID<port_t> cport = (const ID<port_t>)port;
  const ID<exit_status_t> cexit_val = (const ID<exit_status_t>)exit_val;
  const bool csterile = (const bool)sterile;
  const Parameter cparams = (const Parameter)params;
  const std::string cname = (const std::string)name;
  const std::string cip = (const std::string)ip;
  const GeometricRegion* cregion = (const GeometricRegion*)&region;
  const LogicalDataObject* clogical_obj = (const LogicalDataObject*)logical_obj;

  for (int i = 0; i < NUM_COMMANDS;) {
    commands[i++] = new SpawnComputeJobCommand(cname,
                                               cjob_id,
                                               cldo_set_a,
                                               cldo_set_b,
                                               cjob_set_a,
                                               cjob_set_b,
                                               cparent,
                                               cfuture,
                                               csterile,
                                               cparams);
    commands[i++] = new SpawnCopyJobCommand(cjob_id,
                                            cldo_a,
                                            cldo_b,
                                            cjob_set_a,
                                            cjob_set_b,
                                            cparent);
    commands[i++] = new DefineDataCommand(cname,
                                          cldo_a,
                                          cpartition,
                                          cpartition_set,
                                          cjob_id);
    commands[i++] = new HandshakeCommand(cworker,
                                         cip,
                                         cport);
    commands[i++] = new JobDoneCommand(cjob_id,
                                       109.2451,
                                       3.141592654,
                                       8192);
    commands[i++] = new ComputeJobCommand(cname,
                                          cjob_id,
                                          cphy_set_a,
                                          cphy_set_b,
                                          cjob_set_a,
                                          cjob_set_b,
                                          cfuture,
                                          csterile,
                                          cparams);

    commands[i++] = new CreateDataCommand(cjob_id, cname, cldo_a, cphy_a, cjob_set_a);
    commands[i++] = new RemoteCopySendCommand(cjob_id, cjob_id, cphy_a, cworker,
                                              cip, cport, cjob_set_a);
    commands[i++] = new RemoteCopyReceiveCommand(cjob_id, cphy_a, cjob_set_a);
    commands[i++] = new LocalCopyCommand(cjob_id, cphy_a, cphy_b, cjob_set_a);
    commands[i++] = new DefinePartitionCommand(cpartition, *cregion);
    commands[i++] = new LdoAddCommand(clogical_obj);
    commands[i++] = new LdoRemoveCommand(logical_obj);
    commands[i++] = new PartitionAddCommand(cpartition, *cregion);
    commands[i++] = new PartitionRemoveCommand(cpartition);
    commands[i++] = new TerminateCommand(cexit_val);
    commands[i++] = new ProfileCommand(cworker, 20, 40, 60, 21, 41, 61);
  }
}

#endif  // NIMBUS_TEST_MESSAGING_CREATE_COMMANDS_H_
