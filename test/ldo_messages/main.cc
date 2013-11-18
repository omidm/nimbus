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
  * This file tests whether LogicalDataObjects are being serialized
  * and deserialized correctly.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */


#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"
#include "shared/scheduler_command.h"
#include "shared/ldo_add_command.h"
#include "shared/ldo_remove_command.h"
#include "shared/dbg.h"

using namespace nimbus;  // NOLINT

void printLdo(nimbus::LogicalDataObject* obj) {
  printf("**Object - ID: %llu, Name: %s", obj->id(), obj->variable().c_str());
  printf(" region: [%llu+%llu, %llu+%llu, %llu+%llu]\n", obj->region()->x(), obj->region()->dx(), obj->region()->y(), obj->region()->dy(), obj->region()->z(), obj->region()->dz());  // NOLINT
}

int main(int argc, char *argv[]) {
  std::string variable = "pressure";
  GeometricRegion* region = new GeometricRegion(lrand48(), lrand48(),
                                                lrand48(), lrand48(),
                                                lrand48(), lrand48());
  LogicalDataObject* obj = new LogicalDataObject(4399, variable, region);

  LdoAddCommand* addStart = new LdoAddCommand(obj);
  std::string strVal = addStart->toString();
  std::string params;

  SchedulerCommand* command = new SchedulerCommand();
  SchedulerCommand::PrototypeTable* table = new SchedulerCommand::PrototypeTable();
  table->push_back(new LdoAddCommand());

  SchedulerCommand* c;
  bool result = command->GenerateSchedulerCommandChild(strVal,
                                                       table,
                                                       c);
  if (result == false) {
    printf("Failed to parse command type %s\n", strVal.c_str());
    return -1;
  }

  LdoAddCommand* addEnd = static_cast<LdoAddCommand*>(c);
  printf("Testing add command.\n");
  printf("Before: ");
  printLdo(obj);
  printf("After:  ");
  printLdo(addEnd->object());

  LdoRemoveCommand* removeStart = new LdoRemoveCommand(obj);
  strVal = removeStart->toString();

  table->push_back(new LdoRemoveCommand());

  result = command->GenerateSchedulerCommandChild(strVal, table, c);
  LdoRemoveCommand* removeEnd = static_cast<LdoRemoveCommand*>(c);
  printf("Testing remove command.\n");
  printf("Before: ");
  printLdo(obj);
  printf("After:  ");
  printLdo(removeEnd->object());
}
