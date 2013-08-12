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
#include "./simple_scheduler.h"
#include "lib/nimbus.h"
#include "lib/command.h"
#include "lib/parser.h"
int main(int argc, char *argv[]) {
  nimbus::nimbus_initialize();

  std::string str = "createjob name:main id:{0} read:{1,2} write:{1,2} ";
  str += " before:{} after:{1,2,3} type:operation param:t=20,g=6";
  Command cm(str);
  std::cout << cm.toString() << std::endl;

  Log log;
  log.writeToBuffer("**Start of the log file.");
  log.dbg_writeToBuffer("Some DEBUG information in the buffer!", LOG_DEBUG);
  log.dbg_writeToBuffer("Some more DEBUG information in the buffer!",
                        LOG_DEBUG);
  log.writeBufferToFile();
  Log::printLine("Nimbus is up!", LOG_INFO);
  Log::dbg_printLine("DEBUG information will be printed!", LOG_DEBUG);

  SimpleScheduler * s = new SimpleScheduler(NIMBUS_SCHEDULER_PORT);
  s->run();
}

