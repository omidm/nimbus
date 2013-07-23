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
  * The most trivial test application.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include <pthread.h>
#include <iostream>  // NOLINT

#include "lib/scheduler_command.h"

using ::std::cout;
using ::std::endl;

const char* commands[] = {
      "no-op",
      "halt 53",
      "run job3 job0,job1,job2 job5,job6 data4,data5 data5 blah",
      "copy   data4          host34   ",
      "copy       data5 192.244.11.2        ",
      "",
      "newline\n\ntest",
      NULL
};

int main(int argc, char *argv[]) {
  std::cout << "Testing scheduler command class." << std::endl;
  int i = 0;
  while (commands[i] != NULL) {
    cout << "Testing command \'" << commands[i] << std::endl;
    SchedulerCommand* c = new SchedulerCommand(commands[i]);
    cout << "  translated to string \'" << c->toString() << '\'' << std::endl;
    cout << "  translated to tokens ";
    cout << c->getName() << ":";
    CommandParameterList params = c->getParameters();
    CommandParameterList::const_iterator iter = params.begin();
    for (; iter != params.end(); ++iter) {
      std::string param = *iter;
      cout << param << ";";
    }
    cout << std::endl;
    i++;
    delete c;
  }
}
