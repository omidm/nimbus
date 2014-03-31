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
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <boost/thread/thread.hpp>
#include <iostream>  // NOLINT
#include "shared/dbg.h"
#include "shared/nimbus_types.h"
#include "shared/worker_data_exchanger.h"
#include "test/data_exchanger/config.h"

int main() {
  nimbus::WorkerDataExchanger data_exchanger(PORT_1);
  boost::thread data_exchanger_thread(
      boost::bind(&nimbus::WorkerDataExchanger::Run, &data_exchanger));
  data_exchanger.AddContactInfo(TO_WORKER_ID, TO_IP, PORT_2);
  nimbus::SerializedData ser_data;
  for (int index = 0; index < 20; ++index) {
    // Write the data.
    std::ostringstream convert;
    convert << index;      // insert the textual representation of 'Number' in
    ser_data.Parse(convert.str());
    data_exchanger.SendSerializedData(
        INIT_JOB_ID + index, TO_WORKER_ID, ser_data);
    std::cout << index << " sent" << std::endl;
  }
  return 0;
}

