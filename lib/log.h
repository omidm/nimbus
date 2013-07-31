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
  * Nimbus log interface. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_LIB_LOG_H_
#define NIMBUS_LIB_LOG_H_

#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream> // NOLINT
#include <string>

enum LOG_TYPE {
  ERROR,
  WARNING,
  INFO,
  DEBUG,
  NONE
};

std::string getTag(LOG_TYPE type);

class Log {
  public:
    Log();
    explicit Log(std::ostream* os);
    explicit Log(std::string fname);
    Log(std::ostream* os, std::string fname);
    ~Log();


    void setOutputStream(std::ostream* os);

    void setFileName(std::string fname);

    void clearBuffer();

    void clearLogFile();

    void writeToBuffer(std::string buf,
        LOG_TYPE type = NONE, bool flag = true);

    void writeToFile(std::string buf,
        LOG_TYPE type = NONE, bool flag = true);

    void writeToOutputStream(std::string buf,
        LOG_TYPE type = NONE, bool flag = true);

    void writeBufferToFile(bool flag = true);

    void writeBufferToOutputStream(bool flag = true);

    static void printLine(std::string msg,
        LOG_TYPE type = NONE, bool flag = true) {
      if (flag)
        std::cout << getTag(type) << msg << std::endl;
    };

    static void print(std::string msg,
        LOG_TYPE type = NONE, bool flag = true) {
      if (flag)
        std::cout << getTag(type) << msg;
    };

  private:
    std::ostream* output_stream;
    std::stringstream buffer;
    std::string log_file_name;
};











#endif  // NIMBUS_LIB_LOG_H_
