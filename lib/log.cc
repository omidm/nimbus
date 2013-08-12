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

#include "lib/log.h"

Log::Log()
: output_stream(&std::cout),
  log_file_name("log.txt") {
    clearBuffer();
    clearLogFile();
}

Log::Log(std::ostream* os)
: output_stream(os),
  log_file_name("log.txt") {
    clearBuffer();
    clearLogFile();
}

Log::Log(std::string fname)
: output_stream(&std::cout),
  log_file_name(fname) {
    clearBuffer();
    clearLogFile();
}

Log::Log(std::ostream* os, std::string fname)
: output_stream(os),
  log_file_name(fname) {
    clearBuffer();
    clearLogFile();
}

Log::~Log() {
  clearBuffer();
  clearLogFile();
}


void Log::setOutputStream(std::ostream* os) {
  output_stream = os;
}

void Log::setFileName(std::string fname) {
  log_file_name = fname;
}

void Log::clearBuffer() {
  buffer.str("");
}

void Log::clearLogFile() {
  std::stringstream ss;
  ss.str("");
  std::ofstream ofs;
  ofs.open(log_file_name.c_str());
  ofs << ss.str();
  ofs.close();
}

void Log::writeToBuffer(std::string buf, LOG_TYPE type) {
  buffer << getTag(type) << buf << std::endl;
}

void Log::writeToFile(std::string buf, LOG_TYPE type) {
  std::ofstream ofs;
  ofs.open(log_file_name.c_str(), std::ofstream::app);
  ofs << getTag(type) << buf << std::endl;
  ofs.close();
}

void Log::writeToOutputStream(std::string buf, LOG_TYPE type) {
  buffer << getTag(type) << buf << std::endl;
}

void Log::writeBufferToFile() {
  std::ofstream ofs;
  ofs.open(log_file_name.c_str(), std::ofstream::app);
  ofs << buffer.str();
  ofs.close();
}


void Log::writeBufferToOutputStream() {
  *output_stream << buffer.str();
}

void Log::printLine(std::string msg, LOG_TYPE type) {
  std::cout << getTag(type) << msg << std::endl;
}

void Log::print(std::string msg, LOG_TYPE type) {
  std::cout << getTag(type) << msg;
}


std::string getTag(LOG_TYPE type) {
  switch (type) {
    case LOG_ERROR:
      return "ERROR: ";
    case LOG_WARNING:
      return "WARNING: ";
    case LOG_INFO:
      return "INFO: ";
    case LOG_DEBUG:
      return "DEBUG: ";
    case LOG_NONE:
      return "";
    default :
      return "";
  }
}

