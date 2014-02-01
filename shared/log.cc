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

#include "shared/log.h"

Log::Log()
: output_stream_(&std::cout),
  log_file_name_("log.txt") {
    ClearBuffer();
    ClearLogFile();
    InitTime();
    ResetTimer();
}

Log::Log(std::ostream* os)
: output_stream_(os),
  log_file_name_("log.txt") {
    ClearBuffer();
    ClearLogFile();
    InitTime();
    ResetTimer();
}

Log::Log(std::string fname)
: output_stream_(&std::cout),
  log_file_name_(fname) {
    ClearBuffer();
    ClearLogFile();
    InitTime();
    ResetTimer();
}

Log::Log(std::ostream* os, std::string fname)
: output_stream_(os),
  log_file_name_(fname) {
    ClearBuffer();
    ClearLogFile();
    InitTime();
    ResetTimer();
}

Log::~Log() {
}

void Log::set_output_stream(std::ostream* os) {
  output_stream_ = os;
}

void Log::set_file_name(std::string fname) {
  log_file_name_ = fname;
}

void Log::InitTime() {
  gettimeofday(&start_time_, NULL);
  timer_ = 0;
  timer_is_on_ = false;
}

struct timeval* Log::start_time() {
  return &start_time_;
}

void Log::set_start_time(struct timeval* time) {
  start_time_.tv_sec = time->tv_sec;
  start_time_.tv_usec = time->tv_usec;
}

double Log::timer() {
  if (timer_is_on_) {
    struct timeval t;
    gettimeofday(&t, NULL);
    double temp  = (static_cast<double>(t.tv_sec - timer_start_time_.tv_sec)) +
      .000001 * (static_cast<double>(t.tv_usec - timer_start_time_.tv_usec));
    return timer_ + temp;
  } else {
    return timer_;
  }
}

double Log::GetTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  double time  = (static_cast<double>(t.tv_sec - start_time_.tv_sec)) +
  .000001 * (static_cast<double>(t.tv_usec - start_time_.tv_usec));
  return time;
}

void Log::ResetTimer() {
  timer_ = 0;
  timer_is_on_ = false;
}

void Log::StartTimer() {
  timer_ = 0;
  gettimeofday(&timer_start_time_, NULL);
  timer_is_on_ = true;
}

void Log::ResumeTimer() {
  if (timer_is_on_)
    return;
  gettimeofday(&timer_start_time_, NULL);
  timer_is_on_ = true;
}

void Log::StopTimer() {
  if (!timer_is_on_)
    return;
  struct timeval t;
  gettimeofday(&t, NULL);
  timer_  += (static_cast<double>(t.tv_sec - timer_start_time_.tv_sec)) +
  .000001 * (static_cast<double>(t.tv_usec - timer_start_time_.tv_usec));
  timer_is_on_ = false;
}

void Log::ClearBuffer() {
  buffer_.str("");
}

void Log::ClearLogFile() {
  std::stringstream ss;
  ss.str("");
  std::ofstream ofs;
  ofs.open(log_file_name_.c_str());
  ofs << ss.str();
  ofs.close();
}

void Log::WriteToBuffer(std::string buf, LOG_TYPE type) {
  buffer_ << GetTag(type) << buf << std::endl;
}

void Log::WriteToFile(std::string buf, LOG_TYPE type) {
  std::ofstream ofs;
  ofs.open(log_file_name_.c_str(), std::ofstream::app);
  ofs << GetTag(type) << buf << std::endl;
  ofs.close();
}

void Log::WriteToOutputStream(std::string buf, LOG_TYPE type) {
  *output_stream_ << GetTag(type) << buf << std::endl;
}

void Log::WriteBufferToFile() {
  std::ofstream ofs;
  ofs.open(log_file_name_.c_str(), std::ofstream::app);
  ofs << buffer_.str();
  ofs.close();
}


void Log::WriteBufferToOutputStream() {
  *output_stream_ << buffer_.str();
}

void Log::PrintLine(std::string msg, LOG_TYPE type) {
  std::cout << GetTag(type) << msg << std::endl;
}

void Log::Print(std::string msg, LOG_TYPE type) {
  std::cout << GetTag(type) << msg;
}


std::string GetTag(LOG_TYPE type) {
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

