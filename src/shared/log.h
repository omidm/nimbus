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
  * Nimbus log interface. Inoredr to exclude all the log command from the
  * binary file compile with -D_NIMBUS_NO_LOG flag.  
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SHARED_LOG_H_
#define NIMBUS_SRC_SHARED_LOG_H_

#include <boost/thread.hpp>
#include <sys/time.h>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream> // NOLINT
#include <string>

#define LOG_MAX_BUFF_SIZE  256000

#ifndef _NIMBUS_NO_LOG
#define log_Print(...) Print(__VA_ARGS__)
#define log_PrintLine(...) PrintLine(__VA_ARGS__)
#define log_WriteToFile(...) WriteToFile(__VA_ARGS__)
#define log_WriteToBuffer(...) WriteToBuffer(__VA_ARGS__)
#define log_WriteBufferToFile(...) WriteBufferToFile(__VA_ARGS__)
#define log_WriteToOutputStream(...) WriteToOutputStream(__VA_ARGS__)
#define log_WriteBufferToOutputStream(...) WriteBufferToOutputStream(__VA_ARGS__) // NOLINT
#define log_ResetTimer(...) ResetTimer(__VA_ARGS__)
#define log_StartTimer(...) StartTimer(__VA_ARGS__)
#define log_ResumeTimer(...) ResumeTimer(__VA_ARGS__)
#define log_StopTimer(...) StopTimer(__VA_ARGS__)
#define log_AddToTimer(...) AddToTimer(__VA_ARGS__)
#else
#define log_Print(...) none()
#define log_PrintLine(...) none()
#define log_WriteToFile(...) none()
#define log_WriteToBuffer(...) none()
#define log_WriteBufferToFile(...) none()
#define log_WriteToOutputStream(...) none()
#define log_WriteBufferToOutputStream(...) none()
#define log_ResetTimer(...) none()
#define log_StartTimer(...) none()
#define log_ResumeTimer(...) none()
#define log_StopTimer(...) none()
#define log_AddToTimer(...) none()
#endif




enum LOG_TYPE {
  LOG_ERROR,
  LOG_WARNING,
  LOG_INFO,
  LOG_DEBUG,
  LOG_NONE
};

std::string GetTag(LOG_TYPE type);

class Log {
  public:
    enum Type {
      NO_FILE = 0
    };

    Log();
    explicit Log(Type);
    explicit Log(std::ostream* os);
    explicit Log(std::string fname);
    Log(std::ostream* os, std::string fname);
    ~Log();

    void set_output_stream(std::ostream* os);
    void set_file_name(std::string fname);
    void set_start_time(struct timeval* time);

    struct timeval* start_time();
    double timer();

    void InitTime();
    double GetTime();

    void ResetTimer();
    void StartTimer();
    void ResumeTimer();
    void StopTimer();
    void AddToTimer(double elapsed);

    void ClearBuffer();
    void ClearLogFile();

    void WriteToBuffer(std::string buf, LOG_TYPE type = LOG_NONE);
    void WriteToFile(std::string buf, LOG_TYPE type = LOG_NONE);
    void WriteToOutputStream(std::string buf, LOG_TYPE type = LOG_NONE);

    void WriteBufferToFile();
    void WriteBufferToOutputStream();

    static void PrintLine(std::string msg, LOG_TYPE type = LOG_NONE);
    static void Print(std::string msg, LOG_TYPE type = LOG_NONE);

    static double GetRawTime();

    static void none() {}

  private:
    std::ostream* output_stream_;
    std::string log_file_name_;
    std::stringstream buffer_;
    struct timeval start_time_;
    struct timeval timer_start_time_;
    double timer_;
    bool timer_is_on_;

    boost::mutex file_mutex_;
    boost::mutex stream_mutex_;
    boost::mutex buffer_mutex_;
};


#endif  // NIMBUS_SRC_SHARED_LOG_H_
