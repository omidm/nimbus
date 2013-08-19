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
  * This is the code for parsing the log.txt files dumped by Nimbus
  * scheduler/worker.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include <boost/tokenizer.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream> // NOLINT
#include <string>
#include <map>

using boost::tokenizer;
using boost::char_separator;


int main(int argc, char *argv[]) {
  std::ifstream input;
  std::string file_name;
  double useful_time = 0;
  double start_time, end_time;
  std::map<std::string, double> time_profile;

  if (argc <= 1)
    file_name = "log.txt";
  else
    file_name = argv[1];
  std::cout << "file name: " << file_name << std::endl;

  input.open(file_name.c_str());
  if (input.fail()) {
    std::cout << "Could not open file: " << file_name << std::endl;
    return -1;
  }
  while (true) {
    std::string line;
    getline(input, line);
    if (input.fail()) break;

    // std::cout << line << std::endl;
    char_separator<char> separator(" \t");
    tokenizer<char_separator<char> > tokens(line, separator);
    tokenizer<char_separator<char> >::iterator iter = tokens.begin();

    std::string job_name;
    iter++; iter++; iter++; iter++; job_name = *iter;
    if (time_profile.find(job_name) == time_profile.end())
      time_profile[job_name] = 0;

    double time;
    iter++; iter++; iter++; iter++;
    std::stringstream ss(*iter);
    ss >> time;
    time_profile[job_name] = time_profile[job_name] + time;
    // std::cout << job_name << time << std::endl;

    iter++; iter++;
    std::stringstream sss(*iter);
    sss >> end_time;
    if (job_name == "main")
      start_time = end_time;
  }
  input.close();

  std::map<std::string, double>::iterator iter;
  for (iter = time_profile.begin(); iter != time_profile.end(); iter++) {
    printf("name: %25s time(s): %6.3lf\n",
        iter->first.c_str(), iter->second/1000);
    useful_time += iter->second/1000;
  }
  printf("\n%-20s %6.3lf\n", "total time(s):", end_time - start_time);
  printf("%-20s %6.3lf\n", "useful time(s):", useful_time);
  printf("%-20s %6.3lf\n", "nimbus overhead(s):", end_time - start_time - useful_time);
  printf("%-20s %6.3lf\n", "efficiency:", useful_time/(end_time - start_time));
}






