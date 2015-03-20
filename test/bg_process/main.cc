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

#define DEFAULT_THREAD_NUM 1
#define DEFAULT_ARRAY_SIZE 1  // In (K)

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <iostream> // NOLINT
#include <list>
#include <vector>

void Thread(size_t size) {
  std::list<size_t> array;
  for (size_t i = 0; i < size; ++i) {
    array.push_back(i);
  }
  while (true) {
    size_t add = 0;
    std::list<size_t>::iterator iter;
    for (iter = array.begin(); iter != array.end(); ++iter) {
      add += *iter;
    }
    for (iter = array.begin(); iter != array.end(); ++iter) {
      *iter = *iter + add;
    }
  }
}

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t array_size;
  size_t thread_num;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce help message")
    // Optinal arguments
    ("thread_num,t", po::value<size_t>(&thread_num)->default_value(DEFAULT_THREAD_NUM), "number of threads") // NOLINT
    ("array_size,s", po::value<size_t>(&array_size)->default_value(DEFAULT_ARRAY_SIZE), "size of array (K)"); // NOLINT

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  try {
    po::notify(vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return -1;
  }


  std::list<boost::thread*> threads_;

  for (size_t i = 0; i < thread_num; ++i) {
    boost::thread *t = new boost::thread(boost::bind(Thread, array_size * 1000));
    threads_.push_back(t);
  }

  std::list<boost::thread*>::iterator iter = threads_.begin();
  for (; iter != threads_.end(); ++iter) {
    (*iter)->join();
  }

  return 0;
}


