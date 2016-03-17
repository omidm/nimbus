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
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <boost/program_options.hpp>
#include <vector>
#include "src/shared/nimbus.h"
#include "src/shared/nimbus_types.h"

#define DEFAULT_DIMENSION 10
#define DEFAULT_ITERATION_NUM 10
#define DEFAULT_PARTITION_NUM 10
#define DEFAULT_SAMPLE_NUM_M 1

using namespace nimbus; // NOLINT

class Sample {
  public:
    Sample() {}
    ~Sample() {}
    double label_;
    std::vector<double> vector_;
    std::vector<double>* vector() { return &vector_;}
};

double VectorDotProduct(const std::vector<double>* vec1,
                        const std::vector<double>* vec2) {
  assert(vec1->size() == vec2->size());
  double result = 0;
  for (size_t i = 0; i < vec1->size(); ++i) {
    result += vec1->operator[](i) * vec2->operator[](i);
  }
  return result;
}

void VectorAddWithScale(std::vector<double>* acc,
                        const std::vector<double>* add,
                        const double& scale) {
  assert(acc->size() == add->size());
  for (size_t i = 0; i < acc->size(); ++i) {
    acc->operator[](i) += add->operator[](i) * scale;
  }
}



int main(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t dimension;
  size_t iteration_num;
  double sample_num_m;

  po::options_description desc("Gradient Program Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("dimension,d", po::value<std::size_t>(&dimension)->default_value(DEFAULT_DIMENSION), "dimension of the sample vectors") // NOLINT
    ("iteration,i", po::value<std::size_t>(&iteration_num)->default_value(DEFAULT_ITERATION_NUM), "number of iterations") // NOLINT
    ("sn", po::value<double>(&sample_num_m)->default_value(DEFAULT_SAMPLE_NUM_M), "number of samples in Million"); // NOLINT

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  } catch (std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return -1;
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  try {
    po::notify(vm);
  } catch (std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return -1;
  }

  size_t sample_num = (size_t)(sample_num_m * 1e6);

  printf("** dimension    : %lu\n", dimension);
  printf("** iteration_num: %lu\n", iteration_num);
  printf("** sample_num   : %lu\n", sample_num);

  Log log(Log::NO_FILE);
  std::vector<Sample> samples;
  log.log_StartTimer();
  for (size_t i = 0; i < sample_num; ++i) {
    samples.push_back(Sample());
    Sample& s = samples.back();
    // s.label_ = i % 2;
    s.label_ = static_cast<double>((i % 2) *2) - 1;  // for +1 and -1 labels
    for (size_t j = 0; j < dimension; ++j) {
      s.vector_.push_back(13);
    }
  }
  log.log_StopTimer();
  printf("Sample generation took: %2.2f(ms).\n", log.timer() * 1000);

  std::vector<double> weight(dimension, 1);

  double sum_elapsed = 0;
  for (size_t i = 0; i < iteration_num; ++i) {
    log.log_StartTimer();
    std::vector<double> gradient(dimension, 0);
    std::vector<Sample>::iterator iter = samples.begin();
    for (; iter != samples.end(); ++iter) {
      double l = iter->label_;
      double scale = (1 / (1 + exp(l * VectorDotProduct(iter->vector(), &weight))) - 1) * l;
      VectorAddWithScale(&gradient, iter->vector(), scale);
    }
    VectorAddWithScale(&weight, &gradient, 1);
    log.log_StopTimer();
    printf("Iteration %lu elapsed time: %2.2f(ms).\n", i, log.timer() * 1000);
    sum_elapsed += log.timer() * 1000;
  }

  printf("Average gradient for %lu iterations %2.2f(ms).\n\n",
      iteration_num, sum_elapsed/iteration_num);
}

