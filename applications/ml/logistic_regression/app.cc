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
 * Distributed Logistic Regression application
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <boost/program_options.hpp>
#include "./app.h"
#include "./job.h"
#include "./data.h"

#define DEFAULT_DIMENSION 10
#define DEFAULT_ITERATION_NUM 10
#define DEFAULT_PARTITION_NUM 10
#define DEFAULT_SAMPLE_NUM_M 1
#define DEFAULT_REDUCTION_PARTITION_NUM 0

LogisticRegression::LogisticRegression(const size_t& dimension,
                                       const size_t& iteration_num,
                                       const size_t& partition_num,
                                       const size_t& sample_num_m,
                                       const size_t& reduction_partition_num)
  : dimension_(dimension),
    iteration_num_(iteration_num),
    partition_num_(partition_num),
    sample_num_m_(sample_num_m),
    reduction_partition_num_(reduction_partition_num) {
      assert(sizeof(size_t) == 8); // NOLINT
      assert(sizeof(double) == 8); // NOLINT
      assert(((sample_num_m * size_t(1e6)) % partition_num) == 0);
      if (reduction_partition_num != 0) {
        assert((partition_num % reduction_partition_num) == 0);
      }
      sample_num_per_partition_ = (sample_num_m * 1e6) / partition_num;
      std::cout << "**** number of partitions:        " << partition_num_ << std::endl;
      std::cout << "**** sample number per partition: " << sample_num_per_partition_ << std::endl;
      std::cout << "**** number of reduction partitions:        " << reduction_partition_num_ << std::endl; // NOLINT
      reduction_combiner_active_ = true;
}

LogisticRegression::~LogisticRegression() {
};


size_t LogisticRegression::iteration_num() {
  return iteration_num_;
}

size_t LogisticRegression::dimension() {
  return dimension_;
}

size_t LogisticRegression::partition_num() {
  return partition_num_;
}

size_t LogisticRegression::sample_num_m() {
  return sample_num_m_;
}

size_t LogisticRegression::reduction_partition_num() {
  return reduction_partition_num_;
}

size_t LogisticRegression::sample_num_per_partition() {
  return sample_num_per_partition_;
}

bool LogisticRegression::reduction_combiner_active() {
  return reduction_combiner_active_;
}

void LogisticRegression::set_reduction_combiner_active(bool flag) {
  reduction_combiner_active_ = flag;
}

void LogisticRegression::Load() {
  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(INIT_JOB_NAME, new Init(this));
  RegisterJob(LOOP_JOB_NAME, new ForLoop(this));
  RegisterJob(GRADIENT_JOB_NAME, new Gradient(this));
  RegisterJob(REDUCE_JOB_NAME, new Reduce(this));
  RegisterJob(COMBINE_JOB_NAME, new Combine(this));
  RegisterJob(REDUCE_L1_JOB_NAME, new ReduceL1(this));
  RegisterJob(REDUCE_L2_JOB_NAME, new ReduceL2(this));
  RegisterJob(REDUCE_L3_JOB_NAME, new ReduceL3(this));

  RegisterData(WEIGHT_DATA_NAME, new Weight(dimension_, WEIGHT_DATA_NAME));
  RegisterData(SCRATCH_WEIGHT_DATA_NAME, new Weight(dimension_, SCRATCH_WEIGHT_DATA_NAME));
  RegisterData(SAMPLE_BATCH_DATA_NAME, new SampleBatch(dimension_, sample_num_per_partition_));
};

extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t dimension;
  size_t iteration_num;
  size_t partition_num;
  size_t sample_num_m;
  size_t reduction_partition_num;

  po::options_description desc("Logistic Regression Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("dimension,d", po::value<std::size_t>(&dimension)->default_value(DEFAULT_DIMENSION), "dimension of the sample vectors") // NOLINT
    ("iteration,i", po::value<std::size_t>(&iteration_num)->default_value(DEFAULT_ITERATION_NUM), "number of iterations") // NOLINT
    ("pn", po::value<std::size_t>(&partition_num)->default_value(DEFAULT_PARTITION_NUM), "number of partitions") // NOLINT
    ("sn", po::value<std::size_t>(&sample_num_m)->default_value(DEFAULT_SAMPLE_NUM_M), "number of samples in Million") // NOLINT
    ("rpn", po::value<std::size_t>(&reduction_partition_num)->default_value(DEFAULT_REDUCTION_PARTITION_NUM), "number of reduction partitions for manual reduction by application with read/write set, if not specified nimbus handles reduction automatically through scratch/reduce set.") // NOLINT
    ("drc", "deactivate reduction combiner"); // NOLINT

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return NULL;
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return NULL;
  }

  try {
    po::notify(vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return NULL;
  }

  LogisticRegression *app = new LogisticRegression(dimension,
                                                   iteration_num,
                                                   partition_num,
                                                   sample_num_m,
                                                   reduction_partition_num);

  if (vm.count("drc")) {
    app->set_reduction_combiner_active(false);
  }

  return app;
}

