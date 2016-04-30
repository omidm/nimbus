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
 * Distributed K-Means application
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <boost/program_options.hpp>
#include "./app.h"
#include "./job.h"
#include "./data.h"

#define DEFAULT_DIMENSION 10
#define DEFAULT_CLUSTER_NUM 2
#define DEFAULT_ITERATION_NUM 10
#define DEFAULT_PARTITION_NUM 10
#define DEFAULT_SAMPLE_NUM_M 1
#define DEFAULT_SPIN_WAIT_US 0
#define DEFAULT_REDUCTION_PARTITION_NUM 1

KMeans::KMeans(const size_t& dimension,
               const size_t& cluster_num,
               const size_t& iteration_num,
               const size_t& partition_num,
               const double& sample_num_m,
               const size_t& spin_wait_us,
               const size_t& reduction_partition_num)
  : dimension_(dimension),
    cluster_num_(cluster_num),
    iteration_num_(iteration_num),
    partition_num_(partition_num),
    sample_num_m_(sample_num_m),
    spin_wait_us_(spin_wait_us),
    reduction_partition_num_(reduction_partition_num) {
      assert(sizeof(size_t) == 8); // NOLINT
      assert(sizeof(double) == 8); // NOLINT
      assert((size_t(sample_num_m * size_t(1e6)) % partition_num) == 0);
      assert((partition_num % reduction_partition_num) == 0);
      sample_num_per_partition_ = (sample_num_m * 1e6) / partition_num;
      dbg(DBG_APP, "APPLICATION: number of partitions:        %lu\n", partition_num_);
      dbg(DBG_APP, "APPLICATION: sample number per partition: %lu\n", sample_num_per_partition_);
      automatic_reduction_active_ = true;
      reduction_combiner_active_ = true;
}

KMeans::~KMeans() {
};

size_t KMeans::iteration_num() {
  return iteration_num_;
}

size_t KMeans::dimension() {
  return dimension_;
}

size_t KMeans::cluster_num() {
  return cluster_num_;
}

size_t KMeans::partition_num() {
  return partition_num_;
}

double KMeans::sample_num_m() {
  return sample_num_m_;
}

size_t KMeans::spin_wait_us() {
  return spin_wait_us_;
}

size_t KMeans::reduction_partition_num() {
  return reduction_partition_num_;
}

size_t KMeans::sample_num_per_partition() {
  return sample_num_per_partition_;
}

bool KMeans::automatic_reduction_active() {
  return automatic_reduction_active_;
}

bool KMeans::reduction_combiner_active() {
  return reduction_combiner_active_;
}

void KMeans::set_automatic_reduction_active(bool flag) {
  automatic_reduction_active_ = flag;
  if (!automatic_reduction_active_) {
    dbg(DBG_APP, "APPLICATION: automatic reduction deactivated!");
  }
}

void KMeans::set_reduction_combiner_active(bool flag) {
  reduction_combiner_active_ = flag;
  if (!reduction_combiner_active_) {
    dbg(DBG_APP, "APPLICATION: reduction combiner deactivated!");
  }
}

void KMeans::Load() {
  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(INIT_SAMPLES_JOB_NAME, new InitSamples(this));
  RegisterJob(INIT_MEANS_JOB_NAME, new InitMeans(this));
  RegisterJob(LOOP_JOB_NAME, new ForLoop(this));
  RegisterJob(CLUSTER_JOB_NAME, new Cluster(this));

  RegisterJob(REDUCE_JOB_NAME, new Reduce(this));
  RegisterJob(COMBINE_JOB_NAME, new Combine(this));

  RegisterJob(REDUCE_L1_JOB_NAME, new ReduceL1(this));
  RegisterJob(REDUCE_L2_JOB_NAME, new ReduceL2(this));
  RegisterJob(SYNCH_JOB_NAME, new Synch(this));

  RegisterData(MEANS_DATA_NAME, new Means(dimension_, cluster_num_, MEANS_DATA_NAME));
  RegisterData(SAMPLE_BATCH_DATA_NAME, new SampleBatch(dimension_, sample_num_per_partition_));

  RegisterData(SCRATCH_MEANS_DATA_NAME, new Means(dimension_, cluster_num_, SCRATCH_MEANS_DATA_NAME)); // NOLINT
};

extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t dimension;
  size_t cluster_num;
  size_t iteration_num;
  size_t partition_num;
  double sample_num_m;
  size_t spin_wait_us;
  size_t reduction_partition_num;


  po::options_description desc("K-Means Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("dimension,d", po::value<std::size_t>(&dimension)->default_value(DEFAULT_DIMENSION), "dimension of the sample vectors") // NOLINT
    ("iteration,i", po::value<std::size_t>(&iteration_num)->default_value(DEFAULT_ITERATION_NUM), "number of iterations") // NOLINT
    ("sample_num_m,s", po::value<double>(&sample_num_m)->default_value(DEFAULT_SAMPLE_NUM_M), "number of samples in Million") // NOLINT
    ("cluster_num,c", po::value<std::size_t>(&cluster_num)->default_value(DEFAULT_CLUSTER_NUM), "number of clusters") // NOLINT
    ("partition_num,p", po::value<std::size_t>(&partition_num)->default_value(DEFAULT_PARTITION_NUM), "number of partitions") // NOLINT
    ("spin_wait,w", po::value<std::size_t>(&spin_wait_us)->default_value(DEFAULT_SPIN_WAIT_US), "spin wait in micro seconds, if non zero,replaces the gradient operation with fixed spin wait.") // NOLINT
    ("reduction_partition_num,r", po::value<std::size_t>(&reduction_partition_num)->default_value(DEFAULT_REDUCTION_PARTITION_NUM), "number of reduction partitions for manual reduction by application with read/write set.") // NOLINT
    ("dar", "deactivate automatic reduction") // NOLINT
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

  KMeans *app = new KMeans(dimension,
                           cluster_num,
                           iteration_num,
                           partition_num,
                           sample_num_m,
                           spin_wait_us,
                           reduction_partition_num);

  if (vm.count("dar")) {
    app->set_automatic_reduction_active(false);
  }

  if (vm.count("drc")) {
    app->set_reduction_combiner_active(false);
  }

  return app;
}


