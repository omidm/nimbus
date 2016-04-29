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
 * Data classes for logestic regression application.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <stdlib.h>
#include "./data.h"

Sample::Sample(const size_t& dimension)
  : dimension_(dimension),
    vector_(dimension) {
}

Sample::~Sample() {
}

size_t Sample::dimension() const {
  return dimension_;
}

std::vector<double>* Sample::vector() {
  return &vector_;
}


SampleBatch::SampleBatch(const size_t& dimension,
                         const size_t& sample_num)
  : dimension_(dimension),
    sample_num_(sample_num) {
  set_name(SAMPLE_BATCH_DATA_NAME);
}

SampleBatch::~SampleBatch() {
}

void SampleBatch::Create() {
  for (size_t i = 0; i < sample_num_; ++i) {
    samples_.push_back(Sample(dimension_));
  }
}

void SampleBatch::Destroy() {
  samples_.clear();
}

Data* SampleBatch::Clone() {
  return new SampleBatch(dimension_, sample_num_);
}

void SampleBatch::Copy(Data* from) {
  SampleBatch *data = static_cast<SampleBatch*>(from);
  assert(data->name() == SAMPLE_BATCH_DATA_NAME);
  assert(sample_num_ == data->sample_num());
  samples_.clear();
  typename std::vector<Sample>::const_iterator iter = data->samples()->begin();
  for (; iter != data->samples()->end(); ++iter) {
    samples_.push_back(*iter);
  }
}

bool SampleBatch::Serialize(SerializedData* ser_data) {
  // currently we only support double!
  data_msgs::SampleBatchMsg sample_batch_msg;
  typename std::vector<Sample>::iterator iter = samples_.begin();
  for (; iter != samples_.end(); ++iter) {
    data_msgs::SampleMsg *sample_msg = sample_batch_msg.add_samples();
    for (size_t i = 0; i < iter->dimension(); i++) {
      sample_msg->add_elems(iter->vector()->operator[](i));
    }
  }

  std::string str;
  sample_batch_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());

  return true;
}


bool SampleBatch::DeSerialize(const SerializedData& ser_data, Data** result) {
  // currently we only support double!
  data_msgs::SampleBatchMsg sample_batch_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  sample_batch_msg.ParseFromString(str);
  SampleBatch* sb = new SampleBatch(dimension_, sample_num_);
  sb->Create();
  assert(size_t(sample_batch_msg.samples_size()) == sample_num_);
  for (size_t i = 0; i < sample_num_; i++) {
    data_msgs::SampleMsg sample_msg = sample_batch_msg.samples(i);
    assert(size_t(sample_msg.elems_size()) == dimension_);
    for (size_t j =0; j < dimension_; ++j) {
      sb->samples_[i].vector()->operator[](j) = sample_msg.elems(j);
    }
  }

  *result = sb;
  return true;
}

size_t SampleBatch::dimension() const {
  return dimension_;
}

size_t SampleBatch::sample_num() const {
  return sample_num_;
}

std::vector<Sample>* SampleBatch::samples() {
  return &samples_;
}


Mean::Mean(const size_t& dimension)
  : dimension_(dimension),
    vector_(dimension),
    scratch_(dimension) {
  scratch_weight_ = 0;
}

Mean::~Mean() {
}

size_t Mean::dimension() const {
  return dimension_;
}

std::vector<double>* Mean::vector() {
  return &vector_;
}

std::vector<double>* Mean::scratch() {
  return &scratch_;
}

size_t Mean::scratch_weight() const {
  return scratch_weight_;
}

void Mean::set_scratch_weight(const size_t& w) {
  scratch_weight_ = w;
}


Means::Means(const size_t& dimension,
             const size_t& cluster_num,
             const std::string& name)
  : dimension_(dimension),
    cluster_num_(cluster_num) {
  set_name(name);
}

Means::~Means() {
}

void Means::Create() {
  for (size_t i = 0; i < cluster_num_; ++i) {
    means_.push_back(Mean(dimension_));
  }
}

void Means::Destroy() {
  means_.clear();
}

Data* Means::Clone() {
  return new Means(dimension_, cluster_num_, name());
}

void Means::Copy(Data* from) {
  Means *data = static_cast<Means*>(from); // NOLINT
  assert(data->name() == MEANS_DATA_NAME);
  means_.clear();
  typename std::vector<Mean>::const_iterator iter = data->means()->begin();
  for (; iter != data->means()->end(); ++iter) {
    means_.push_back(*iter);
  }
}

bool Means::Serialize(SerializedData* ser_data) {
  // currently we only support double!
  data_msgs::MeansMsg means_msg;
  typename std::vector<Mean>::iterator iter = means_.begin();
  for (; iter != means_.end(); ++iter) {
    data_msgs::MeanMsg *mean_msg = means_msg.add_means();
    mean_msg->set_scratch_weight(iter->scratch_weight());
    for (size_t i = 0; i < iter->dimension(); i++) {
      mean_msg->add_vector_elems(iter->vector()->operator[](i));
      mean_msg->add_scratch_elems(iter->scratch()->operator[](i));
    }
  }
  std::string str;
  means_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool Means::DeSerialize(const SerializedData& ser_data, Data** result) {
  // currently we only support double!
  data_msgs::MeansMsg means_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  means_msg.ParseFromString(str);
  Means* m = new Means(dimension_, cluster_num_, name());
  m->Create();
  assert(size_t(means_msg.means_size()) == cluster_num_);
  for (size_t i = 0; i < cluster_num_; i++) {
    data_msgs::MeanMsg mean_msg = means_msg.means(i);
    assert(size_t(mean_msg.vector_elems_size()) == dimension_);
    assert(size_t(mean_msg.scratch_elems_size()) == dimension_);
    m->means_[i].set_scratch_weight(mean_msg.scratch_weight());
    for (size_t j =0; j < dimension_; ++j) {
      m->means_[i].vector()->operator[](j) = mean_msg.vector_elems(j);
      m->means_[i].scratch()->operator[](j) = mean_msg.scratch_elems(j);
    }
  }

  *result = m;
  return true;
}

size_t Means::dimension() const {
  return dimension_;
}

size_t Means::cluster_num() const {
  return cluster_num_;
}

std::vector<Mean>* Means::means() {
  return &means_;
}

// void Means::set_vector(const std::vector<double>& vector) {
//   vector_ = vector;
// }
//
// void Means::set_gradient(const std::vector<double>& gradient) {
//   gradient_ = gradient;
// }

