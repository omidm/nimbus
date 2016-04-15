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

double Sample::label() const {
  return label_;
}

std::vector<double>* Sample::vector() {
  return &vector_;
}

void Sample::set_label(double label) {
  label_ = label;
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
  SampleBatch *data = reinterpret_cast<SampleBatch*>(from); // NOLINT
  assert(data);
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
    sample_msg->set_label(iter->label());
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
    sb->samples_[i].set_label(sample_msg.label());
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


Weight::Weight(const size_t& dimension,
               const std::string& name)
  : dimension_(dimension) {
  set_name(name);
}

Weight::~Weight() {
}

void Weight::Create() {
  vector_.resize(dimension_, 1);
  gradient_.resize(dimension_, 0);
}

void Weight::Destroy() {
  vector_.clear();
  gradient_.clear();
}

Data* Weight::Clone() {
  return new Weight(dimension_, name());
}

void Weight::Copy(Data* from) {
  Weight *data = reinterpret_cast<Weight*>(from); // NOLINT
  assert(data);
  vector_.clear();
  {
    typename std::vector<double>::const_iterator iter = data->vector()->begin();
    for (; iter != data->vector()->end(); ++iter) {
      vector_.push_back(*iter);
    }
  }
  gradient_.clear();
  {
    typename std::vector<double>::const_iterator iter = data->gradient()->begin();
    for (; iter != data->gradient()->end(); ++iter) {
      gradient_.push_back(*iter);
    }
  }
}

bool Weight::Serialize(SerializedData* ser_data) {
  // currently we only support double!
  data_msgs::WeightMsg weight_msg;
  for (size_t i = 0; i < dimension_; i++) {
    weight_msg.add_vector_elems(vector_[i]);
  }
  for (size_t i = 0; i < dimension_; i++) {
    weight_msg.add_gradient_elems(gradient_[i]);
  }
  std::string str;
  weight_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool Weight::DeSerialize(const SerializedData& ser_data, Data** result) {
  // currently we only support double!
  data_msgs::WeightMsg weight_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  weight_msg.ParseFromString(str);
  Weight* w = new Weight(dimension_, name());
  w->Create();
  assert(size_t(weight_msg.vector_elems_size()) == dimension_);
  for (size_t i = 0; i < dimension_; i++) {
     w->vector_[i] = weight_msg.vector_elems(i);
  }
  assert(size_t(weight_msg.gradient_elems_size()) == dimension_);
  for (size_t i = 0; i < dimension_; i++) {
     w->gradient_[i] = weight_msg.gradient_elems(i);
  }

  *result = w;
  return true;
}

size_t Weight::dimension() const {
  return dimension_;
}

std::vector<double>* Weight::vector() {
  return &vector_;
}

std::vector<double>* Weight::gradient() {
  return &gradient_;
}

void Weight::set_vector(const std::vector<double>& vector) {
  vector_ = vector;
}

void Weight::set_gradient(const std::vector<double>& gradient) {
  gradient_ = gradient;
}

