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

float Sample::label() const {
  return label_;
}

std::vector<float>* Sample::vector() {
  return &vector_;
}

void Sample::set_label(float label) {
  label_ = label;
}

SampleBatch::SampleBatch(const size_t& dimension,
                         const size_t& sample_num)
  : dimension_(dimension),
    sample_num_(sample_num) {
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
  // currently we only support float!
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
  // currently we only support float!

  data_msgs::SampleBatchMsg sample_batch_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  sample_batch_msg.ParseFromString(str);
  SampleBatch* sb = new SampleBatch(dimension_, sample_num_);
  sb->Create();
  assert(size_t(sample_batch_msg.samples_size()) == sample_num_);
  for (size_t i = 0; i < dimension_; i++) {
    data_msgs::SampleMsg sample_msg = sample_batch_msg.samples(i);
    Sample sample(dimension_);
    sample.set_label(sample_msg.label());
    assert(size_t(sample_msg.elems_size()) == dimension_);
    for (size_t j =0; j < dimension_; ++j) {
      sample.vector()->operator[](j) = sample_msg.elems(j);
    }
    sb->samples_[i] = sample;
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


Weight::Weight(const size_t& dimension)
  : dimension_(dimension) {
}

Weight::~Weight() {
}

void Weight::Create() {
  vector_.resize(dimension_);
}

void Weight::Destroy() {
  vector_.clear();
}

Data* Weight::Clone() {
  return new Weight(dimension_);
}

void Weight::Copy(Data* from) {
  Weight *data = reinterpret_cast<Weight*>(from); // NOLINT
  assert(data);
  vector_.clear();
  typename std::vector<float>::const_iterator iter = data->vector()->begin();
  for (; iter != data->vector()->end(); ++iter) {
    vector_.push_back(*iter);
  }
}

bool Weight::Serialize(SerializedData* ser_data) {
  // currently we only support float!
  data_msgs::VectorMsg vec_msg;
  for (size_t i = 0; i < dimension_; i++) {
    vec_msg.add_elems(vector_[i]);
  }
  std::string str;
  vec_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool Weight::DeSerialize(const SerializedData& ser_data, Data** result) {
  // currently we only support float!
  data_msgs::VectorMsg vec_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  vec_msg.ParseFromString(str);
  Weight* w = new Weight(dimension_);
  w->Create();
  assert(size_t(vec_msg.elems_size()) == dimension_);
  for (size_t i = 0; i < dimension_; i++) {
     w->vector_[i] = vec_msg.elems(i);
  }

  *result = w;
  return true;
}

size_t Weight::dimension() const {
  return dimension_;
}

const std::vector<float>* Weight::vector() const {
  return &vector_;
}

void Weight::set_vector(const std::vector<float>& vector) {
  vector_ = vector;
}







