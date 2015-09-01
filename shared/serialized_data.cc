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
  * Class that represent the serialized data in Nimbus.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/serialized_data.h"

using namespace nimbus; // NOLINT

SerializedData::SerializedData()
  :size_(0) {
}

SerializedData::SerializedData(std::string str)
: data_ptr_(boost::shared_array<char>(new char[str.size()])), size_(str.size()) {
    memcpy(data_ptr_.get(), str.c_str(), size_);
}

SerializedData::SerializedData(const boost::shared_array<char>& data_ptr, const size_t& size)
: data_ptr_(data_ptr), size_(size) {
}

SerializedData::SerializedData(const SerializedData& other)
: data_ptr_(other.data_ptr_),
  size_(other.size_),
  header_(other.header_) {
}

SerializedData::~SerializedData() {
}

boost::shared_array<char> SerializedData::data_ptr() const {
  return data_ptr_;
}

char* SerializedData::data_ptr_raw() const {
  return data_ptr_.get();
}

std::string SerializedData::header() const {
  return header_;
}

void  SerializedData::set_data_ptr(boost::shared_array<char> ptr) {
  data_ptr_ = ptr;
}

void  SerializedData::set_data_ptr(char* ptr) {
  data_ptr_ = boost::shared_array<char> (ptr);
}

size_t SerializedData:: size() const {
  return size_;
}

void SerializedData:: set_size(size_t size) {
  size_ = size;
}

void SerializedData::set_header(const std::string str) const {
  header_ = str;
}

void SerializedData::set_header(char *str, size_t size) const {
  header_ = std::string(str, size);
}

bool SerializedData::Parse(const std::string& input) {
  std::string str = input;
  size_ = str.length();
  data_ptr_ = boost::shared_array<char>(new char[size_]);
  memcpy(data_ptr_.get(), str.c_str(), size_);
  return true;
}

std::string SerializedData::ToNetworkData() {
  if (size_ == 0) {
    std::string str = "empty";
    return str;
  } else {
    std::string str(data_ptr_.get(), size_);
    return str;
  }
}

SerializedData& SerializedData::operator= (const SerializedData& right) {
  data_ptr_ = right.data_ptr_;
  size_ = right.size_;
  return *this;
}

