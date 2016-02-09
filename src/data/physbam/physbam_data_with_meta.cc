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
  * Author: Hang Qu <quhang@stanford.edu>
  */
#include <boost/functional/hash.hpp>
#include <string>
#include "src/data/physbam/physbam_data.h"
#include "src/data/physbam/protobuf_compiled/pd_message.pb.h"
#include "src/shared/dbg.h"
#include "src/data/physbam/physbam_data_with_meta.h"

namespace nimbus {

PhysBAMDataWithMeta::PhysBAMDataWithMeta() {
  ResetMetaData();
}

Data* PhysBAMDataWithMeta::Clone() {
  PhysBAMDataWithMeta* d = new PhysBAMDataWithMeta();
  d->set_buffer(NULL, size_);
  return d;
}

void PhysBAMDataWithMeta::Create() {
  PhysBAMData::Create();
  ResetMetaData();
}

void PhysBAMDataWithMeta::Destroy() {
  PhysBAMData::Destroy();
  ResetMetaData();
}

void PhysBAMDataWithMeta::Copy(Data *from) {
  PhysBAMData::Copy(from);
  PhysBAMDataWithMeta* pfrom = dynamic_cast<PhysBAMDataWithMeta*>(from);  // NOLINT
  assert(pfrom != NULL);
  has_meta_data_ = pfrom->has_meta_data();
  meta_data_size_ = pfrom->meta_data_size();
  meta_data_hash_ = pfrom->meta_data_hash();
}

bool PhysBAMDataWithMeta::Serialize(SerializedData* ser_data) {
  nimbus_message::pd_message pd; // NOLINT
  if (buffer_) {
    std::string buf(buffer_, size_);
    pd.set_buffer(buf);
  }
  if (size_) {
    pd.set_size(size_);
  } else {
    pd.set_size(0);
  }
  pd.set_hash(0);
  if (has_meta_data_) {
    pd.set_meta_data_size(meta_data_size_);
    pd.set_meta_data_hash(meta_data_hash_);
  }
  std::string ser;
  bool success = pd.SerializeToString(&ser);
  char *buf = new char[ser.length() + 1];
  memcpy(buf, ser.c_str(), sizeof(char) * (ser.length() + 1)); // NOLINT
  if (!success)
    return success;
  ser_data->set_data_ptr(buf);
  ser_data->set_size(sizeof(char) * ser.length() + 1); // NOLINT
  return success;
}

bool PhysBAMDataWithMeta::DeSerialize(const SerializedData& ser_data,
                                      Data** result) {
  const char *buf = ser_data.data_ptr_raw();
  const int buf_size = ser_data.size();
  if (buf_size <= 0)
    return false;
  assert(buf);

  nimbus_message::pd_message pd; // NOLINT
  std::string temp(buf, buf_size-1); // NOLINT
  bool success = pd.ParseFromString(temp);
  if (!success)
    return success;

  PhysBAMDataWithMeta *data = new PhysBAMDataWithMeta();
  (*result) = static_cast<Data *>(data);
  if (pd.has_buffer()) {
    int size = pd.size();
    char *buffer = new char[size];
    memcpy(buffer, pd.buffer().c_str(), sizeof(char) * size); // NOLINT
    data->set_buffer(buffer, size);
  } else {
    data->set_buffer(NULL, pd.size());
  }

  data->ResetMetaData();
  if (!pd.has_meta_data_size()) {
    dbg(DBG_WARN, "No meta data received.\n");
  }
  if (pd.has_meta_data_size()) {
    data->set_has_meta_data();
    data->set_meta_data_size(pd.meta_data_size());
  }
  if (pd.has_meta_data_hash()) {
    data->set_meta_data_hash(pd.meta_data_hash());
  }
  return success;
}

void PhysBAMDataWithMeta::ResetMetaData() {
  has_meta_data_ = false;
  meta_data_size_ = 0;
  meta_data_hash_ = 0;
}

void PhysBAMDataWithMeta::MarkMetaDataInTempBuffer() {
  meta_data_size_ = temp_buffer_->tellp();
  has_meta_data_ = true;
  std::size_t temp = HASH_SEED;
  /*
  const std::string& temp_str = temp_buffer_->str();
  const char* pointer = temp_str.c_str();
  if (meta_data_size_ != 0) {
    boost::hash_range(temp, pointer, pointer + meta_data_size_);
  }
  */
  meta_data_hash_ = static_cast<int64_t>(temp);
}

}  // namespace nimbus
