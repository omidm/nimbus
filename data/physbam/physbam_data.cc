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

/***********************************************************************
 * AUTHOR: Philip Levis <pal>
 *   FILE: .//physbam_data.cc
 *   DATE: Thu Dec 12 17:15:37 2013
 *  DESCR:
 ***********************************************************************/
#include <string>
#include "data/physbam/physbam_data.h"
#include "data/physbam/protobuf_compiled/pd_message.pb.h"
#include "shared/dbg.h"

namespace nimbus {
/**
 * \fn nimbus::PhysBAMData::PhysBAMData()
 * \brief Brief description.
 * \return
*/
PhysBAMData::PhysBAMData(): size_(0), buffer_(0), temp_buffer_(0) {}

uint32_t PhysBAMData::HashCode() {
  if (buffer_ == NULL) {
    return 0;
  }
  uint32_t hash = 0;
  for (uint32_t i = 0; i < size_; ++i) {
    hash += buffer_[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }
  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);
  return hash;
}

char* PhysBAMData::buffer() const {
    return buffer_;
}

void PhysBAMData::set_buffer(char *b, int_dimension_t s) {
    buffer_ = b;
    size_ = s;
}

int_dimension_t PhysBAMData::size() const {
    return size_;
}

void PhysBAMData::set_size(int_dimension_t s) {
    size_ = s;
}

/**
 * \fn Data * nimbus::PhysBAMData::Clone()
 * \brief Brief description.
 * \return
*/
Data * PhysBAMData::Clone() {
  PhysBAMData* d = new PhysBAMData();
  d->set_buffer(NULL, size_);
  return d;
}


/**
 * \fn void nimbus::PhysBAMData::Create()
 * \brief Brief description.
 * \return
*/
void PhysBAMData::Create() {
  if (size_ && !buffer_) {
      buffer_ = static_cast<char*>(malloc(size_));
      memset(buffer_, 0, size_);
  }
}


/**
 * \fn void nimbus::PhysBAMData::Destroy()
 * \brief Brief description.
 * \return
*/
void PhysBAMData::Destroy() {
  if (buffer_) {
    delete [] buffer_;
    buffer_ = NULL;
  }
  size_ = 0;
}


/**
 * \fn void nimbus::PhysBAMData::Copy(Data *from)
 * \brief Brief description.
 * \param from
 * \return
*/
void PhysBAMData::Copy(Data *from) {
  Destroy();
  PhysBAMData* pfrom = static_cast<PhysBAMData*>(from);
  size_ = pfrom->size();
  buffer_ = static_cast<char*>(malloc(size_));
  memcpy(buffer_, pfrom->buffer(), size_);
  hash = pfrom->hash;
}


/**
 * \fn bool nimbus::PhysBAMData::Serialize(SerializedData *ser_data)
 * \brief Brief description.
 * \param ser_data
 * \return
*/
bool PhysBAMData::Serialize(SerializedData *ser_data) {
  nimbus_message::pd_message pd; // NOLINT
  if (buffer_) {
      std::string buf(buffer_, size_);
      pd.set_buffer(buf);
  }
  if (size_)
      pd.set_size(size_);
  else
      pd.set_size(0);
  if (hash != this->HashCode()) {
    dbg(DBG_ERROR, "Data %s got corrupted somewhere before serialization!!\n", name().c_str());
    //assert(false);
  }
  pd.set_hash(hash);
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


/**
 * \fn bool nimbus::PhysBAMData::DeSerialize(const SerializedData &ser_data,
                                 Data **result)
 * \brief This function does not free buffer. Destroy should be called to free
 * buffer first, otherwise there will be a memory leak.
 * \param ser_data
 * \param result
 * \return
*/
bool PhysBAMData::DeSerialize(const SerializedData &ser_data,
                              Data **result) {
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
    PhysBAMData *data = new PhysBAMData();
    (*result) = static_cast<Data *>(data);
    if (pd.has_buffer()) {
        int size = pd.size();
        char *buffer = new char[size];
        memcpy(buffer, pd.buffer().c_str(), sizeof(char) * size); // NOLINT
        data->set_buffer(buffer, size);
    } else {
        data->set_buffer(NULL, pd.size());
    }
    if (pd.hash() != data->HashCode()) {
      dbg(DBG_ERROR, "For %s data sent != data received!!!\n", name().c_str());
      //assert(false);
    } else {
      hash = pd.hash();
    }
    return success;
}

/** Clear out the data from the temporary buffer. Note that this will
 * lose any uncommitted data. */
void PhysBAMData::ClearTempBuffer() {
  if (temp_buffer_) {
    delete temp_buffer_;
  }
  temp_buffer_ = new std::stringstream();
}

bool PhysBAMData::AddToTempBuffer(char* buffer, int len) {
  temp_buffer_->write(buffer, len);
  return true;
}

int PhysBAMData::CommitTempBuffer() {
  int len = temp_buffer_->tellp();
  if (buffer_)
    delete buffer_;
  size_ = len;
  buffer_ = static_cast<char*>(malloc(len));
  temp_buffer_->read(buffer_, len);
  if (temp_buffer_->eof()) {
    dbg(DBG_WARN, "When copying a temporary buffer into the permanent buffer in a PhysBAMData object, the read was incomplete.\n");  // NOLINT
  }
  return len;
}

}  // namespace nimbus
