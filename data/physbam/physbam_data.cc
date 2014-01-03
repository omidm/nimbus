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
#include "data/physbam/physbam_data.h"

namespace nimbus {
/**
 * \fn nimbus::PhysBAMData::PhysBAMData()
 * \brief Brief description.
 * \return
*/
PhysBAMData::PhysBAMData(): size_(0), buffer_(0), temp_buffer_(0) {}


/**
 * \fn Data * nimbus::PhysBAMData::Clone()
 * \brief Brief description.
 * \return
*/
Data * PhysBAMData::Clone() {
  PhysBAMData* d = new PhysBAMData();
  char *buf = NULL;
  if (size_)
    buf = static_cast<char*>(malloc(size_));
  if (size_ && buffer_)
    memcpy(buf, buffer_, size_);
  d->set_buffer(buf, size_);
  return d;
}


/**
 * \fn void nimbus::PhysBAMData::Create()
 * \brief Brief description.
 * \return
*/
void PhysBAMData::Create() {
  if (size_ && !buffer_)
      buffer_ = static_cast<char*>(malloc(size_));
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
}


/**
 * \fn bool nimbus::PhysBAMData::Serialize(SerializedData *ser_data)
 * \brief Brief description.
 * \param ser_data
 * \return
*/
bool PhysBAMData::Serialize(SerializedData *ser_data) {
  return false;
}


/**
 * \fn bool nimbus::PhysBAMData::DeSerialize(const SerializedData &ser_data,
                                 Data **result)
 * \brief Brief description.
 * \param ser_data
 * \param result
 * \return
*/
bool PhysBAMData::DeSerialize(const SerializedData &ser_data,
                              Data **result) {
  return false;
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
  // The usage of read/write/tellp seems incorrect.  -quhang
  // int len = temp_buffer_->tellp();
  int len = temp_buffer_->str().size();
  if (buffer_)
    delete buffer_;
  size_ = len;
  buffer_ = static_cast<char*>(malloc(len));
  temp_buffer_->read(buffer_, len);
  temp_buffer_->clear();
  if (!temp_buffer_->str().empty()) {
    dbg(DBG_WARN, "When copying a temporary buffer into the permanent buffer in a PhysBAMData object, the read was incomplete. %i bytes remaining.\n", temp_buffer_->tellp());  // NOLINT
  }
  return len;
}

}  // namespace nimbus
