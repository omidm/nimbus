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
  * A PhysBAM physical data instance. This is simply a buffer of bytes.
  * An application uses a TranslatorPhysBAM to read from a PhysBAMData
  * to generate a PhysBAM object at the beginning of a job, and uses a
  * TranslatorPhysBAM to store a PhysBAM object into a PhysBAMData at the
  * end of a job.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_DATA_PHYSBAM_PHYSBAM_DATA_H_
#define NIMBUS_DATA_PHYSBAM_PHYSBAM_DATA_H_

#include "shared/cluster.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"
#include <sstream>  // NOLINT

namespace nimbus {

class PhysBAMData: public Data {
 private:
  int_dimension_t size_;
  char* buffer_;
  std::stringstream* temp_buffer_;

 public:
  PhysBAMData();
  virtual ~PhysBAMData() {}

  virtual Data* Clone();
  virtual void Create();
  virtual void Destroy();

  virtual void Copy(Data* from);
  virtual bool Serialize(SerializedData* ser_data);
  virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

  virtual char* buffer() {return buffer_;}
  virtual void set_buffer(char* b, int_dimension_t s);

  virtual int_dimension_t size();
  virtual void set_size(int_dimension_t s);

  virtual void ClearTempBuffer();
  virtual bool AddToTempBuffer(char* buffer, int len);
  virtual int CommitTempBuffer();


  // Not implemented, not clear what these interfaces mean. -pal
  virtual void duplicate(Computer source, Computer destination) {}
  virtual void migrate(Computer sourcer, Computer destination) {}
  virtual void split(Data *, Data *) {}
  virtual void merge(Data *, Data *) {}
};

}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_PHYSBAM_DATA_H_
