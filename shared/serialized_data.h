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

#ifndef NIMBUS_SHARED_SERIALIZED_DATA_H_
#define NIMBUS_SHARED_SERIALIZED_DATA_H_


#include <boost/shared_ptr.hpp>
#include <list>
#include <map>
#include <string>
#include "shared/escaper.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class SerializedData {
  public:
    SerializedData();
    SerializedData(const boost::shared_ptr<char>& data_ptr, const size_t& size);
    SerializedData(const SerializedData& other);
    ~SerializedData();

    size_t size() const;
    char* data_ptr_raw() const;
    boost::shared_ptr<char> data_ptr() const;

    void set_size(size_t size);
    void set_data_ptr(char* ptr);
    void set_data_ptr(boost::shared_ptr<char> ptr);

    std::string toString();
    SerializedData& operator= (const SerializedData& right);

  private:
    boost::shared_ptr<char> data_ptr_;
    size_t size_;
};

typedef std::list<SerializedData*> SerializedDataList;
typedef std::map<job_id_t, SerializedData*> SerializedDataMap;

}  // namespace nimbus

#endif  // NIMBUS_SHARED_SERIALIZED_DATA_H_
