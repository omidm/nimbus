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
  * Object representation of the parameter field.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_PARAMETER_H_
#define NIMBUS_SHARED_PARAMETER_H_

#include <boost/tokenizer.hpp>
#include <string>
#include "shared/nimbus_types.h"
#include "shared/idset.h"
#include "shared/serialized_data.h"



namespace nimbus {

class Parameter {
  public:
    Parameter();
    explicit Parameter(const SerializedData& ser_data);
    Parameter(const Parameter& other);
    virtual ~Parameter();

    SerializedData ser_data();

    void set_ser_data(SerializedData ser_data);

    std::string ToNetworkData();
    Parameter& operator= (const Parameter& right);

  private:
    SerializedData ser_data_;
    IDSet<param_id_t> idset_;
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_PARAMETER_H_
