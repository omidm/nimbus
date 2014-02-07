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
 * Scalar data.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <sstream>
#include <string>
#include "data/scalar_data.h"

namespace nimbus {

    template<class T> ScalarData<T>::ScalarData() {}

    template<class T> ScalarData<T>::ScalarData(std::string name) {
        set_name(name);
    }

    template<class T> Data* ScalarData<T>::Clone() {
        return new ScalarData<T>(name());
    }

    template<class T> void ScalarData<T>::Create() {}

    template<class T> void ScalarData<T>::Destroy() {}

    template<class T> void ScalarData<T>::Copy(Data* from) {
        ScalarData<T> *sfrom = static_cast<ScalarData<T> * >(from);
        scalar_ = sfrom->scalar();
    }

    template<class T> bool ScalarData<T>::Serialize(SerializedData* ser_data) {
        std::ostringstream ser;
        ser << scalar_;
        const char *temp = ser.str().c_str();
        size_t len = ser.str().size() + 1;
        char *buf = new char[len];
        memcpy(buf, temp, len);
        ser_data->set_data_ptr(buf);
        ser_data->set_size(ser.str().size()+1); // NOLINT
        return true;
    }

    template <class T> bool ScalarData<T>::
    DeSerialize(const SerializedData &ser_data, Data **result) {
        std::string str(ser_data.data_ptr_raw(), ser_data.size());
        std::istringstream ser(str);
        ScalarData<T> * sd  = new ScalarData<T>();
        *result = sd;
        T temp;
        ser >> temp;
        sd->set_scalar(temp);
        return true;
    }

    template<class T> void ScalarData<T>::set_scalar(T scalar) {
        scalar_ = scalar;
    }

    template<class T> T ScalarData<T>::scalar() {
        return scalar_;
    }

}  // namespace nimbus

template class nimbus::ScalarData<int>;
