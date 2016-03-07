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
 * Data classes for lofestic regression application.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#ifndef NIMBUS_APPLICATIONS_ML_LR_DATA_H_
#define NIMBUS_APPLICATIONS_ML_LR_DATA_H_

#include <assert.h>
#include <vector>
#include <string>
#include <iostream> // NOLINT
#include "src/shared/nimbus.h"
#include "protobuf_compiled/data_msgs.pb.h"

#define WEIGHT_DATA_NAME "wight"
#define SCRATCH_WEIGHT_DATA_NAME "scratch_wight"
#define SAMPLE_BATCH_DATA_NAME "sample_batch"

using namespace nimbus; // NOLINT


class Sample {
  public:
    explicit Sample(const size_t& dimension);
    virtual ~Sample();

    size_t dimension() const;
    double label() const;
    std::vector<double>* vector();

    void set_label(double label);

  private:
    size_t dimension_;
    double label_;
    std::vector<double> vector_;
};

class SampleBatch : public Data {
  public:
    SampleBatch(const size_t& dimension,
                const size_t& sample_num);
    virtual ~SampleBatch();

    virtual void Create();
    virtual void Destroy();
    virtual Data * Clone();
    virtual void Copy(Data* from);
    virtual bool Serialize(SerializedData* ser_data);
    virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

    size_t dimension() const;
    size_t sample_num() const;
    std::vector<Sample>* samples();

  private:
    size_t dimension_;
    size_t sample_num_;
    std::vector<Sample> samples_;
};

class Weight : public Data {
  public:
    explicit Weight(const size_t& dimension,
                    const std::string& name);
    virtual ~Weight();

    virtual void Create();
    virtual void Destroy();
    virtual Data * Clone();
    virtual void Copy(Data* from);
    virtual bool Serialize(SerializedData* ser_data);
    virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

    size_t dimension() const;
    std::vector<double>* vector();
    std::vector<double>* gradient();

    void set_vector(const std::vector<double>& vector);
    void set_gradient(const std::vector<double>& gradient);

  private:
    size_t dimension_;
    std::vector<double> vector_;
    std::vector<double> gradient_;
};


#endif  // NIMBUS_APPLICATIONS_ML_LR_DATA_H_
