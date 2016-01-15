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

#ifndef NIMBUS_APPLICATION_K_MEANS_DATA_H_
#define NIMBUS_APPLICATION_K_MEANS_DATA_H_

#include <assert.h>
#include <vector>
#include <string>
#include <iostream> // NOLINT
#include "shared/nimbus.h"
#include "protobuf_compiled/data_msgs.pb.h"

#define MEANS_DATA_NAME "means"
#define SAMPLE_BATCH_DATA_NAME "sample_batch"

using namespace nimbus; // NOLINT


class Sample {
  public:
    explicit Sample(const size_t& dimension);
    virtual ~Sample();

    size_t dimension() const;
    std::vector<double>* vector();

  private:
    size_t dimension_;
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


class Mean {
  public:
    explicit Mean(const size_t& dimension);
    virtual ~Mean();

    size_t dimension() const;
    std::vector<double>* vector();
    std::vector<double>* scratch();
    size_t scratch_weight() const;

    void set_scratch_weight(const size_t& w);

  private:
    size_t dimension_;
    std::vector<double> vector_;
    std::vector<double> scratch_;
    size_t scratch_weight_;
};

class Means : public Data {
  public:
    explicit Means(const size_t& dimension,
                   const size_t& cluster_num);
    virtual ~Means();

    virtual void Create();
    virtual void Destroy();
    virtual Data * Clone();
    virtual void Copy(Data* from);
    virtual bool Serialize(SerializedData* ser_data);
    virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

    size_t dimension() const;
    size_t cluster_num() const;
    std::vector<Mean>* means();

  private:
    size_t dimension_;
    size_t cluster_num_;
    std::vector<Mean> means_;
};


#endif  // NIMBUS_APPLICATION_K_MEANS_DATA_H_
