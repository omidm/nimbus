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
 * An Example application that is meant to run over single worker.
 * It is simply applying a stencil over a one dimensional array.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_STENCIL_1D_SINGLE_APP_H_
#define NIMBUS_APPLICATION_STENCIL_1D_SINGLE_APP_H_

#include <boost/tokenizer.hpp>
#include <string>
#include <iostream> // NOLINT
#include "worker/application.h"
#include "worker/job.h"
#include "worker/data.h"
#include "shared/parser.h"
#include "shared/nimbus_types.h"
#include "protobuf_compiled/vector_msg.pb.h"

#define ML 4
#define GL 1
#define LOOP_COUNTER 15
#define LOOP_CONDITION 0

using nimbus::Job;
using nimbus::Data;
using nimbus::Application;
using boost::tokenizer;
using boost::char_separator;


class Main : public Job {
  public:
    explicit Main(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Init : public Job {
  public:
    Init();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Print : public Job {
  public:
    Print();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class ApplyLeft : public Job {
  public:
    ApplyLeft();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class ApplyRight : public Job {
  public:
    ApplyRight();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class UpdateLeft : public Job {
  public:
    UpdateLeft();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class UpdateRight : public Job {
  public:
    UpdateRight();
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class ForLoop : public Job {
  public:
    explicit ForLoop(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job * Clone();
};

class Vec : public Data {
  public:
    explicit Vec(int size);
    virtual ~Vec();

    virtual void Create();
    virtual void Destroy();
    virtual Data * Clone();
    virtual void Copy(Data* from);
    virtual bool Serialize(SerializedData* ser_data);
    virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

    int size();
    int* arr();

  private:
    int size_;
    int *arr_;
};

class App : public Application {
  public:
    App();
    virtual void Load();
};











#endif  // NIMBUS_APPLICATION_STENCIL_1D_SINGLE_APP_H_
