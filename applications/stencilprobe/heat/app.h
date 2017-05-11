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
 * Application interface.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#ifndef NIMBUS_APPLICATIONS_SCAFFOLD_APP_H_
#define NIMBUS_APPLICATIONS_SCAFFOLD_APP_H_

#include <iostream> // NOLINT
#include <string>
#include "src/worker/application.h"
#include "src/shared/nimbus_types.h"
#include "src/application_utils/partition_handler.h"
#include "./app_data.h"

using nimbus::Application;

extern nimbus::PartitionHandler ph;
extern AppDataVec AppDataVecPrototype;
extern GeometricRegion GlobalRegion;

class Heat : public Application {
  public:
    Heat(size_t iter_num,
         int bw,
         int nx, int ny, int nz,
         int pnx, int pny, int pnz);
    ~Heat();
    virtual void Load();

    size_t iter_num_;
    int bw_;
    int nx_;
    int ny_;
    int nz_;
    int pnx_;
    int pny_;
    int pnz_;
};



#endif  // NIMBUS_APPLICATIONS_SCAFFOLD_APP_H_
