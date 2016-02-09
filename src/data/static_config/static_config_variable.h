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
 * The base class for a static config variable. A static config variable is
 * a simulation variable that is not changed after initialization, so Nimbus
 * tries to cache static config variables so that they are constructed only once
 * during a simulation.
 *
 * Author: Hang Qu <quhang@stanford.edu>
 */

#ifndef NIMBUS_SRC_DATA_STATIC_CONFIG_STATIC_CONFIG_VARIABLE_H_
#define NIMBUS_SRC_DATA_STATIC_CONFIG_STATIC_CONFIG_VARIABLE_H_

#include "src/shared/geometric_region.h"
#include "src/shared/nimbus_types.h"

namespace nimbus {

class StaticConfigVariable {
 public:
  virtual StaticConfigVariable* CreateNew(const GeometricRegion& local_region)
      const = 0;
  virtual void Destroy() = 0;
  virtual int ChangeUserCount(int delta) {
    users_ += delta;
    return users_;
  }
 private:
  int users_;
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_DATA_STATIC_CONFIG_STATIC_CONFIG_VARIABLE_H_
