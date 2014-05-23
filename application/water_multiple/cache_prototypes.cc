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
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <string>

#include "application/water_multiple/cache_data_include.h"
#include "application/water_multiple/cache_prototypes.h"
#include "application/water_multiple/parameters.h"

namespace application {

nimbus::GeometricRegion zero_reg;

CacheFaceArray<T> kCacheFaceVel(kDefaultRegion, zero_reg);
CacheFaceArray<T> kCacheFaceVelGhost(kDefaultRegion, zero_reg, kGhostNum);
CacheFaceArray<bool> kCachePsiN(kDefaultRegion, zero_reg, 1);

CacheScalarArray<T> kCachePhi3(kDefaultRegion, zero_reg, 3);
CacheScalarArray<T> kCachePhi7(kDefaultRegion, zero_reg, 7);
CacheScalarArray<T> kCachePhi8(kDefaultRegion, zero_reg, 8);
CacheScalarArray<bool> kCachePsiD(kDefaultRegion, zero_reg, 1);

CacheParticleLevelsetEvolution<float> kCachePLE(kDefaultRegion, zero_reg, 3);

} // namespace application

