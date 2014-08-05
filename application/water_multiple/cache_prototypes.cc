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

#include "application/water_multiple/cache_prototypes.h"
#include "application/water_multiple/parameters.h"

namespace application {

CacheFaceArray<T> kCacheFaceVel(kDefaultRegion, 0, true);
CacheFaceArray<T> kCacheFaceVelGhost(kDefaultRegion, 3, true);
CacheFaceArray<bool> kCachePsiN(kDefaultRegion, 1, true);

CacheScalarArray<T> kCachePhi3(kDefaultRegion, 3, true);
CacheScalarArray<T> kCachePhi7(kDefaultRegion, 7, true);
CacheScalarArray<T> kCachePhi8(kDefaultRegion, 8, true);
CacheScalarArray<bool> kCachePsiD(kDefaultRegion, 1, true);

// Varibales for projection.
CacheScalarArray<T> kCachePressure(kDefaultRegion, 1, true);
CacheScalarArray<int> kCacheColors(kDefaultRegion, 1, true);
CacheScalarArray<T> kCacheDivergence(kDefaultRegion, 1, true);
// TODO(quhang): this cache variable is questionable, because it cannot be
// deleted if meta_p is being used.
CacheScalarArray<int> kCacheArrayC2M(kDefaultRegion, 0, true);

CacheParticleLevelsetEvolution<float> kCachePLE(kDefaultRegion, 3, true);

CacheSparseMatrix kCacheSparseMatrixA(kDefaultRegion, true);
CacheSparseMatrix kCacheSparseMatrixC(kDefaultRegion, true);

CacheArrayM2C kCacheArrayM2C(kDefaultRegion, true);

CacheCompressedScalarArray<float> kCacheMetaP(kDefaultRegion, 1, true);

CacheVector kCacheVectorB(kDefaultRegion, true);
CacheVector kCacheVectorPressure(kDefaultRegion, true);
CacheVector kCacheVectorZ(kDefaultRegion, true);
CacheVector kCacheVectorTemp(kDefaultRegion, true);
} // namespace application

