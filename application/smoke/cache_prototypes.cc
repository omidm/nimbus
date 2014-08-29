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
 * Modifier for smoke: Andrew Lim <alim16@stanford.edu> 
 */

#include "application/smoke/cache_prototypes.h"
#include "application/smoke/parameters.h"

namespace application {

CacheFaceArray<T> kCacheFaceVel(kDefaultRegion, 0, true);
CacheFaceArray<T> kCacheFaceVelGhost(kDefaultRegion, 3, true);
CacheFaceArray<bool> kCachePsiN(kDefaultRegion, 1, true);

CacheScalarArray<T> kCacheDensity(kDefaultRegion, 0, true);
CacheScalarArray<T> kCacheDensityGhost(kDefaultRegion, 3, true);

CacheScalarArray<bool> kCachePsiD(kDefaultRegion, 1, true);

// Varibales for projection.
CacheScalarArray<T> kCachePressure(kDefaultRegion, 1, true);
CacheScalarArray<T> kCacheVectorPGridFormat(kDefaultRegion, 1, true);
CacheScalarArray<int> kCacheColors(kDefaultRegion, 1, true);
CacheScalarArray<T> kCacheDivergence(kDefaultRegion, 1, true);

CacheSparseMatrix kCacheSparseMatrixA(kDefaultRegion, true);
CacheSparseMatrix kCacheSparseMatrixC(kDefaultRegion, true);

CacheArrayM2C kCacheArrayM2C(kDefaultRegion, true);

} // namespace application

