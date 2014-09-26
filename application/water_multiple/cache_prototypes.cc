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

CacheFaceArray<T> kCacheFaceVel(kDefaultRegion, 0, true, "face_vel");
CacheFaceArray<T> kCacheFaceVelGhost(kDefaultRegion, 3, true, "face_vel_ghost");
CacheFaceArray<bool> kCachePsiN(kDefaultRegion, 1, true, "psi_n");

CacheScalarArray<T> kCachePhi3(kDefaultRegion, 3, true, "phi_3");
CacheScalarArray<T> kCachePhi7(kDefaultRegion, 7, true, "phi_7");
CacheScalarArray<T> kCachePhi8(kDefaultRegion, 8, true, "phi_8");
CacheScalarArray<bool> kCachePsiD(kDefaultRegion, 1, true, "psi_d");

// Varibales for projection.
CacheScalarArray<T> kCachePressure(kDefaultRegion, 1, true, "pressure");
CacheScalarArray<int> kCacheColors(kDefaultRegion, 1, true,
                                   "filled_region_colors");
CacheScalarArray<T> kCacheDivergence(kDefaultRegion, 1, true, "divergence");
// TODO(quhang): this cache variable is questionable, because it cannot be
// deleted if meta_p is being used.
CacheRawGridArray kCacheIndexC2M(kDefaultRegion, true, "index_c2m");

CacheParticleLevelsetEvolution<float> kCachePLE(kDefaultRegion, 3, true,
                                                "particle_container");

CacheSparseMatrix kCacheSparseMatrixA(kDefaultRegion, true, "matrix_c");
CacheSparseMatrix kCacheSparseMatrixC(kDefaultRegion, true, "matrix_a");

CacheArrayM2C kCacheArrayM2C(kDefaultRegion, true, "index_m2c");

CacheCompressedScalarArray<float> kCacheMetaP(kDefaultRegion, 1, true,
                                              "vector_p_meta_format");

CacheVector kCacheVectorB(kDefaultRegion, true, "vector_b");
CacheVector kCacheVectorPressure(kDefaultRegion, true, "vector_pressure");
CacheVector kCacheVectorZ(kDefaultRegion, true, "vector_z");
CacheVector kCacheVectorTemp(kDefaultRegion, true, "vector_temp");
} // namespace application

