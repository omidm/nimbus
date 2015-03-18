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
 * TranslatorPhysBAM is a class for translating Nimbus physical objects into
 * PhysBAM data objects. It is templated to make it a little easier
 * to handle PhysBAM's generality, taking a simple template parameter
 * of a PhysBAM VECTOR. This VECTOR is typically a 2D or 3D float: it
 * represents a point in space. The class derives the scalar type
 * (typically float) from this VECTOR, as well as the dimensionality.
 *
 * Author: Philip Levis <pal@cs.stanford.edu>,
 *         Chinmayee Shah <chshah@stanford.edu>,
 *         Hang Qu <quhang@stanford.edu>
 */

#ifndef NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
#define NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_

#include <sys/syscall.h>
#include <algorithm>
#include <cmath>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "data/physbam/physbam_include.h"
#include "data/physbam/physbam_data.h"
#include "data/physbam/physbam_data_with_meta.h"

#include "shared/log.h"
#include "shared/fast_log.hh"
#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "worker/physical_data_instance.h"
#include "worker/worker.h"

#define TRANSLATE_LOG_H "[Translator]"

namespace nimbus {

template <class TS> class TranslatorPhysBAM {
    public:
        typedef typename PhysBAM::VECTOR<TS, 3> TV;
        typedef typename PhysBAM::VECTOR<int, 3> TV_INT;
        typedef typename PhysBAM::GRID<TV> Grid;
        typedef typename PhysBAM::FACE_INDEX<TV::dimension> FaceIndex;
        typedef typename PhysBAM::VECTOR<int_dimension_t, 3> Dimension3Vector;

        // Container class for particles and removed particles.
        typedef typename PhysBAM::PARTICLE_LEVELSET_UNIFORM<Grid> ParticleContainer;
        // Particle bucket. Each node has a linked list of particle buckets.
        typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> ParticleBucket;
        // Particle array, indexed by node.
        typedef typename PhysBAM::ARRAY<ParticleBucket*, TV_INT> ParticleArray;

        typedef typename PhysBAM::PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>
            RemovedParticleBucket;
        typedef typename PhysBAM::ARRAY<RemovedParticleBucket*, TV_INT>
            RemovedParticleArray;

        typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleLevelset;
        // typedef typename PhysBAM::ARRAY<TS, TV_INT> ScalarArray;

        enum {
            X_COORD = 1,
            Y_COORD = 2,
            Z_COORD = 3
        };

        explicit TranslatorPhysBAM() {}
        virtual ~TranslatorPhysBAM() {}

        // Data structures used to format particles in PhysBAMData.
        // Should be changed to protocol buffer later for compatibility.
        // TODO(quhang) .
        struct ParticleInternal {
            TS position[3];
            TS radius;
            uint16_t quantized_collision_distance;
            int32_t id;
        };

        struct RemovedParticleInternal : public ParticleInternal {
            TS v[3];
        };

        static Log *log;

        /** Take a FaceArray described by region and read its data from the
         *  PhysicalDataInstance objects in the objects array.
         */
        template<typename T> static void ReadFaceArray(
                const GeometricRegion &region,
                const GeometricRegion &inner,
                const Coord &shift,
                const DataArray &read_set,
                typename PhysBAM::ARRAY<T, FaceIndex>* fa) {
            timer::StartTimer(timer::k1);
            ReadFaceArrayInner(region, inner, shift, read_set, fa);
            timer::StopTimer(timer::k1);
        }

        /** Take a FaceArray described by region and read its data from the
         *  PhysicalDataInstance objects in the objects array.
         */
        template<typename T> static void ReadFaceArrayInner(
                const GeometricRegion &region,
                const GeometricRegion &inner,
                const Coord &shift,
                const DataArray &read_set,
                typename PhysBAM::ARRAY<T, FaceIndex>* fa) {
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Face Array (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (read_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Read Face Array (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                return;
            }

            size_t read_size = read_set.size();

            // mark data that must be read in phase 1 and phase 2 for each
            // dimension
            std::vector<bool> read_flag(4*read_size, false);
            for (size_t i = 0; i < read_set.size(); ++i) {
                Data *data = read_set[i];
                const GeometricRegion dregion = data->region();
                Dimension3Vector overlap = GetOverlapSize(dregion, region);
                if (HasOverlap(overlap)) {
                    read_flag[4*i] = true;
                    // hard coded ghost width here to determine phases
                    read_flag[4*i+X_COORD] = (dregion.dx() != 3);  // read vx in phase 1
                    read_flag[4*i+Y_COORD] = (dregion.dy() != 3);  // read vy in phase 1
                    read_flag[4*i+Z_COORD] = (dregion.dz() != 3);  // read vz in phase 1
                }
            }

            for (size_t i = 0; i < read_set.size(); ++i) {
                // read anything in phase 1?
                if (read_flag[4*i] &&
                        (read_flag[4*i+X_COORD] || read_flag[4*i+Y_COORD] || read_flag[4*i+Z_COORD])) {  // NOLINT
                    // setup
                    PhysBAMData *data = static_cast<PhysBAMData *>(read_set[i]);
                    T* buffer = reinterpret_cast<T*>(data->buffer());
                    const GeometricRegion dregion = data->region();
                    const int_dimension_t ddx = dregion.dx();
                    const int_dimension_t ddy = dregion.dy();
                    const int_dimension_t ddz = dregion.dz();
                    const Dimension3Vector src = GetOffset(dregion, region);
                    const Dimension3Vector dest = GetOffset(region, dregion);
                    const Dimension3Vector overlap = GetOverlapSize(dregion, region);
                    for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                        if (!read_flag[4*i+dim])
                            continue;  // read in phase 2
                        // setup for vx/ vy/ vz
                        int range_x, range_y, range_z;
                        int mult_x, mult_y;
                        int src_offset;
                        switch (dim) {
                            case X_COORD:
                                range_x = overlap(X_COORD) + 1;
                                range_y = overlap(Y_COORD);
                                range_z = overlap(Z_COORD);
                                mult_x = ddy * ddz;
                                mult_y = ddz;
                                src_offset = 0;  // offset into buffer
                                break;
                            case Y_COORD:
                                range_x = overlap(X_COORD);
                                range_y = overlap(Y_COORD) + 1;
                                range_z = overlap(Z_COORD);
                                mult_x = (ddy + 1) * ddz;
                                mult_y = ddz;
                                src_offset = (ddx + 1) * ddy * ddz;  // offset into buffer
                                break;
                            case Z_COORD:
                                range_x = overlap(X_COORD);
                                range_y = overlap(Y_COORD);
                                range_z = overlap(Z_COORD) + 1;
                                mult_x = ddy * (ddz + 1);
                                mult_y = ddz + 1;
                                src_offset = (ddx + 1) * ddy * ddz +
                                             ddx * (ddy + 1) * ddz;  // offset into buffer
                                break;
                        }  // switch dim
                        for (int x = 0; x < range_x; ++x) {
                            for (int y = 0; y < range_y; ++y) {
                                int source_x = x + src(X_COORD);
                                int source_y = y + src(Y_COORD);
                                int source_index = source_x * mult_x +
                                                   source_y * mult_y +
                                                   src_offset;
                                int dest_x = x + dest(X_COORD) + region.x() - shift.x;
                                int dest_y = y + dest(Y_COORD) + region.y() - shift.y;
                                int dest_z = 0 + dest(Z_COORD) + region.z() - shift.z;
                                const TV_INT destination_index(dest_x, dest_y, dest_z);
                                memcpy(&((*fa)(dim, destination_index)),  // NOLINT
                                       &(buffer[source_index]),
                                       sizeof(T) * range_z);

                            }  // read y for loop
                        }  // read x for loop
                    }  // for dim
                }  //  if read anything
            }  //  outermost for loop over entired read set

            const int cx1 = inner.x();
            const int cx2 = inner.x() + inner.dx();
            const int cy1 = inner.y();
            const int cy2 = inner.y() + inner.dy();
            const int cz1 = inner.z();
            const int cz2 = inner.z() + inner.dz();

            for (size_t i = 0; i < read_set.size(); ++i) {
                // read anything in phase 1?
                if (read_flag[4*i] &&
                        (!read_flag[4*i+X_COORD] || !read_flag[4*i+Y_COORD] || !read_flag[4*i+Z_COORD])) {  // NOLINT
                    // setup
                    PhysBAMData *data = static_cast<PhysBAMData *>(read_set[i]);
                    T* buffer = reinterpret_cast<T*>(data->buffer());
                    const GeometricRegion dregion = data->region();
                    const Dimension3Vector src = GetOffset(dregion, region);
                    const Dimension3Vector dest = GetOffset(region, dregion);
                    const Dimension3Vector overlap = GetOverlapSize(dregion, region);
                    for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                        if (read_flag[4*i+dim])
                            continue;  // already read in phase 1
                        // case for vx/ vy/ vz
                        // exapnded out for performance -- I can't figure out a
                        // better way to handle averaging of common face values
                        // -- Chinmayee
                        switch (dim) {
                            case X_COORD: {
                                const int_dimension_t ddy = dregion.dy();
                                const int_dimension_t ddz = dregion.dz();
                                assert(overlap(X_COORD) == 3);
                                const int range_x = 3 + 1;
                                const int range_y = overlap(Y_COORD);
                                const int range_z = overlap(Z_COORD);
                                const int mult_x = ddy * ddz;
                                const int mult_y = ddz;
                                const int mult_z = 1;
                                const int src_offset = 0;  // offset into buffer
                                for (int x = 0; x < range_x; ++x) {
                                    for (int y = 0; y < range_y; ++y) {
                                        for (int z = 0; z < range_z; ++z) {
                                            int source_x = x + src(X_COORD);
                                            int source_y = y + src(Y_COORD);
                                            int source_z = z + src(Y_COORD);
                                            int source_index = source_x * mult_x +
                                                               source_y * mult_y +
                                                               source_z * mult_z +
                                                               src_offset;
                                            int loc_x = x + dest(X_COORD) + region.x();
                                            int loc_y = y + dest(Y_COORD) + region.y();
                                            int loc_z = z + dest(Z_COORD) + region.z();
                                            int dest_x = loc_x - shift.x;
                                            int dest_y = loc_y - shift.y;
                                            int dest_z = loc_z - shift.z;
                                            const TV_INT destination_index(dest_x, dest_y, dest_z);  // NOLINT
                                            if (loc_x == cx1 || loc_x == cx2) {
                                                (*fa)(dim, destination_index) += buffer[source_index];  // NOLINT
                                                (*fa)(dim, destination_index) /= 2;
                                            } else {
                                                (*fa)(dim, destination_index) = buffer[source_index];  // NOLINT
                                            }  // branch for averaging
                                        }  // read z for loop
                                    }  // read y for loop
                                }  // read x for loop
                                break;
                            }
                            case Y_COORD: {
                                const int_dimension_t ddx = dregion.dx();
                                const int_dimension_t ddy = 3;
                                const int_dimension_t ddz = dregion.dz();
                                assert(overlap(Y_COORD) == 3);
                                const int range_x = overlap(X_COORD);
                                const int range_y = 3 + 1;
                                const int range_z = overlap(Z_COORD);
                                const int mult_x = (ddy + 1) * ddz;
                                const int mult_y = ddz;
                                const int mult_z = 1;
                                const int src_offset = (ddx + 1) * ddy * ddz;  // offset into buffer
                                for (int x = 0; x < range_x; ++x) {
                                    for (int y = 0; y < range_y; ++y) {
                                        for (int z = 0; z < range_z; ++z) {
                                            int source_x = x + src(X_COORD);
                                            int source_y = y + src(Y_COORD);
                                            int source_z = z + src(Y_COORD);
                                            int source_index = source_x * mult_x +
                                                               source_y * mult_y +
                                                               source_z * mult_z +
                                                               src_offset;
                                            int loc_x = x + dest(X_COORD) + region.x();
                                            int loc_y = y + dest(Y_COORD) + region.y();
                                            int loc_z = z + dest(Z_COORD) + region.z();
                                            int dest_x = loc_x - shift.x;
                                            int dest_y = loc_y - shift.y;
                                            int dest_z = loc_z - shift.z;
                                            const TV_INT destination_index(dest_x, dest_y, dest_z);  // NOLINT
                                            if (loc_y == cy1 || loc_y == cy2) {
                                                (*fa)(dim, destination_index) += buffer[source_index];  // NOLINT
                                                (*fa)(dim, destination_index) /= 2;
                                            } else {
                                                (*fa)(dim, destination_index) = buffer[source_index];  // NOLINT
                                            }  // branch for averaging
                                        }  // read z for loop
                                    }  // read y for loop
                                }  // read x for loop
                                break;
                            }
                            case Z_COORD: {
                                const int_dimension_t ddx = dregion.dx();
                                const int_dimension_t ddy = dregion.dy();
                                const int_dimension_t ddz = 3;
                                assert(overlap(Z_COORD) == 3);
                                const int range_x = overlap(X_COORD);
                                const int range_y = overlap(Y_COORD);
                                const int range_z = 3 + 1;
                                const int mult_x = ddy * (ddz + 1);
                                const int mult_y = ddz + 1;
                                const int mult_z = 1;
                                const int src_offset = (ddx + 1) * ddy * ddz +
                                                       ddx * (ddy + 1) * ddz;  // offset into buffer
                                for (int x = 0; x < range_x; ++x) {
                                    for (int y = 0; y < range_y; ++y) {
                                        for (int z = 0; z < range_z; ++z) {
                                            int source_x = x + src(X_COORD);
                                            int source_y = y + src(Y_COORD);
                                            int source_z = z + src(Y_COORD);
                                            int source_index = source_x * mult_x +
                                                               source_y * mult_y +
                                                               source_z * mult_z +
                                                               src_offset;
                                            int loc_x = x + dest(X_COORD) + region.x();
                                            int loc_y = y + dest(Y_COORD) + region.y();
                                            int loc_z = z + dest(Z_COORD) + region.z();
                                            int dest_x = loc_x - shift.x;
                                            int dest_y = loc_y - shift.y;
                                            int dest_z = loc_z - shift.z;
                                            const TV_INT destination_index(dest_x, dest_y, dest_z);  // NOLINT
                                            if (loc_z == cz1 || loc_z == cz2) {
                                                (*fa)(dim, destination_index) += buffer[source_index];  // NOLINT
                                                (*fa)(dim, destination_index) /= 2;
                                            } else {
                                                (*fa)(dim, destination_index) = buffer[source_index];  // NOLINT
                                            }  // branch for averaging
                                        }  // read z for loop
                                    }  // read y for loop
                                }  // read x for loop
                                break;
                            }
                        }  // switch dim
                    }  // for dim
                }  //  if read anything
            }  //  outermost for loop over entired read set

            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Face Array (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
        }

        /** Take a FaceArray described by region and read its data from the
         *  PhysicalDataInstance objects in the objects array.
         */
        static void ReadFaceArrayInner(
                const GeometricRegion &region,
                const GeometricRegion &inner,
                const Coord &shift,
                const DataArray &read_set,
                typename PhysBAM::ARRAY<bool, FaceIndex>* fa) {
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Face Array (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (read_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Read Face Array (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                return;
            }

            size_t read_size = read_set.size();

            // mark data that must be read in phase 1 and phase 2 for each
            // dimension
            std::vector<bool> read_flag(4*read_size, false);
            for (size_t i = 0; i < read_set.size(); ++i) {
                Data *data = read_set[i];
                const GeometricRegion dregion = data->region();
                Dimension3Vector overlap = GetOverlapSize(dregion, region);
                if (HasOverlap(overlap)) {
                    read_flag[4*i] = true;
                    // hard coded ghost width here to determine phases
                    read_flag[4*i+X_COORD] = (dregion.dx() != 3);  // read vx in phase 1
                    read_flag[4*i+Y_COORD] = (dregion.dy() != 3);  // read vy in phase 1
                    read_flag[4*i+Z_COORD] = (dregion.dz() != 3);  // read vz in phase 1
                }
            }

            for (size_t i = 0; i < read_set.size(); ++i) {
                // read anything in phase 1?
                if (read_flag[4*i] &&
                        (read_flag[4*i+X_COORD] || read_flag[4*i+Y_COORD] || read_flag[4*i+Z_COORD])) {  // NOLINT
                    // setup
                    PhysBAMData *data = static_cast<PhysBAMData *>(read_set[i]);
                    bool* buffer = reinterpret_cast<bool*>(data->buffer());
                    const GeometricRegion dregion = data->region();
                    const int_dimension_t ddx = dregion.dx();
                    const int_dimension_t ddy = dregion.dy();
                    const int_dimension_t ddz = dregion.dz();
                    const Dimension3Vector src = GetOffset(dregion, region);
                    const Dimension3Vector dest = GetOffset(region, dregion);
                    const Dimension3Vector overlap = GetOverlapSize(dregion, region);
                    for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                        if (!read_flag[4*i+dim])
                            continue;  // read in phase 2
                        // setup for vx/ vy/ vz
                        int range_x, range_y, range_z;
                        int mult_x, mult_y;
                        int src_offset;
                        switch (dim) {
                            case X_COORD:
                                range_x = overlap(X_COORD) + 1;
                                range_y = overlap(Y_COORD);
                                range_z = overlap(Z_COORD);
                                mult_x = ddy * ddz;
                                mult_y = ddz;
                                src_offset = 0;  // offset into buffer
                                break;
                            case Y_COORD:
                                range_x = overlap(X_COORD);
                                range_y = overlap(Y_COORD) + 1;
                                range_z = overlap(Z_COORD);
                                mult_x = (ddy + 1) * ddz;
                                mult_y = ddz;
                                src_offset = (ddx + 1) * ddy * ddz;  // offset into buffer
                                break;
                            case Z_COORD:
                                range_x = overlap(X_COORD);
                                range_y = overlap(Y_COORD);
                                range_z = overlap(Z_COORD) + 1;
                                mult_x = ddy * (ddz + 1);
                                mult_y = ddz + 1;
                                src_offset = (ddx + 1) * ddy * ddz +
                                             ddx * (ddy + 1) * ddz;  // offset into buffer
                                break;
                        }  // switch dim
                        for (int x = 0; x < range_x; ++x) {
                            for (int y = 0; y < range_y; ++y) {
                                int source_x = x + src(X_COORD);
                                int source_y = y + src(Y_COORD);
                                int source_index = source_x * mult_x +
                                                   source_y * mult_y +
                                                   src_offset;
                                int dest_x = x + dest(X_COORD) + region.x() - shift.x;
                                int dest_y = y + dest(Y_COORD) + region.y() - shift.y;
                                int dest_z = 0 + dest(Z_COORD) + region.z() - shift.z;
                                const TV_INT destination_index(dest_x, dest_y, dest_z);
                                memcpy(&((*fa)(dim, destination_index)),  // NOLINT
                                       &(buffer[source_index]),
                                       sizeof(bool) * range_z);

                            }  // read y for loop
                        }  // read x for loop
                    }  // for dim
                }  //  if read anything
            }  //  outermost for loop over entired read set

            const int cx1 = inner.x();
            const int cx2 = inner.x() + inner.dx();
            const int cy1 = inner.y();
            const int cy2 = inner.y() + inner.dy();
            const int cz1 = inner.z();
            const int cz2 = inner.z() + inner.dz();

            for (size_t i = 0; i < read_set.size(); ++i) {
                // read anything in phase 1?
                if (read_flag[4*i] &&
                        (!read_flag[4*i+X_COORD] || !read_flag[4*i+Y_COORD] || !read_flag[4*i+Z_COORD])) {  // NOLINT
                    // setup
                    PhysBAMData *data = static_cast<PhysBAMData *>(read_set[i]);
                    bool* buffer = reinterpret_cast<bool*>(data->buffer());
                    const GeometricRegion dregion = data->region();
                    const Dimension3Vector src = GetOffset(dregion, region);
                    const Dimension3Vector dest = GetOffset(region, dregion);
                    const Dimension3Vector overlap = GetOverlapSize(dregion, region);
                    for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                        if (read_flag[4*i+dim])
                            continue;  // already read in phase 1
                        // case for vx/ vy/ vz
                        // exapnded out for performance -- I can't figure out a
                        // better way to handle averaging of common face values
                        // -- Chinmayee
                        switch (dim) {
                            case X_COORD: {
                                const int_dimension_t ddy = dregion.dy();
                                const int_dimension_t ddz = dregion.dz();
                                assert(overlap(X_COORD) == 3);
                                const int range_x = 3 + 1;
                                const int range_y = overlap(Y_COORD);
                                const int range_z = overlap(Z_COORD);
                                const int mult_x = ddy * ddz;
                                const int mult_y = ddz;
                                const int mult_z = 1;
                                const int src_offset = 0;  // offset into buffer
                                for (int x = 0; x < range_x; ++x) {
                                    for (int y = 0; y < range_y; ++y) {
                                        for (int z = 0; z < range_z; ++z) {
                                            int source_x = x + src(X_COORD);
                                            int source_y = y + src(Y_COORD);
                                            int source_z = z + src(Y_COORD);
                                            int source_index = source_x * mult_x +
                                                               source_y * mult_y +
                                                               source_z * mult_z +
                                                               src_offset;
                                            int loc_x = x + dest(X_COORD) + region.x();
                                            int loc_y = y + dest(Y_COORD) + region.y();
                                            int loc_z = z + dest(Z_COORD) + region.z();
                                            int dest_x = loc_x - shift.x;
                                            int dest_y = loc_y - shift.y;
                                            int dest_z = loc_z - shift.z;
                                            const TV_INT destination_index(dest_x, dest_y, dest_z);  // NOLINT
                                            if (!(loc_x == cx1 || loc_x == cx2)) {
                                                (*fa)(dim, destination_index) = buffer[source_index];  // NOLINT
                                            }  // branch for averaging
                                        }  // read z for loop
                                    }  // read y for loop
                                }  // read x for loop
                                break;
                            }
                            case Y_COORD: {
                                const int_dimension_t ddx = dregion.dx();
                                const int_dimension_t ddy = 3;
                                const int_dimension_t ddz = dregion.dz();
                                assert(overlap(Y_COORD) == 3);
                                const int range_x = overlap(X_COORD);
                                const int range_y = 3 + 1;
                                const int range_z = overlap(Z_COORD);
                                const int mult_x = (ddy + 1) * ddz;
                                const int mult_y = ddz;
                                const int mult_z = 1;
                                const int src_offset = (ddx + 1) * ddy * ddz;  // offset into buffer
                                for (int x = 0; x < range_x; ++x) {
                                    for (int y = 0; y < range_y; ++y) {
                                        for (int z = 0; z < range_z; ++z) {
                                            int source_x = x + src(X_COORD);
                                            int source_y = y + src(Y_COORD);
                                            int source_z = z + src(Y_COORD);
                                            int source_index = source_x * mult_x +
                                                               source_y * mult_y +
                                                               source_z * mult_z +
                                                               src_offset;
                                            int loc_x = x + dest(X_COORD) + region.x();
                                            int loc_y = y + dest(Y_COORD) + region.y();
                                            int loc_z = z + dest(Z_COORD) + region.z();
                                            int dest_x = loc_x - shift.x;
                                            int dest_y = loc_y - shift.y;
                                            int dest_z = loc_z - shift.z;
                                            const TV_INT destination_index(dest_x, dest_y, dest_z);  // NOLINT
                                            if (!(loc_y == cy1 || loc_y == cy2)) {
                                                (*fa)(dim, destination_index) = buffer[source_index];  // NOLINT
                                            }  // branch for averaging
                                        }  // read z for loop
                                    }  // read y for loop
                                }  // read x for loop
                                break;
                            }
                            case Z_COORD: {
                                const int_dimension_t ddx = dregion.dx();
                                const int_dimension_t ddy = dregion.dy();
                                const int_dimension_t ddz = 3;
                                assert(overlap(Z_COORD) == 3);
                                const int range_x = overlap(X_COORD);
                                const int range_y = overlap(Y_COORD);
                                const int range_z = 3 + 1;
                                const int mult_x = ddy * (ddz + 1);
                                const int mult_y = ddz + 1;
                                const int mult_z = 1;
                                const int src_offset = (ddx + 1) * ddy * ddz +
                                                       ddx * (ddy + 1) * ddz;  // offset into buffer
                                for (int x = 0; x < range_x; ++x) {
                                    for (int y = 0; y < range_y; ++y) {
                                        for (int z = 0; z < range_z; ++z) {
                                            int source_x = x + src(X_COORD);
                                            int source_y = y + src(Y_COORD);
                                            int source_z = z + src(Y_COORD);
                                            int source_index = source_x * mult_x +
                                                               source_y * mult_y +
                                                               source_z * mult_z +
                                                               src_offset;
                                            int loc_x = x + dest(X_COORD) + region.x();
                                            int loc_y = y + dest(Y_COORD) + region.y();
                                            int loc_z = z + dest(Z_COORD) + region.z();
                                            int dest_x = loc_x - shift.x;
                                            int dest_y = loc_y - shift.y;
                                            int dest_z = loc_z - shift.z;
                                            const TV_INT destination_index(dest_x, dest_y, dest_z);  // NOLINT
                                            if (!(loc_z == cz1 || loc_z == cz2)) {
                                                (*fa)(dim, destination_index) = buffer[source_index];  // NOLINT
                                            }  // branch for averaging
                                        }  // read z for loop
                                    }  // read y for loop
                                }  // read x for loop
                                break;
                            }
                        }  // switch dim
                    }  // for dim
                }  //  if read anything
            }  //  outermost for loop over entired read set

            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Face Array (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
        }

        /** Take a FaceArray described by region and write it out to the
         *  PhysicalDataInstance objects in the objects array. */
        template<typename T> static void WriteFaceArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                typename PhysBAM::ARRAY<T, FaceIndex>* fa) {
            timer::StartTimer(timer::k2);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Face Array (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (write_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Write Face Array (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k2);
                return;
            }
            DataArray::const_iterator iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                const GeometricRegion dregion = data->region();
                const int_dimension_t ddx = dregion.dx();
                const int_dimension_t ddy = dregion.dy();
                const int_dimension_t ddz = dregion.dz();

                Dimension3Vector overlap = GetOverlapSize(dregion, region);
                if (!HasOverlap(overlap)) {continue;}

                T* buffer = reinterpret_cast<T*>(data->buffer());

                Dimension3Vector src = GetOffset(region, dregion);
                Dimension3Vector dest = GetOffset(dregion, region);

                //  x, y and z values are stored separately due to the
                // difference in number of x, y and z values in face arrays
                for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                    // setup for vx/ vy/ vz
                    int range_x, range_y, range_z;
                    int mult_x, mult_y;
                    int dst_offset;
                    switch (dim) {
                        case X_COORD:
                            range_x = overlap(X_COORD) + 1;
                            range_y = overlap(Y_COORD);
                            range_z = overlap(Z_COORD);
                            mult_x = ddy * ddz;
                            mult_y = ddz;
                            dst_offset = 0;  // offset into buffer
                            break;
                        case Y_COORD:
                            range_x = overlap(X_COORD);
                            range_y = overlap(Y_COORD) + 1;
                            range_z = overlap(Z_COORD);
                            mult_x = (ddy + 1) * ddz;
                            mult_y = ddz;
                            dst_offset = (ddx + 1) * ddy * ddz;  // offset into buffer
                            break;
                        case Z_COORD:
                            range_x = overlap(X_COORD);
                            range_y = overlap(Y_COORD);
                            range_z = overlap(Z_COORD) + 1;
                            mult_x = ddy * (ddz + 1);
                            mult_y = ddz + 1;
                            dst_offset = (ddx + 1) * ddy * ddz +
                                         ddx * (ddy + 1) * ddz;  // offset into buffer
                            break;
                    }
                    // printf("Region %s\n", dregion.ToNetworkData().c_str());
                    for (int x = 0; x < range_x; ++x) {
                        for (int y = 0; y < range_y; ++y) {
                            int dest_x = x + dest(X_COORD);
                            int dest_y = y + dest(Y_COORD);
                            int destination_index = dest_x * mult_x +
                                                  + dest_y * mult_y +
                                                  + dst_offset;
                            int source_x = x + src(X_COORD) + region.x() - shift.x;
                            int source_y = y + src(Y_COORD) + region.y() - shift.y;
                            int source_z = 0 + src(Z_COORD) + region.z() - shift.z;
                            TV_INT source_index(source_x, source_y, source_z);

                            // The PhysBAM FACE_ARRAY object abstracts away whether
                            // the data is stored in struct of array or array of struct
                            // form (in practice, usually struct of arrays
                            assert(destination_index + range_z  <= (int)data->size() / (int) sizeof(T) && destination_index >= 0);  // NOLINT
                            memcpy(&(buffer[destination_index]),
                                   &((*fa)(dim, source_index)),  // NOLINT
                                   range_z * sizeof(T));
                        }  // write for y
                    }  // write for x
                }  // for dim
            }  // write set for loop
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Face Array (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k2);
        }

        /*
         * Interface for deleting particles - delete particles corresponding to
         * the regions covered by da. May change DataArray to a region array
         * instead.
         */
        static void DeleteParticles(
                const Coord &shift,
                const std::vector<GeometricRegion> &regions,
                ParticleContainer *particle_container,
                const int_dimension_t scale,
                bool positive) {
            timer::StartTimer(timer::k10);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (regions.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Delete Particles (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k10);
                return;
            }

            ParticleArray *particles;
            if (positive) {
                particles = &particle_container->positive_particles;
            } else {
                particles = &particle_container->negative_particles;
            }

            Coord neg_shift(-shift.x, -shift.y, -shift.z);

            for (size_t r = 0; r < regions.size(); ++r) {
                GeometricRegion region = regions[r];
                GeometricRegion pregion = regions[r];
                pregion.Translate(neg_shift);

                for (int x = pregion.x(); x <= pregion.x() + pregion.dx(); ++x) {
                    for (int y = pregion.y(); y <= pregion.y() + pregion.dy() ; ++y) {
                        for (int z = pregion.z(); z <= pregion.z() + pregion.dz(); ++z) {
                            TV_INT bucket_index(x, y, z);

                            if (!(x == pregion.x() || x == pregion.x() + pregion.dx() ||
                                  y == pregion.y() || y == pregion.y() + pregion.dy() ||
                                  z == pregion.z() || z == pregion.z() + pregion.dz())) {
                                particle_container->Free_Particle_And_Clear_Pointer(
                                        (*particles)(bucket_index));
                            } else {
                                // TODO(later): can optimize this further
                                ParticleBucket *particle_bucket = (*particles)(bucket_index);
                                if (!particle_bucket)
                                    continue;

                                ParticleBucket *particle_bucket_history = particle_bucket;
                                (*particles)(bucket_index) = NULL;
                                bool new_particles_alloc = false;
                                ParticleBucket *particle_new_bucket = NULL;

                                while (particle_bucket) {
                                    for (int i = 1;
                                            i <= particle_bucket->array_collection->Size();
                                            ++i) {
                                        TV particle_position = particle_bucket->X(i);
                                        TV absolute_position = particle_position *
                                            static_cast<float>(scale) + 1.0;
                                        if (absolute_position.x < region.x() ||
                                            absolute_position.x >= region.x() + region.dx() ||
                                            absolute_position.y < region.y() ||
                                            absolute_position.y >= region.y() + region.dy() ||
                                            absolute_position.z < region.z() ||
                                            absolute_position.z >= region.z() + region.dz()) {
                                            if (!new_particles_alloc) {
                                                new_particles_alloc = true;
                                                (*particles)(bucket_index) =
                                                    particle_container->
                                                    Allocate_Particles(particle_container->
                                                            template_particles);
                                                particle_new_bucket = (*particles)(bucket_index);
                                            }

                                            int index = particle_container->
                                                Add_Particle(particle_new_bucket);
                                            particle_new_bucket->X(index) = particle_bucket->X(i);
                                            particle_new_bucket->radius(index) =
                                                particle_bucket->radius(i);
                                            particle_new_bucket->quantized_collision_distance(index) = // NOLINT
                                                particle_bucket->quantized_collision_distance(i);
                                            if (particle_container->store_unique_particle_id) {
                                                PhysBAM::ARRAY_VIEW<int> *new_id =
                                                    particle_new_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                PhysBAM::ARRAY_VIEW<int> *id =
                                                    particle_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                (*new_id)(index) = (*id)(i);
                                            }
                                        }
                                    }
                                    particle_bucket = particle_bucket->next;
                                }
                                // Delete the original bucket.
                                particle_container->
                                    Free_Particle_And_Clear_Pointer(particle_bucket_history);
                            }
                        }
                    }
                }
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k10);
        }

        static void DeleteRemovedParticles(
                const Coord &shift,
                const std::vector<GeometricRegion> &regions,
                ParticleContainer *particle_container,
                const int_dimension_t scale,
                bool positive) {
            timer::StartTimer(timer::k11);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Removed Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (regions.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Delete Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k11);
                return;
            }

            RemovedParticleArray *particles;
            if (positive) {
                particles = &particle_container->removed_positive_particles;
            } else {
                particles = &particle_container->removed_negative_particles;
            }

            Coord neg_shift(-shift.x, -shift.y, -shift.z);

            for (size_t r = 0; r < regions.size(); ++r) {
                GeometricRegion region = regions[r];
                GeometricRegion pregion = regions[r];
                pregion.Translate(neg_shift);

                for (int z = pregion.z(); z <= pregion.z() + pregion.dz(); ++z) {
                    for (int y = pregion.y(); y <= pregion.y() + pregion.dy() ; ++y) {
                        for (int x = pregion.x(); x <= pregion.x() + pregion.dx(); ++x) {
                            TV_INT bucket_index(x, y, z);

                            if (!(x == pregion.x() || x == pregion.x() + pregion.dx() ||
                                  y == pregion.y() || y == pregion.y() + pregion.dy() ||
                                  z == pregion.z() || z == pregion.z() + pregion.dz())) {
                                particle_container->Free_Particle_And_Clear_Pointer(
                                        (*particles)(bucket_index));
                            } else {
                                // TODO(later): can optimize this further
                                RemovedParticleBucket *particle_bucket = (*particles)(bucket_index); // NOLINT
                                if (!particle_bucket)
                                    continue;

                                RemovedParticleBucket *particle_bucket_history = particle_bucket;
                                (*particles)(bucket_index) = NULL;
                                bool new_particles_alloc = false;
                                RemovedParticleBucket *particle_new_bucket = NULL;

                                while (particle_bucket) {
                                    for (int i = 1;
                                            i <= particle_bucket->array_collection->Size();
                                            ++i) {
                                        TV particle_position = particle_bucket->X(i);
                                        TV absolute_position = particle_position *
                                            static_cast<float>(scale) + 1.0;
                                        if (absolute_position.x < region.x() ||
                                            absolute_position.x >= region.x() + region.dx() ||
                                            absolute_position.y < region.y() ||
                                            absolute_position.y >= region.y() + region.dy() ||
                                            absolute_position.z < region.z() ||
                                            absolute_position.z >= region.z() + region.dz()) {
                                            if (!new_particles_alloc) {
                                                new_particles_alloc = true;
                                                (*particles)(bucket_index) =
                                                    particle_container->
                                                    Allocate_Particles(particle_container->
                                                            template_removed_particles);
                                                particle_new_bucket = (*particles)(bucket_index);
                                            }

                                            int index = particle_container->
                                                Add_Particle(particle_new_bucket);
                                            particle_new_bucket->X(index) = particle_bucket->X(i);
                                            particle_new_bucket->radius(index) =
                                                particle_bucket->radius(i);
                                            particle_new_bucket->quantized_collision_distance(index) = // NOLINT
                                                particle_bucket->quantized_collision_distance(i);
                                            if (particle_container->store_unique_particle_id) {
                                                PhysBAM::ARRAY_VIEW<int> *new_id =
                                                    particle_new_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                PhysBAM::ARRAY_VIEW<int> *id =
                                                    particle_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                (*new_id)(index) = (*id)(i);
                                            }
                                            particle_new_bucket->V(index) = particle_bucket->V(i);
                                        }
                                    }
                                    particle_bucket = particle_bucket->next;
                                }
                                // Delete the original bucket.
                                particle_container->
                                    Free_Particle_And_Clear_Pointer(particle_bucket_history);
                            }
                        }
                    }
                }
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k11);
        }

        /*
         * Interface for deleting particles - delete particles outside inner,
         * inside outer.
         */
        static void DeleteParticles(
                const Coord &shift,
                GeometricRegion &inner,
                GeometricRegion &outer,
                ParticleContainer *particle_container,
                const int_dimension_t scale,
                bool positive) {
            timer::StartTimer(timer::k12);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            ParticleArray *particles;
            if (positive) {
                particles = &particle_container->positive_particles;
            } else {
                particles = &particle_container->negative_particles;
            }
            Coord neg_shift(-shift.x, -shift.y, -shift.z);
            GeometricRegion pouter = outer;
            pouter.Translate(neg_shift);
            GeometricRegion pinner = inner;
            pinner.Translate(neg_shift);

            const int_dimension_t pix = pinner.x();
            const int_dimension_t piy = pinner.y();
            const int_dimension_t piz = pinner.z();
            const int_dimension_t pidx = pinner.dx();
            const int_dimension_t pidy = pinner.dy();
            const int_dimension_t pidz = pinner.dz();
            const int_dimension_t pox = pouter.x();
            const int_dimension_t poy = pouter.y();
            const int_dimension_t poz = pouter.z();
            const int_dimension_t podx = pouter.dx();
            const int_dimension_t pody = pouter.dy();
            const int_dimension_t podz = pouter.dz();

            const int_dimension_t ix = inner.x();
            const int_dimension_t iy = inner.y();
            const int_dimension_t iz = inner.z();
            const int_dimension_t idx = inner.dx();
            const int_dimension_t idy = inner.dy();
            const int_dimension_t idz = inner.dz();

            for (int z = poz; z <= poz + podz; ++z) {
                for (int y = poy; y <= poy + pody ; ++y) {
                    for (int x = pox; x <= pox + podx; ++x) {
                        TV_INT bucket_index(x, y, z);

                        if ((x < pix || x > pix + pidx ||
                             y < piy || y > piy + pidy ||
                             z < piz || z > piz + pidz)) {
                            particle_container->Free_Particle_And_Clear_Pointer(
                                    (*particles)(bucket_index));
                        } else {
                            // TODO(later): can optimize this further
                            if (x == pix || x == pix + pidx ||
                                y == piy || y == piy + pidy ||
                                z == piz || z == piz + pidz) {
                                ParticleBucket *particle_bucket = (*particles)(bucket_index);
                                if (!particle_bucket)
                                    continue;

                                ParticleBucket *particle_bucket_history = particle_bucket;
                                (*particles)(bucket_index) = NULL;
                                bool new_particles_alloc = false;
                                ParticleBucket *particle_new_bucket = NULL;

                                while (particle_bucket) {
                                    for (int i = 1;
                                            i <= particle_bucket->array_collection->Size();
                                            ++i) {
                                        TV particle_position = particle_bucket->X(i);
                                        TV absolute_position = particle_position *
                                            static_cast<float>(scale) + 1.0;
                                        if (absolute_position.x >= ix &&
                                            absolute_position.x < ix + idx &&
                                            absolute_position.y >= iy &&
                                            absolute_position.y < iy + idy &&
                                            absolute_position.z >= iz &&
                                            absolute_position.z < iz + idz) {
                                            if (!new_particles_alloc) {
                                                new_particles_alloc = true;
                                                (*particles)(bucket_index) =
                                                    particle_container->
                                                    Allocate_Particles(particle_container->
                                                            template_particles);
                                                particle_new_bucket = (*particles)(bucket_index);
                                            }

                                            int index = particle_container->
                                                Add_Particle(particle_new_bucket);
                                            particle_new_bucket->X(index) = particle_bucket->X(i);
                                            particle_new_bucket->radius(index) =
                                                particle_bucket->radius(i);
                                            particle_new_bucket->quantized_collision_distance(index) = // NOLINT
                                                particle_bucket->quantized_collision_distance(i);
                                            if (particle_container->store_unique_particle_id) {
                                                PhysBAM::ARRAY_VIEW<int> *new_id =
                                                    particle_new_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                PhysBAM::ARRAY_VIEW<int> *id =
                                                    particle_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                (*new_id)(index) = (*id)(i);
                                            }
                                        }
                                    }
                                    particle_bucket = particle_bucket->next;
                                }
                                // Delete the original bucket.
                                particle_container->
                                    Free_Particle_And_Clear_Pointer(particle_bucket_history);
                            }
                        }
                    }
                }
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k12);
        }

        /*
         * Interface for deleting removed particles - delete particles outside inner,
         * inside outer.
         */
        static void DeleteRemovedParticles(
                const Coord &shift,
                GeometricRegion &inner,
                GeometricRegion &outer,
                ParticleContainer *particle_container,
                const int_dimension_t scale,
                bool positive) {
            timer::StartTimer(timer::k13);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Removed Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            RemovedParticleArray *particles;
            if (positive) {
                particles = &particle_container->removed_positive_particles;
            } else {
                particles = &particle_container->removed_negative_particles;
            }
            Coord neg_shift(-shift.x, -shift.y, -shift.z);
            GeometricRegion pouter = outer;
            pouter.Translate(neg_shift);
            GeometricRegion pinner = inner;
            pinner.Translate(neg_shift);

            for (int z = pouter.z(); z <= pouter.z() + pouter.dz(); ++z) {
                for (int y = pouter.y(); y <= pouter.y() + pouter.dy() ; ++y) {
                    for (int x = pouter.x(); x <= pouter.x() + pouter.dx(); ++x) {
                        TV_INT bucket_index(x, y, z);

                        if ((x < pinner.x() || x > pinner.x() + pinner.dx() ||
                             y < pinner.y() || y > pinner.y() + pinner.dy() ||
                             z < pinner.z() || z > pinner.z() + pinner.dz())) {
                            particle_container->Free_Particle_And_Clear_Pointer(
                                    (*particles)(bucket_index));
                        } else {
                            // TODO(later): can optimize this further
                            if (x == pinner.x() || x == pinner.x() + pinner.dx() ||
                                y == pinner.y() || y == pinner.y() + pinner.dy() ||
                                z == pinner.z() || z == pinner.z() + pinner.dz()) {
                                RemovedParticleBucket *particle_bucket = (*particles)(bucket_index);
                                if (!particle_bucket)
                                    continue;

                                RemovedParticleBucket *particle_bucket_history = particle_bucket;
                                (*particles)(bucket_index) = NULL;
                                bool new_particles_alloc = false;
                                RemovedParticleBucket *particle_new_bucket = NULL;

                                while (particle_bucket) {
                                    for (int i = 1;
                                            i <= particle_bucket->array_collection->Size();
                                            ++i) {
                                        TV particle_position = particle_bucket->X(i);
                                        TV absolute_position = particle_position *
                                            static_cast<float>(scale) + 1.0;
                                        if (absolute_position.x >= inner.x() &&
                                            absolute_position.x < inner.x() + inner.dx() &&
                                            absolute_position.y >= inner.y() &&
                                            absolute_position.y < inner.y() + inner.dy() &&
                                            absolute_position.z >= inner.z() &&
                                            absolute_position.z < inner.z() + inner.dz()) {
                                            if (!new_particles_alloc) {
                                                new_particles_alloc = true;
                                                (*particles)(bucket_index) =
                                                    particle_container->
                                                    Allocate_Particles(particle_container->
                                                            template_removed_particles);
                                                particle_new_bucket = (*particles)(bucket_index);
                                            }

                                            int index = particle_container->
                                                Add_Particle(particle_new_bucket);
                                            particle_new_bucket->X(index) = particle_bucket->X(i);
                                            particle_new_bucket->radius(index) =
                                                particle_bucket->radius(i);
                                            particle_new_bucket->quantized_collision_distance(index) = // NOLINT
                                                particle_bucket->quantized_collision_distance(i);
                                            if (particle_container->store_unique_particle_id) {
                                                PhysBAM::ARRAY_VIEW<int> *new_id =
                                                    particle_new_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                PhysBAM::ARRAY_VIEW<int> *id =
                                                    particle_bucket->array_collection->
                                                    template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID); // NOLINT
                                                (*new_id)(index) = (*id)(i);
                                            }
                                            particle_new_bucket->V(index) = particle_bucket->V(i);
                                        }
                                    }
                                    particle_bucket = particle_bucket->next;
                                }
                                // Delete the original bucket.
                                particle_container->
                                    Free_Particle_And_Clear_Pointer(particle_bucket_history);
                            }
                        }
                    }
                }
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Delete Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k13);
        }

        /* Reads the particle data from the PhysicalDataInstances "instances",
         * limited by the corresponding global region calculated from shifting the
         * local region "region" by the offset "shift", into the local region
         * "region" of the PhysBAM particle container "particle_container->.
         *
         * "positive" option specifies whether to work on positive particles or
         * negative particles.
         * "merge" option specifies whether to keep original particle data in
         * "particle_container->.
         */
        static void ReadParticles(const GeometricRegion &region,
                const Coord &shift,
                const DataArray &read_set,
                ParticleContainer *particle_container,
                const int_dimension_t kScale,
                bool positive,
                bool merge = false) {
            timer::StartTimer(timer::k3);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            ParticleArray* particles;
            if (positive) {
                particles = &particle_container->positive_particles;
            } else {
                particles = &particle_container->negative_particles;
            }

            Coord neg_shift(-shift.x, -shift.y, -shift.z);
            GeometricRegion pregion = region;
            pregion.Translate(neg_shift);

            // Checks whether the geometric region in the particle array is valid, and
            // clears corresponding buckets inside the geometric region if necessary.
            if (!merge) {
                for (int z = pregion.z(); z <= pregion.z() + pregion.dz(); ++z) {
                    for (int y = pregion.y(); y <= pregion.y() + pregion.dy(); ++y) {
                        for (int x = pregion.x(); x <= pregion.x() + pregion.dx(); ++x) {
                            TV_INT bucket_index(x, y, z);
                            particle_container->Free_Particle_And_Clear_Pointer(
                                    (*particles)(bucket_index));
                        }
                    }
                }
            }

            DataArray::const_iterator iter = read_set.begin();
            for (; iter != read_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                const GeometricRegion dregion = data->region();
                const size_t dsize = data->size();
                ParticleInternal* buffer =
                    reinterpret_cast<ParticleInternal*>(data->buffer());
                ParticleInternal* buffer_end = buffer + static_cast<int>(dsize)
                    / static_cast<int>(sizeof(ParticleInternal));

                for (ParticleInternal* p = buffer; p != buffer_end; ++p) {
                    TV absolute_position;
                    absolute_position.x = p->position[0];
                    absolute_position.y = p->position[1];
                    absolute_position.z = p->position[2];
                    // TODO(quhang) Needs to deal with the particles that lies exactly on
                    // the boundary.
                    if (absolute_position.x >= region.x() &&
                        absolute_position.x < region.x() + region.dx() &&
                        absolute_position.y >= region.y() &&
                        absolute_position.y < region.y() + region.dy() &&
                        absolute_position.z >= region.z() &&
                        absolute_position.z < region.z() + region.dz()) {
                        TV_INT bucket_index(round(absolute_position.x - shift.x),
                                round(absolute_position.y - shift.y),
                                round(absolute_position.z - shift.z));
                        // assert(particles->Valid_Index(bucket_index));
                        if (!(*particles)(bucket_index)) {
                            (*particles)(bucket_index) = particle_container->
                                Allocate_Particles(particle_container->template_particles);
                        }
                        ParticleBucket* particle_bucket =
                            (*particles)(bucket_index);

                        // Note that Add_Particle traverses a linked list of particle
                        // buckets, so it's O(N^2) time. Blech.
                        int index = particle_container->Add_Particle(particle_bucket);
                        particle_bucket->X(index) =
                            (absolute_position - 1.0) / (float) kScale; // NOLINT
                        particle_bucket->radius(index) = p->radius;
                        particle_bucket->quantized_collision_distance(index) =
                            p->quantized_collision_distance;
                        if (particle_container->store_unique_particle_id) {
                            PhysBAM::ARRAY_VIEW<int>* id = particle_bucket->array_collection->
                                template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                            (*id)(index) = p->id;
                        }
                    }
                }  // End the loop for buffer.
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k3);
        }


        /* Writes the particle data from PhysBAM particle container
         * "particle_container->, limited the local region "region", into the
         * corresponding global region of physical instances "instances" after
         * performing the coordinate shifting specified by "shift".
         *
         * "positive" option specifies whether to work on positive particles or
         * negative particles.
         */
        static void WriteParticles(const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                ParticleContainer *particle_container,
                const int_dimension_t kScale,
                bool positive) {
            timer::StartTimer(timer::k4);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (write_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Write Particles (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k4);
                return;
            }

            DataArray::const_iterator iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                data->ClearTempBuffer();
            }

            ParticleArray* particles;
            if (positive) {
                particles = &particle_container->positive_particles;
            } else {
                particles = &particle_container->negative_particles;
            }

            const Coord neg_shift(-shift.x, -shift.y, -shift.z);

            // Iterate across instances
            iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                const GeometricRegion data_region = data->region();
                GeometricRegion pregion = data_region;
                pregion.Translate(neg_shift);

                // Loop through each particle bucket in the specified region.
                for (int z = pregion.z();
                     z <= pregion.z() + pregion.dz(); ++z) {
                    for (int y = pregion.y();
                         y <= pregion.y() + pregion.dy(); ++y) {
                        for (int x = pregion.x();
                             x <= pregion.x() + pregion.dx(); ++x) {
                            TV_INT bucket_index(x, y, z);
                            ParticleBucket* particle_bucket = (*particles)(bucket_index);
                            while (particle_bucket != NULL) {
                                for (int i = 1;
                                     i <= particle_bucket->array_collection->Size();
                                     ++i) {
                                    TV particle_position = particle_bucket->X(i);
                                    TV absolute_position =
                                        particle_position * (float) kScale + 1.0; // NOLINT
                                    // TODO(quhang) Needs to deal with the case when the particle
                                    // lies exactly on the boundary.
                                    // If it's inside the region of the physical data instance.
                                    if (absolute_position.x >=
                                        data_region.x() &&
                                        absolute_position.x <
                                        (data_region.x() + data_region.dx()) &&
                                        absolute_position.y >=
                                        data_region.y() &&
                                        absolute_position.y <
                                        (data_region.y() + data_region.dy()) &&
                                        absolute_position.z >=
                                        data_region.z() &&
                                        absolute_position.z <
                                        (data_region.z() + data_region.dz())) {
                                        ParticleInternal particle_buffer;
                                        particle_buffer.position[0] = absolute_position.x;
                                        particle_buffer.position[1] = absolute_position.y;
                                        particle_buffer.position[2] = absolute_position.z;
                                        particle_buffer.radius = particle_bucket->radius(i);
                                        particle_buffer.quantized_collision_distance =
                                            particle_bucket->quantized_collision_distance(i);
                                        if (particle_container->store_unique_particle_id) {
                                            PhysBAM::ARRAY_VIEW<int>* id =
                                                particle_bucket->array_collection->
                                                template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                                            particle_buffer.id = (*id)(i);
                                        }
                                        data->AddToTempBuffer(
                                                reinterpret_cast<char*>(&particle_buffer),
                                                sizeof(particle_buffer));
                                    }
                                }  // Finish looping through all particles.
                                particle_bucket = particle_bucket->next;
                            }
                        }
                    }
                }
                // commit the result
                data->CommitTempBuffer();
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k4);
        }

        /* Reads the removed particle data from the PhysicalDataInstances
         * "instances", limited by the corresponding global region calculated from
         * shifting the local region "region" by the offset "shift", into the local
         * region "region" of the PhysBAM particle container "particle_container->.
         *
         * "positive" option specifies whether to work on positive particles or
         * negative particles.
         * "merge" option specifies whether to keep original particle data in
         * "particle_container->.
         */
        static void ReadRemovedParticles(const GeometricRegion &region,
                const Coord &shift,
                const DataArray &read_set,
                ParticleContainer *particle_container,
                const int_dimension_t kScale,
                bool positive,
                bool merge = false) {
            timer::StartTimer(timer::k14);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Removed Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            RemovedParticleArray* particles;
            if (positive) {
                particles = &particle_container->removed_positive_particles;
            } else {
                particles = &particle_container->removed_negative_particles;
            }

            Coord neg_shift(-shift.x, -shift.y, -shift.z);
            GeometricRegion pregion = region;
            pregion.Translate(neg_shift);

            // Checks whether the geometric region in the particle array is valid, and
            // clears corresponding buckets inside the geometric region if necessary.
            if (!merge) {
                for (int z = pregion.z(); z <= pregion.z() + pregion.dz(); ++z) {
                    for (int y = pregion.y(); y <= pregion.y() + pregion.dy(); ++y) {
                        for (int x = pregion.x(); x <= pregion.x() + pregion.dx(); ++x) {
                            TV_INT bucket_index(x, y, z);
                            if ((*particles)(bucket_index)) {
                                delete (*particles)(bucket_index);
                                (*particles)(bucket_index) = NULL;
                            }
                        }
                    }
                }
            }

            if (read_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Read Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                dbg(DBG_WARN, "Physical data instances are empty.\n");
                timer::StopTimer(timer::k14);
                return;
            }

            DataArray::const_iterator iter = read_set.begin();
            for (; iter != read_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                const GeometricRegion dregion = data->region();
                const size_t dsize = data->size();
                RemovedParticleInternal* buffer =
                    reinterpret_cast<RemovedParticleInternal*>(data->buffer());
                RemovedParticleInternal* buffer_end = buffer
                    + static_cast<int>(dsize)
                    / static_cast<int>(sizeof(RemovedParticleInternal));

                for (RemovedParticleInternal* p = buffer; p != buffer_end; ++p) {
                    TV absolute_position;
                    absolute_position.x = p->position[0];
                    absolute_position.y = p->position[1];
                    absolute_position.z = p->position[2];
                    if (absolute_position.x >= region.x() &&
                        absolute_position.x < region.x() + region.dx() &&
                        absolute_position.y >= region.y() &&
                        absolute_position.y < region.y() + region.dy() &&
                        absolute_position.z >= region.z() &&
                        absolute_position.z < region.z() + region.dz()) {
                        TV_INT bucket_index(round(absolute_position.x - shift.x),
                                round(absolute_position.y - shift.y),
                                round(absolute_position.z - shift.z));
                        // assert(particles->Valid_Index(bucket_index));
                        // NOTE(By Chinmayee): Please comment out these changes and don't
                        // delete them when pushing any updates, till we verify that the
                        // code works correctly, and does not give any assertion failure or
                        // seg fault on both Linux and Mac.
                        if (!(*particles)(bucket_index)) {
                            (*particles)(bucket_index) =
                                particle_container->Allocate_Particles(
                                        particle_container->template_removed_particles);
                        }
                        RemovedParticleBucket* particle_bucket =
                            (*particles)(bucket_index);

                        // Note that Add_Particle traverses a linked list of particle
                        // buckets, so it's O(N^2) time. Blech.
                        int index = particle_container->Add_Particle(particle_bucket);
                        particle_bucket->X(index) =
                            (absolute_position - 1.0) / (float) kScale; // NOLINT
                        particle_bucket->radius(index) = p->radius;
                        particle_bucket->quantized_collision_distance(index) =
                            p->quantized_collision_distance;
                        if (particle_container->store_unique_particle_id) {
                            PhysBAM::ARRAY_VIEW<int>* id = particle_bucket->array_collection->
                                template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                            (*id)(index) = p->id;
                        }
                        particle_bucket->V(index) = TV(p->v[0],
                                p->v[1],
                                p->v[2]);
                    }
                }  // End the loop for buffer.
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k14);
        }

        /* Writes the removed particle data from PhysBAM particle container
         * "particle_container->, limited the local region "region", into the
         * corresponding global region of physical instances "instances" after
         * performing the coordinate shifting specified by "shift".
         *
         * "positive" option specifies whether to work on positive particles or
         * negative particles.
         */
        // TODO(quhang) The similar optimization in WriteParticles can be put here.
        static void WriteRemovedParticles(const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                ParticleContainer *particle_container,
                const int_dimension_t kScale,
                bool positive
                ) {
            timer::StartTimer(timer::k15);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Removed Particles (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (write_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Write Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k15);
                return;
            }
            DataArray::const_iterator iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                data->ClearTempBuffer();
            }

            RemovedParticleArray* particles;
            if (positive) {
                particles = &particle_container->removed_positive_particles;
            } else {
                particles = &particle_container->removed_negative_particles;
            }

            Coord neg_shift(-shift.x, -shift.y, -shift.z);

            // Iterate across instances
            iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                const GeometricRegion data_region = data->region();
                GeometricRegion pregion = data_region;
                pregion.Translate(neg_shift);

                // Loop through each particle bucket in the specified region.
                for (int z = pregion.z();
                     z <= pregion.z() + pregion.dz(); ++z) {
                    for (int y = pregion.y();
                         y <= pregion.y() + pregion.dy(); ++y) {
                        for (int x = pregion.x();
                             x <= pregion.x() + pregion.dx(); ++x) {
                            TV_INT bucket_index(x, y, z);
                            RemovedParticleBucket* particle_bucket = (*particles)(bucket_index);
                            while (particle_bucket != NULL) {
                                for (int i = 1;
                                    i <= particle_bucket->array_collection->Size();
                                    ++i) {
                                    TV particle_position = particle_bucket->X(i);
                                    TV absolute_position =
                                        particle_position * (float) kScale + 1.0; // NOLINT
                                    // If it's inside the region of the physical data instance.
                                    if (absolute_position.x >=
                                        data_region.x() &&
                                        absolute_position.x <
                                        (data_region.x() + data_region.dx()) &&
                                        absolute_position.y >=
                                        data_region.y() &&
                                        absolute_position.y <
                                        (data_region.y() + data_region.dy()) &&
                                        absolute_position.z >=
                                        data_region.z() &&
                                        absolute_position.z <
                                        (data_region.z() + data_region.dz())) {
                                        RemovedParticleInternal particle_buffer;
                                        particle_buffer.position[0] = absolute_position.x;
                                        particle_buffer.position[1] = absolute_position.y;
                                        particle_buffer.position[2] = absolute_position.z;
                                        particle_buffer.radius = particle_bucket->radius(i);
                                        particle_buffer.quantized_collision_distance =
                                            particle_bucket->quantized_collision_distance(i);
                                        if (particle_container->store_unique_particle_id) {
                                            PhysBAM::ARRAY_VIEW<int>* id =
                                                particle_bucket->array_collection->
                                                template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                                            particle_buffer.id = (*id)(i);
                                        }
                                        particle_buffer.v[0] = particle_bucket->V(i).x;
                                        particle_buffer.v[1] = particle_bucket->V(i).y;
                                        particle_buffer.v[2] = particle_bucket->V(i).z;
                                        data->AddToTempBuffer(
                                                reinterpret_cast<char*>(&particle_buffer),
                                                sizeof(particle_buffer));
                                    }
                                }  // Finish looping through all particles.
                                particle_bucket = particle_bucket->next;
                            }
                        }
                    }
                }
                // commit the result
                data->CommitTempBuffer();
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Removed Particles (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
                timer::StopTimer(timer::k15);
        }

        /* Read scalar array from PhysicalDataInstances specified by instances,
         * limited by the GeometricRegion specified by region, into the
         * ScalarArray specified by dest. This allocates a new scalar array. */
        template<typename T> static void ReadScalarArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &read_set,
                typename PhysBAM::ARRAY<T, TV_INT>* sa) {
            timer::StartTimer(timer::k5);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Scalar Array (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (read_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Read Scalar Array (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k5);
                return;
            }
            sa->hash_code = 0;
            for (size_t i = 0; i < read_set.size(); ++i) {
                PhysBAMData* data = static_cast<PhysBAMData*>(read_set[i]);
                const GeometricRegion dregion = data->region();
                Dimension3Vector overlap = GetOverlapSize(dregion, region);

                if (HasOverlap(overlap)) {
                    T* buffer  = reinterpret_cast<T*>(data->buffer());

                    Dimension3Vector src  = GetOffset(dregion, region);
                    Dimension3Vector dest = GetOffset(region, dregion);

                    for (int x = 0; x < overlap(X_COORD); ++x) {
                        for (int y = 0; y < overlap(Y_COORD); ++y) {
                            int source_x = x + src(X_COORD);
                            int source_y = y + src(Y_COORD);
                            int source_z = 0 + src(Z_COORD);
                            int source_index =
                                (source_x * (dregion.dy() * dregion.dz())) +
                                (source_y * (dregion.dz())) +
                                source_z;
                            int dest_x = x + dest(X_COORD) + region.x() - shift.x;
                            int dest_y = y + dest(Y_COORD) + region.y() - shift.y;
                            int dest_z = 0 + dest(Z_COORD) + region.z() - shift.z;
                            TV_INT destination_index(dest_x, dest_y, dest_z);
                            memcpy(&((*sa)(destination_index)),  // NOLINT
                                   &(buffer[source_index]),
                                   sizeof(T) * overlap(Z_COORD));
                        }
                    }
                }
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Read Scalar Array (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k5);
        }

        /* Write scalar array data into PhysicalDataInstances specified by instances,
         * limited by the GeometricRegion region. This frees the physbam scalar array. */
        template<typename T> static void WriteScalarArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                typename PhysBAM::ARRAY<T, TV_INT>* sa) {
            timer::StartTimer(timer::k6);
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Scalar Array (New Translator) start : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            if (write_set.empty()) {
                if (log) {
                    std::stringstream msg;
                    pid_t tid = syscall(SYS_gettid);
                    msg << "### TID: " << tid << "  Write Scalar Array (New Translator) end : " << log->GetTime(); // NOLINT
                    log->WriteToFile(msg.str());
                }
                timer::StopTimer(timer::k6);
                return;
            }
            for (size_t i = 0; i < write_set.size(); ++i) {
                PhysBAMData* data = static_cast<PhysBAMData*>(write_set[i]);
                const GeometricRegion dregion = data->region();
                Dimension3Vector overlap = GetOverlapSize(dregion, region);

                if (HasOverlap(overlap)) {
                    T* buffer  = reinterpret_cast<T*>(data->buffer());

                    Dimension3Vector src  = GetOffset(region, dregion);
                    Dimension3Vector dest = GetOffset(dregion, region);

                    for (int x = 0; x < overlap(X_COORD); ++x) {
                        for (int y = 0; y < overlap(Y_COORD); ++y) {
                            int dest_x = x + dest(X_COORD);
                            int dest_y = y + dest(Y_COORD);
                            int dest_z = 0 + dest(Z_COORD);
                            int destination_index =
                                (dest_x * (dregion.dy() * dregion.dz())) +
                                (dest_y * (dregion.dz())) +
                                dest_z;
                            int source_x = x + src(X_COORD) + region.x() - shift.x;
                            int source_y = y + src(Y_COORD) + region.y() - shift.y;
                            int source_z = 0 + src(Z_COORD) + region.z() - shift.z;
                            TV_INT source_index(source_x, source_y, source_z);
                            // assert(destination_index < dsize / (int) sizeof(T) && destination_index >= 0); // NOLINT
                            memcpy(&(buffer[destination_index]),
                                   &((*sa)(source_index)),  // NOLINT
                                   overlap(Z_COORD) * sizeof(T));
                        }
                    }
                }
            }
            if (log) {
                std::stringstream msg;
                pid_t tid = syscall(SYS_gettid);
                msg << "### TID: " << tid << "  Write Scalar Array (New Translator) end : " << log->GetTime(); // NOLINT
                log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k6);
        }

        template<typename T> static void ReadCompressedScalarArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &read_set,
                PhysBAM::VECTOR_ND<T>* data,
                const int_dimension_t data_length,
                const PhysBAM::ARRAY<int, TV_INT>& index_data) {
          timer::StartTimer(timer::k7);
          if (data->Size() != data_length) {
            data->Resize(data_length);
          }
          if (log) {
            std::stringstream msg;
            msg << "### Read Compressed Scalar Array (New Translator) start : "
                << log->GetTime();
            log->WriteToFile(msg.str());
          }
          if (read_set.empty()) {
            if (log) {
              std::stringstream msg;
              msg << "### Read Compressed Scalar Array (New Translator) end : "
                  << log->GetTime();
              log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k7);
            return;
          }
          for (size_t i = 0; i < read_set.size(); ++i) {
            PhysBAMDataWithMeta* nimbus_data =
                dynamic_cast<PhysBAMDataWithMeta*>(read_set[i]);  // NOLINT
            const GeometricRegion nimbus_dregion = nimbus_data->region();
            // assert(nimbus_data != NULL);
            GeometricRegion inter_region = GeometricRegion::GetIntersection(
                nimbus_dregion, region);
            char* buffer = nimbus_data->buffer();
            if (buffer == NULL || nimbus_data->size() == 0) {
              continue;
            }
            // assert(nimbus_data->has_meta_data());
            int_dimension_t elements =
                nimbus_data->meta_data_size() / 3 / sizeof(int_dimension_t);
            char* real_data_buffer = buffer + nimbus_data->meta_data_size();
            for (int_dimension_t rank = 0; rank < elements; ++rank) {
              int_dimension_t x = *reinterpret_cast<int_dimension_t*>(buffer);
              buffer += sizeof(int_dimension_t);
              int_dimension_t y = *reinterpret_cast<int_dimension_t*>(buffer);
              buffer += sizeof(int_dimension_t);
              int_dimension_t z = *reinterpret_cast<int_dimension_t*>(buffer);
              buffer += sizeof(int_dimension_t);
              T value = *reinterpret_cast<T*>(real_data_buffer);
              real_data_buffer += sizeof(T);
              if (x >= inter_region.x() &&
                  x < inter_region.x() + inter_region.dx() &&
                  y >= inter_region.y() &&
                  y < inter_region.y() + inter_region.dy() &&
                  z >= inter_region.z() &&
                  z < inter_region.z() + inter_region.dz()) {
                int m_index = index_data(x - shift.x, y - shift.y, z - shift.z);
                // assert(m_index >= 1);
                // assert(m_index <= data_length);
                (*data)(m_index) = value;
              }
            }
          }
          if (log) {
            std::stringstream msg;
            msg << "### Read Compressed Scalar Array (New Translator) end : "
                << log->GetTime();
            log->WriteToFile(msg.str());
          }
          timer::StopTimer(timer::k7);
        }

        template<typename T> static void WriteCompressedScalarArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                const PhysBAM::VECTOR_ND<T>& data,
                const int_dimension_t data_length,
                const PhysBAM::ARRAY<int, TV_INT>& index_data) {
          timer::StartTimer(timer::k8);
          if (log) {
            std::stringstream msg;
            msg << "### Write Compressed Scalar Array (New Translator) start : "
                << log->GetTime();
            log->WriteToFile(msg.str());
          }
          if (write_set.empty()) {
            if (log) {
              std::stringstream msg;
              msg << "### Write Compressed Scalar Array (New Translator) end : "
                  << log->GetTime();
              log->WriteToFile(msg.str());
            }
            timer::StopTimer(timer::k8);
            return;
          }
          for (size_t i = 0; i < write_set.size(); ++i) {
            PhysBAMDataWithMeta* nimbus_data =
                dynamic_cast<PhysBAMDataWithMeta*>(write_set[i]);  // NOLINT
            // assert(nimbus_data != NULL);
            const GeometricRegion nimbus_dregion = nimbus_data->region();
            GeometricRegion inter_region = GeometricRegion::GetIntersection(
                nimbus_dregion, region);
            if (!inter_region.NoneZeroArea()) {
              continue;
            }
            nimbus_data->ClearTempBuffer();
            std::list<T> buffer;
            for (int_dimension_t x = inter_region.x();
                 x < inter_region.x() + inter_region.dx();
                 ++x)
              for (int_dimension_t y = inter_region.y();
                   y < inter_region.y() + inter_region.dy();
                   ++y)
                for (int_dimension_t z = inter_region.z();
                     z < inter_region.z() + inter_region.dz();
                     ++z) {
                  int m_index = index_data(TV_INT(
                          x - shift.x, y - shift.y, z - shift.z));
                  if (m_index >= 1) {
                    nimbus_data->AddToTempBuffer(
                        reinterpret_cast<char*>(&x),
                        sizeof(x));
                    nimbus_data->AddToTempBuffer(
                        reinterpret_cast<char*>(&y),
                        sizeof(y));
                    nimbus_data->AddToTempBuffer(
                        reinterpret_cast<char*>(&z),
                        sizeof(z));
                    buffer.push_back(data(m_index));
                  }
                }
            nimbus_data->MarkMetaDataInTempBuffer();
            if (!buffer.empty()) {
              for (typename std::list<T>::iterator iter = buffer.begin();
                   iter != buffer.end();
                   ++iter) {
                T value = *iter;
                nimbus_data->AddToTempBuffer(
                    reinterpret_cast<char*>(&value), sizeof(value));
              }
            }
            nimbus_data->CommitTempBuffer();
          }
          if (log) {
            std::stringstream msg;
            msg << "### Write Compressed Scalar Array (New Translator) end : "
                << log->GetTime();
            log->WriteToFile(msg.str());
          }
          timer::StopTimer(timer::k8);
        }

    private:
        /* Return a vector describing what the offset of b
           within a, such that a.x + offset = b.x. If
           offset is negative, return 0. */
        static Dimension3Vector GetOffset(
                const GeometricRegion &a,
                const GeometricRegion &b) {
            Dimension3Vector result;

            // If source is > than dest, its offset is zero (it's contained),
            // otherwise the offset is the difference between the values.
            int_dimension_t x = b.x() - a.x();
            int_dimension_t y = b.y() - a.y();
            int_dimension_t z = b.z() - a.z();
            result(X_COORD) = (x >= 0)? x:0;
            result(Y_COORD) = (y >= 0)? y:0;
            result(Z_COORD) = (z >= 0)? z:0;

            return result;
        }

        static Dimension3Vector GetOverlapSize(
                const GeometricRegion &src,
                const GeometricRegion &dest) {
            Dimension3Vector result;

            int_dimension_t x_start = std::max(src.x(), dest.x());
            int_dimension_t x_end   = std::min(src.x() + src.dx(),
                    dest.x() + dest.dx());
            int_dimension_t x_size = x_end - x_start;

            int_dimension_t y_start = std::max(src.y(), dest.y());
            int_dimension_t y_end   = std::min(src.y() + src.dy(),
                    dest.y() + dest.dy());
            int_dimension_t y_size = y_end - y_start;

            int_dimension_t z_start = std::max(src.z(), dest.z());
            int_dimension_t z_end   = std::min(src.z() + src.dz(),
                    dest.z() + dest.dz());
            int_dimension_t z_size = z_end - z_start;

            result(X_COORD) = (x_size >= 0)? x_size:0;
            result(Y_COORD) = (y_size >= 0)? y_size:0;
            result(Z_COORD) = (z_size >= 0)? z_size:0;

            return result;
        }

        static bool HasOverlap(Dimension3Vector overlapSize) {
            return (overlapSize(X_COORD) > 0 &&
                    overlapSize(Y_COORD) > 0 &&
                    overlapSize(Z_COORD) > 0);
        }
};

template class TranslatorPhysBAM<float>;
template <class TS> Log *TranslatorPhysBAM<TS>::log = NULL;
}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
