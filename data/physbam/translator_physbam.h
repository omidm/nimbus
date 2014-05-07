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

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "data/physbam/physbam_include.h"
#include "data/physbam/physbam_data.h"

#include "shared/log.h"
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

        /** Take a FaceArray described by region and read its data from the
         *  PhysicalDataInstance objects in the objects array.
         */
        template<typename T> static void ReadFaceArray(
                const GeometricRegion &region,
                const GeometricRegion &inner,
                const Coord &shift,
                const DataArray &read_set,
                typename PhysBAM::ARRAY<T, FaceIndex>* fa) {
            if (read_set.empty())
                return;
            PhysBAM::ARRAY<T, FaceIndex> flag;
            flag = *fa;
            flag.Fill(0);
            DataArray read_inner;
            DataArray read_outer;
            for (size_t i = 0; i < read_set.size(); ++i) {
                GeometricRegion r = read_set[i]->region();
                if (inner.Covers(&r))
                    read_inner.push_back(read_set[i]);
                else
                    read_outer.push_back(read_set[i]);
            }
            for (size_t i = 0; i < read_inner.size(); ++i) {
                PhysBAMData *data = static_cast<PhysBAMData*>(read_inner[i]);
                Dimension3Vector overlap = GetOverlapSize(data->region(), region);
                if (HasOverlap(overlap)) {
                    T* buffer = reinterpret_cast<T*>(data->buffer());
                    Dimension3Vector src = GetOffset(data->region(), region);
                    Dimension3Vector dest = GetOffset(region, data->region());
                    //  x, y and z values are stored separately due to the
                    // difference in number of x, y and z values in face arrays
                    for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                        int mult_x = 1;
                        int mult_y = data->region().dx();
                        int mult_z = data->region().dy() * data->region().dx();
                        int range_x = overlap(X_COORD);
                        int range_y = overlap(Y_COORD);
                        int range_z = overlap(Z_COORD);
                        int src_offset = 0;
                        switch (dim) {
                            case X_COORD:
                                range_x += 1;
                                mult_y  += 1;
                                mult_z  += data->region().dy();
                                break;
                            case Y_COORD:
                                range_y += 1;
                                mult_z  += data->region().dx();
                                src_offset += (data->region().dx() + 1) *
                                    (data->region().dy()) *
                                    (data->region().dz());
                                break;
                            case Z_COORD:
                                range_z += 1;
                                src_offset += ((data->region().dx()) *
                                        (data->region().dy() + 1) *
                                        (data->region().dz())) +
                                    ((data->region().dx() + 1) *
                                     (data->region().dy()) *
                                     (data->region().dz()));
                                break;
                        }
                        for (int z = 0; z < range_z; ++z) {
                            for (int y = 0; y < range_y; ++y) {
                                for (int x = 0; x < range_x; ++x) {
                                    int source_x = x + src(X_COORD);
                                    int source_y = y + src(Y_COORD);
                                    int source_z = z + src(Z_COORD);
                                    int source_index = source_x * mult_x +
                                        source_y * mult_y +
                                        source_z * mult_z;
                                    source_index += src_offset;
                                    int dest_x = x + dest(X_COORD) + region.x() - shift.x;
                                    int dest_y = y + dest(Y_COORD) + region.y() - shift.y;
                                    int dest_z = z + dest(Z_COORD) + region.z() - shift.z;
                                    typename PhysBAM::VECTOR<int, 3>
                                        destinationIndex(dest_x, dest_y, dest_z);
                                    assert(source_index < data->size() / (int) sizeof(T) && source_index >= 0); // NOLINT
                                    if (flag(dim, destinationIndex) == 0) {
                                        (*fa)(dim, destinationIndex) = buffer[source_index];
                                        flag(dim, destinationIndex) = 1;
                                    } else if (flag(dim, destinationIndex) == 1) {
                                        if ((*fa)(dim, destinationIndex)
                                                != buffer[source_index]) {
                                            (*fa)(dim, destinationIndex) += buffer[source_index];
                                            (*fa)(dim, destinationIndex) /= 2;
                                        }
                                        flag(dim, destinationIndex) = 2;
                                    } else {
                                        // TODO(quhang) needs a more elegant solution.
                                        assert(false);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            int cx1 = inner.x();
            int cx2 = inner.x() + inner.dx();
            int cy1 = inner.y();
            int cy2 = inner.y() + inner.dy();
            int cz1 = inner.z();
            int cz2 = inner.z() + inner.dz();
            for (size_t i = 0; i < read_outer.size(); ++i) {
                PhysBAMData *data = static_cast<PhysBAMData*>(read_outer[i]);
                Dimension3Vector overlap = GetOverlapSize(data->region(), region);
                if (HasOverlap(overlap)) {
                    T* buffer = reinterpret_cast<T*>(data->buffer());
                    Dimension3Vector src = GetOffset(data->region(), region);
                    Dimension3Vector dest = GetOffset(region, data->region());
                    //  x, y and z values are stored separately due to the
                    // difference in number of x, y and z values in face arrays
                    for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                        int mult_x = 1;
                        int mult_y = data->region().dx();
                        int mult_z = data->region().dy() * data->region().dx();
                        int range_x = overlap(X_COORD);
                        int range_y = overlap(Y_COORD);
                        int range_z = overlap(Z_COORD);
                        int src_offset = 0;
                        switch (dim) {
                            case X_COORD:
                                range_x += 1;
                                mult_y  += 1;
                                mult_z  += data->region().dy();
                                break;
                            case Y_COORD:
                                range_y += 1;
                                mult_z  += data->region().dx();
                                src_offset += (data->region().dx() + 1) *
                                    (data->region().dy()) *
                                    (data->region().dz());
                                break;
                            case Z_COORD:
                                range_z += 1;
                                src_offset += ((data->region().dx()) *
                                        (data->region().dy() + 1) *
                                        (data->region().dz())) +
                                    ((data->region().dx() + 1) *
                                     (data->region().dy()) *
                                     (data->region().dz()));
                                break;
                        }
                        for (int z = 0; z < range_z; ++z) {
                            for (int y = 0; y < range_y; ++y) {
                                for (int x = 0; x < range_x; ++x) {
                                    int source_x = x + src(X_COORD);
                                    int source_y = y + src(Y_COORD);
                                    int source_z = z + src(Z_COORD);
                                    int source_index = source_x * mult_x +
                                        source_y * mult_y +
                                        source_z * mult_z;
                                    source_index += src_offset;
                                    int loc_x = x + dest(X_COORD) + region.x();
                                    int loc_y = y + dest(Y_COORD) + region.y();
                                    int loc_z = z + dest(Z_COORD) + region.z();
                                    int dest_x = loc_x - shift.x;
                                    int dest_y = loc_y - shift.y;
                                    int dest_z = loc_z - shift.z;
                                    typename PhysBAM::VECTOR<int, 3>
                                        destinationIndex(dest_x, dest_y, dest_z);
                                    assert(source_index < data->size() / (int) sizeof(T) && source_index >= 0); // NOLINT
                                    if ( (dim == X_COORD && (loc_x == cx1 || loc_x == cx2)) ||
                                         (dim == Y_COORD && (loc_y == cy1 || loc_y == cy2)) ||
                                         (dim == Z_COORD && (loc_z == cz1 || loc_z == cz2)) ) {
                                            typename PhysBAM::VECTOR<int, 3>
                                                destinationIndex(dest_x, dest_y, dest_z);
                                            assert(source_index < data->size() / (int) sizeof(T) && source_index >= 0); // NOLINT
                                            (*fa)(dim, destinationIndex) +=
                                                buffer[source_index];
                                            (*fa)(dim, destinationIndex) /= 2;
                                            flag(dim, destinationIndex) = 2;
                                    } else {
                                        if (flag(dim, destinationIndex) == 0) {
                                            (*fa)(dim, destinationIndex) = buffer[source_index];
                                            flag(dim, destinationIndex) = 1;
                                        } else {
                                            if ((flag(dim, destinationIndex) == 1)) {
                                                typename PhysBAM::VECTOR<int, 3>
                                                    destinationIndex(dest_x, dest_y, dest_z);
                                                assert(source_index < data->size() / (int) sizeof(T) && source_index >= 0); // NOLINT
                                                (*fa)(dim, destinationIndex) +=
                                                    buffer[source_index];
                                                (*fa)(dim, destinationIndex) /= 2;
                                            } else {
                                                assert(false);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /** Take a FaceArray described by region and write it out to the
         *  PhysicalDataInstance objects in the objects array. */
        template<typename T> static void WriteFaceArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                typename PhysBAM::ARRAY<T, FaceIndex>* fa) {
            if (write_set.empty())
                return;
            int_dimension_t region_size = 0;
            region_size += (region.dx() + 1) * region.dy() * region.dz();
            region_size += region.dx() * (region.dy() + 1) * region.dz();
            region_size += region.dx() * region.dy() * (region.dz() + 1);
            if (region_size != fa->buffer_size) {
                dbg(DBG_WARN, "WARN: writing a face array of size %i for a region of size %i and the two sizes should be equal. This check is wrong so you can ignore this warning. I need to determine correct check. -pal\n", fa->buffer_size, region_size);  // NOLINT
            }
            DataArray::const_iterator iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);

                Dimension3Vector overlap = GetOverlapSize(data->region(), region);
                if (!HasOverlap(overlap)) {continue;}

                T* buffer = reinterpret_cast<T*>(data->buffer());

                Dimension3Vector src = GetOffset(region, data->region());
                Dimension3Vector dest = GetOffset(data->region(), region);

                //  x, y and z values are stored separately due to the
                // difference in number of x, y and z values in face arrays
                for (int dim = X_COORD; dim <= Z_COORD; ++dim) {
                    int mult_x = 1;
                    int mult_y = data->region().dx();
                    int mult_z = data->region().dy() * data->region().dx();
                    int range_x = overlap(X_COORD);
                    int range_y = overlap(Y_COORD);
                    int range_z = overlap(Z_COORD);
                    int dst_offset = 0;
                    switch (dim) {
                        case X_COORD:
                            range_x += 1;
                            mult_y  += 1;
                            mult_z  += data->region().dy();
                            break;
                        case Y_COORD:
                            range_y += 1;
                            mult_z  += data->region().dx();
                            dst_offset += (data->region().dx() + 1) *
                                (data->region().dy()) *
                                (data->region().dz());
                            break;
                        case Z_COORD:
                            range_z += 1;
                            dst_offset += ((data->region().dx()) *
                                    (data->region().dy() + 1) *
                                    (data->region().dz())) +
                                ((data->region().dx() + 1) *
                                 (data->region().dy()) *
                                 (data->region().dz()));
                            break;
                    }
                    for (int z = 0; z < range_z; ++z) {
                        for (int y = 0; y < range_y; ++y) {
                            for (int x = 0; x < range_x; ++x) {
                                int dest_x = x + dest(X_COORD);
                                int dest_y = y + dest(Y_COORD);
                                int dest_z = z + dest(Z_COORD);

                                int destination_index = dest_x * mult_x +
                                    dest_y * mult_y +
                                    dest_z * mult_z;
                                destination_index += dst_offset;

                                int source_x = x + src(X_COORD) + region.x() - shift.x;
                                int source_y = y + src(Y_COORD) + region.y() - shift.y;
                                int source_z = z + src(Z_COORD) + region.z() - shift.z;

                                typename PhysBAM::VECTOR<int, 3>
                                    sourceIndex(source_x, source_y, source_z);

                                // The PhysBAM FACE_ARRAY object abstracts away whether
                                // the data is stored in struct of array or array of struct
                                // form (in practice, usually struct of arrays
                                assert(destination_index < data->size() / (int) sizeof(T) && destination_index >= 0); // NOLINT
                                buffer[destination_index] = (*fa)(dim, sourceIndex);
                            }
                        }
                    }
                }
            }
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
            ParticleArray *particles;
            if (positive) {
                particles = &particle_container->positive_particles;
            } else {
                particles = &particle_container->negative_particles;
            }


            for (size_t r = 0; r < regions.size(); ++r) {
                GeometricRegion region = regions[r];

                for (int z = region.z() + 1; z < region.z() + region.dz() - 1; ++z) {
                    for (int y = region.y() + 1; y < region.y() + region.dy() - 1; ++y) {
                        for (int x = region.x() + 1; x < region.x() + region.dx() - 1; ++x) {
                            TV_INT bucket_index(x, y, z);
                            particle_container->Free_Particle_And_Clear_Pointer(
                                    (*particles)(bucket_index));
                        }
                    }
                }

                int xr[2] = {region.x(), region.x() + region.dx() - 1};
                int yr[2] = {region.y(), region.y() + region.dy() - 1};
                int zr[2] = {region.z(), region.z() + region.dz() - 1};

                for (size_t zd = 0; zd < 2; ++zd) {
                    for (size_t yd = 0; yd < 2; ++yd) {
                        for (size_t xd = 0; xd < 2; ++xd) {
                            TV_INT bucket_index(xr[xd], yr[yd], zr[zd]);
                            ParticleBucket *particle_bucket = (*particles)(bucket_index);
                            if (!particle_bucket)
                                continue;

                            ParticleBucket *particle_bucket_history = particle_bucket;
                            (*particles)(bucket_index) = NULL;
                            bool new_particles_alloc = false;
                            ParticleBucket *particle_new_bucket = NULL;

                            while (particle_bucket) {
                                for (int i = 1;
                                        i < particle_bucket->array_collection->Size();
                                        ++i) {
                                    TV particle_position = particle_bucket->X(i);
                                    TV absolute_position = particle_position *
                                        static_cast<float>(scale) + 1.0;
                                    if (absolute_position.x >= region.x() &&
                                        absolute_position.x < region.x() + region.dx() &&
                                        absolute_position.y >= region.y() &&
                                        absolute_position.y < region.y() + region.dy() &&
                                        absolute_position.z >= region.z() &&
                                        absolute_position.z < region.z() + region.dz()) {
                                        if (!new_particles_alloc) {
                                            new_particles_alloc = true;
                                            (*particles)(bucket_index) =
                                                particle_container->
                                                Allocate_Particles(particle_container->template_particles); // NOLINT
                                            particle_new_bucket = (*particles)(bucket_index);
                                        }

                                        int index = particle_container->
                                            Add_Particle(particle_new_bucket);
                                        particle_new_bucket->X(index) = particle_bucket->X(i);
                                        particle_new_bucket->radius(index) =
                                            particle_bucket->radius(i);
                                        particle_new_bucket->quantized_collision_distance(index) =
                                            particle_bucket->quantized_collision_distance(i);
                                        if (particle_container->store_unique_particle_id) {
                                            PhysBAM::ARRAY_VIEW<int> *new_id =
                                                particle_new_bucket->array_collection->
                                                template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                                            PhysBAM::ARRAY_VIEW<int> *id =
                                                particle_bucket->array_collection->
                                                template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
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

        static void DeleteRemovedParticles(
                const Coord &shift,
                const DataArray &da,
                ParticleContainer *particle_container,
                bool positive) {
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
         *
         * More explanation for the region calculation, if,
         *     region = {-1, -1, -1, 12, 12, 12}, shift = {10, 10, 10},
         * then we will copy the data in "instances" limited by global region:
         *     {9, 9, 9, 22, 22, 22}
         * into "particle_container-> for the local region:
         *     {-1, -1, -1, 12, 12, 12}.
         */
        static void ReadParticles(const GeometricRegion &region,
                const Coord &shift,
                const DataArray &read_set,
                ParticleContainer *particle_container,
                const int_dimension_t kScale,
                bool positive,
                bool merge = false) {
            if (read_set.empty())
                return;
            ParticleArray* particles;
            if (positive) {
                particles = &particle_container->positive_particles;
            } else {
                particles = &particle_container->negative_particles;
            }

            // Checks whether the geometric region in the particle array is valid, and
            // clears corresponding buckets inside the geometric region if necessary.
            if (!merge) {
                for (int z = region.z(); z <= region.z() + region.dz(); ++z) {
                    for (int y = region.y(); y <= region.y() + region.dy(); ++y) {
                        for (int x = region.x(); x <= region.x() + region.dx(); ++x) {
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
                ParticleInternal* buffer =
                    reinterpret_cast<ParticleInternal*>(data->buffer());
                ParticleInternal* buffer_end = buffer + static_cast<int>(data->size())
                    / static_cast<int>(sizeof(ParticleInternal));

                for (ParticleInternal* p = buffer; p != buffer_end; ++p) {
                    TV absolute_position;
                    absolute_position.x = p->position[0];
                    absolute_position.y = p->position[1];
                    absolute_position.z = p->position[2];
                    // TODO(quhang) Needs to deal with the particles that lies exactly on
                    // the boundary.
                    if (absolute_position.x >= region.x() + shift.x &&
                            absolute_position.x < region.x() + region.dx() + shift.x &&
                            absolute_position.y >= region.y() + shift.y &&
                            absolute_position.y < region.y() + region.dy() + shift.y &&
                            absolute_position.z >= region.z() + shift.z &&
                            absolute_position.z < region.z() + region.dz() + shift.z) {
                        TV_INT bucket_index(round(absolute_position.x - shift.x),
                                round(absolute_position.y - shift.y),
                                round(absolute_position.z - shift.z));
                        assert(particles->Valid_Index(bucket_index));
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
        }


        /* Writes the particle data from PhysBAM particle container
         * "particle_container->, limited the local region "region", into the
         * corresponding global region of physical instances "instances" after
         * performing the coordinate shifting specified by "shift".
         *
         * "positive" option specifies whether to work on positive particles or
         * negative particles.
         *
         * More explanation for the region calculation, if,
         *     region = {-1, -1, -1, 12, 12, 12}, shift = {10, 10, 10},
         * then we will copy the data of "particle_container-> in the local region:
         *     {-1, -1, -1, 12, 12, 12}
         * into the global region of corresponding "instances":
         *     {9, 9, 9, 22, 22, 22}.
         */
        static void WriteParticles(const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                ParticleContainer *particle_container,
                const int_dimension_t kScale,
                bool positive) {
            if (write_set.empty())
                return;

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

            // Loop through each particle bucket in the specified region.
            for (int z = region.z(); z <= region.z() + region.dz(); ++z)
                for (int y = region.y(); y <= region.y() + region.dy(); ++y)
                    for (int x = region.x(); x <= region.x() + region.dx(); ++x) {
                        TV_INT bucket_index(x, y, z);
                        DataArray::const_iterator iter = write_set.begin();
                        for (; iter != write_set.end(); ++iter) {
                            // Iterate across instances, checking each one.
                            PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                            GeometricRegion data_region = data->region();
                            // TODO(quhang) needs to double check the margin setting.
                            const int_dimension_t kMargin = 1;
                            if (x + shift.x <
                                    data_region.x() - kMargin ||
                                    x + shift.x >
                                    data_region.x() + data_region.dx() + kMargin ||
                                    y + shift.y <
                                    data_region.y() - kMargin ||
                                    y + shift.y >
                                    data_region.y() + data_region.dy() + kMargin ||
                                    z + shift.z <
                                    data_region.z() - kMargin ||
                                    z + shift.z >
                                    data_region.z() + data_region.dz() + kMargin) {
                                continue;
                            }
                            ParticleBucket* particle_bucket = (*particles)(bucket_index);
                            while (particle_bucket) {
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
                                }
                                particle_bucket = particle_bucket->next;
                            }  // Finish looping through all particles.
                        }
                    }

            // Now that we've copied particles into the temporary data buffers,
            // commit the results.
            iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                data->CommitTempBuffer();
            }
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
            if (read_set.empty())
                return;
            RemovedParticleArray* particles;
            if (positive) {
                particles = &particle_container->removed_positive_particles;
            } else {
                particles = &particle_container->removed_negative_particles;
            }

            // Checks whether the geometric region in the particle array is valid, and
            // clears corresponding buckets inside the geometric region if necessary.
            if (!merge) {
                for (int z = region.z(); z <= region.z() + region.dz(); ++z) {
                    for (int y = region.y(); y <= region.y() + region.dy(); ++y) {
                        for (int x = region.x(); x <= region.x() + region.dx(); ++x) {
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
                dbg(DBG_WARN, "Physical data instances are empty.\n");
                return;
            }

            DataArray::const_iterator iter = read_set.begin();
            for (; iter != read_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                RemovedParticleInternal* buffer =
                    reinterpret_cast<RemovedParticleInternal*>(data->buffer());
                RemovedParticleInternal* buffer_end = buffer
                    + static_cast<int>(data->size())
                    / static_cast<int>(sizeof(RemovedParticleInternal));

                for (RemovedParticleInternal* p = buffer; p != buffer_end; ++p) {
                    TV absolute_position;
                    absolute_position.x = p->position[0];
                    absolute_position.y = p->position[1];
                    absolute_position.z = p->position[2];
                    if (absolute_position.x >= region.x() + shift.x &&
                            absolute_position.x < region.x() + region.dx() + shift.x &&
                            absolute_position.y >= region.y() + shift.y &&
                            absolute_position.y < region.y() + region.dy() + shift.y &&
                            absolute_position.z >= region.z() + shift.z &&
                            absolute_position.z < region.z() + region.dz() + shift.z) {
                        TV_INT bucket_index(round(absolute_position.x - shift.x),
                                round(absolute_position.y - shift.y),
                                round(absolute_position.z - shift.z));
                        assert(particles->Valid_Index(bucket_index));
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
            if (write_set.empty())
                return;
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

            // Loop through each particle bucket in the specified region.
            for (int z = region.z(); z <= region.z() + region.dz(); ++z)
                for (int y = region.y(); y <= region.y() + region.dy(); ++y)
                    for (int x = region.x(); x <= region.x() + region.dx(); ++x) {
                        TV_INT bucket_index(x, y, z);
                        RemovedParticleBucket* particle_bucket = (*particles)(bucket_index);
                        while (particle_bucket) {
                            for (int i = 1;
                                    i <= particle_bucket->array_collection->Size();
                                    ++i) {
                                TV particle_position = particle_bucket->X(i);
                                TV absolute_position =
                                    particle_position * (float) kScale + 1.0; // NOLINT
                                DataArray::const_iterator iter = write_set.begin();
                                // Iterate across instances, checking each one.
                                for (; iter != write_set.end(); ++iter) {
                                    PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                                    GeometricRegion data_region = data->region();
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
                                }
                            }
                            particle_bucket = particle_bucket->next;
                        }
                    }

            // Now that we've copied particles into the temporary data buffers,
            // commit the results.
            iter = write_set.begin();
            for (; iter != write_set.end(); ++iter) {
                PhysBAMData* data = static_cast<PhysBAMData*>(*iter);
                data->CommitTempBuffer();
            }
        }

        /* Read scalar array from PhysicalDataInstances specified by instances,
         * limited by the GeometricRegion specified by region, into the
         * ScalarArray specified by dest. This allocates a new scalar array. */
        template<typename T> static void ReadScalarArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &read_set,
                typename PhysBAM::ARRAY<T, TV_INT>* sa) {
            if (read_set.empty())
                return;
            for (size_t i = 0; i < read_set.size(); ++i) {
                PhysBAMData* data = static_cast<PhysBAMData*>(read_set[i]);
                Dimension3Vector overlap = GetOverlapSize(data->region(), region);

                if (HasOverlap(overlap)) {
                    T* buffer  = reinterpret_cast<T*>(data->buffer());

                    Dimension3Vector src  = GetOffset(data->region(), region);
                    Dimension3Vector dest = GetOffset(region, data->region());

                    for (int z = 0; z < overlap(Z_COORD); ++z) {
                        for (int y = 0; y < overlap(Y_COORD); ++y) {
                            for (int x = 0; x < overlap(X_COORD); ++x) {
                                int source_x = x + src(X_COORD);
                                int source_y = y + src(Y_COORD);
                                int source_z = z + src(Z_COORD);
                                int source_index =
                                    (source_z * (data->region().dy() * data->region().dx())) +
                                    (source_y * (data->region().dx())) +
                                    source_x;
                                int dest_x = x + dest(X_COORD) + region.x() - shift.x;
                                int dest_y = y + dest(Y_COORD) + region.y() - shift.y;
                                int dest_z = z + dest(Z_COORD) + region.z() - shift.z;
                                TV_INT destination_index(dest_x, dest_y, dest_z);
                                assert(source_index < data->size() / (int) sizeof(T) && source_index >= 0); // NOLINT
                                (*sa)(destination_index) = buffer[source_index];
                            }
                        }
                    }
                }
            }
        }

        /* Write scalar array data into PhysicalDataInstances specified by instances,
         * limited by the GeometricRegion region. This frees the physbam scalar array. */
        template<typename T> static void WriteScalarArray(
                const GeometricRegion &region,
                const Coord &shift,
                const DataArray &write_set,
                typename PhysBAM::ARRAY<T, TV_INT>* sa) {
            if (write_set.empty())
                return;
            for (size_t i = 0; i < write_set.size(); ++i) {
                PhysBAMData* data = static_cast<PhysBAMData*>(write_set[i]);
                GeometricRegion temp = data->region();
                Dimension3Vector overlap = GetOverlapSize(temp, region);

                if (HasOverlap(overlap)) {
                    T* buffer  = reinterpret_cast<T*>(data->buffer());

                    Dimension3Vector src  = GetOffset(region, data->region());
                    Dimension3Vector dest = GetOffset(data->region(), region);

                    for (int z = 0; z < overlap(Z_COORD); ++z) {
                        for (int y = 0; y < overlap(Y_COORD); ++y) {
                            for (int x = 0; x < overlap(X_COORD); ++x) {
                                int dest_x = x + dest(X_COORD);
                                int dest_y = y + dest(Y_COORD);
                                int dest_z = z + dest(Z_COORD);
                                int destination_index =
                                    (dest_z * (data->region().dy() * data->region().dx())) +
                                    (dest_y * (data->region().dx())) +
                                    dest_x;
                                int source_x = x + src(X_COORD) + region.x() - shift.x;
                                int source_y = y + src(Y_COORD) + region.y() - shift.y;
                                int source_z = z + src(Z_COORD) + region.z() - shift.z;
                                TV_INT source_index(source_x, source_y, source_z);
                                assert(destination_index < data->size() / (int) sizeof(T) && destination_index >= 0); // NOLINT
                                buffer[destination_index] = (*sa)(source_index);
                            }
                        }
                    }
                }
            }
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
}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
