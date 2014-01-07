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
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
#define NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_

#include <algorithm>
#include <cmath>
#include <string>

#include "data/physbam/physbam_include.h"
#include "data/physbam/physbam_data.h"

#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "stdio.h" // NOLINT
#include "worker/physical_data_instance.h"
#include "worker/worker.h"

namespace nimbus {

  template <class VECTOR_TYPE> class TranslatorPhysBAM {
  public:
    typedef VECTOR_TYPE TV;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PhysBAM::GRID<TV> Grid;
    typedef typename TV::SCALAR scalar_t;
    typedef PhysBAM::VECTOR<int_dimension_t, 3> Dimension3Vector;
    typedef PhysBAM::VECTOR<int, 3> Int3Vector;
    typedef typename PhysBAM::FACE_INDEX<TV::dimension> FaceIndex;
    typedef typename PhysBAM::ARRAY<scalar_t, FaceIndex> FaceArray;

    // typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> Particles;
    // typedef typename PhysBAM::ARRAY<Particles*, PhysBAM::VECTOR<int, 3> > ParticlesArray;
    //  typedef typename PhysBAM::ARRAY<Particles*, int> ParticlesArray;

    // Container class for particles and removed particles.
    typedef typename PhysBAM::PARTICLE_LEVELSET_UNIFORM<Grid> ParticleContainer;
    // Particle bucket. Each node has a linked list of particle buckets.
    typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> ParticleBucket;
    // Particle array, indexed by node.
    typedef typename PhysBAM::ARRAY<ParticleBucket*, TV_INT> ParticleArray;

    // TODO(quhang) Change the name.
    typedef typename PhysBAM::PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>
        RemovedParticleBucket;
    typedef typename PhysBAM::ARRAY<RemovedParticleBucket*, TV_INT>
        RemovedParticleArray;

    typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleLevelset;
    typedef typename PhysBAM::ARRAY<scalar_t, Int3Vector> ScalarArray;

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
      int_dimension_t index[3];
      scalar_t delta[3];
      scalar_t radius;
      uint16_t quantized_collision_distance;
      int32_t id;
    };

    struct RemovedParticleInternal : public ParticleInternal {
      scalar_t v[3];
    };

    virtual void ReadFaceArray(const GeometricRegion* region,
                               const PdiVector* objects,
                               FaceArray* fa) {
      if (objects != NULL) {
        PdiVector::const_iterator iter = objects->begin();
        for (; iter != objects->end(); ++iter) {
          const PhysicalDataInstance* obj = *iter;
          Dimension3Vector overlap = GetOverlapSize(obj->region(), region);
          if (HasOverlap(overlap)) {
            std::string reg_str = region->toString();
            dbg(DBG_TRANSLATE, "Incorporating physical object %lu into FaceArray for region %s.\n", obj->id(), reg_str.c_str()); // NOLINT`
            PhysBAMData* data = static_cast<PhysBAMData*>(obj->data());
            scalar_t* buffer = reinterpret_cast<scalar_t*>(data->buffer());

            Dimension3Vector src = GetOffset(obj->region(), region);
            Dimension3Vector dest = GetOffset(region, obj->region());

            //  x, y and z values are stored separately due to the
            // difference in number of x, y and z values in face arrays
            for (int dim = X_COORD; dim <= Z_COORD; dim++) {
              int mult_x = 1;
              int mult_y = obj->region()->dx();
              int mult_z = obj->region()->dy() * obj->region()->dx();
              int range_x = overlap(X_COORD);
              int range_y = overlap(Y_COORD);
              int range_z = overlap(Z_COORD);
              int src_offset = 0;
              switch (dim) {
                  case X_COORD:
                    range_x += 1;
                    mult_y  += 1;
                    mult_z  += obj->region()->dy();
                    break;
                  case Y_COORD:
                    range_y += 1;
                    mult_z  += obj->region()->dx();
                    src_offset += (obj->region()->dx() + 1) *
                                  (obj->region()->dy()) *
                                  (obj->region()->dz());
                    break;
                  case Z_COORD:
                    range_z += 1;
                    src_offset += ((obj->region()->dx()) *
                                  (obj->region()->dy()+1) *
                                  (obj->region()->dz())) +
                                  ((obj->region()->dx() + 1) *
                                  (obj->region()->dy()) *
                                  (obj->region()->dz()));
                    break;
                }
              for (int z = 0; z < range_z; z++) {
                for (int y = 0; y < range_y; y++) {
                  for (int x = 0; x < range_x; x++) {
                    int source_x = x + src(X_COORD);
                    int source_y = y + src(Y_COORD);
                    int source_z = z + src(Z_COORD);

                    int source_index = source_x * mult_x +
                                       source_y * mult_y +
                                       source_z * mult_z;
                    source_index += src_offset;

                    int dest_x = x + dest(X_COORD) + region->x();
                    int dest_y = y + dest(Y_COORD) + region->y();
                    int dest_z = z + dest(Z_COORD) + region->z();

                    typename PhysBAM::VECTOR<int, 3>
                      destinationIndex(dest_x, dest_y, dest_z);

                    // The PhysBAM FACE_ARRAY object abstracts away whether
                    // the data is stored in struct of array or array of struct
                    // form (in practice, usually struct of arrays.
                    assert(source_index < data->size() && source_index >= 0);
                    (*fa)(dim, destinationIndex) = buffer[source_index];
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
    virtual bool WriteFaceArray(const GeometricRegion* region,
                                PdiVector* objects,
                                FaceArray* fa) {
      int_dimension_t region_size = 0;
      region_size += (region->dx() + 1) * region->dy() * region->dz();
      region_size += region->dx() * (region->dy() + 1) * region->dz();
      region_size += region->dx() * region->dy() * (region->dz() + 1);
      if (region_size != fa->buffer_size) {
        dbg(DBG_WARN, "WARN: writing a face array of size %i for a region of size %i and the two sizes should be equal. This check is wrong so you can ignore this warning. I need to determine correct check. -pal\n", fa->buffer_size, region_size);  // NOLINT
        //  return false;
      }

      if (objects != NULL) {
        PdiVector::iterator iter = objects->begin();

        // Loop over the Nimbus objects, copying the relevant PhysBAM data
        // into each one
        for (; iter != objects->end(); ++iter) {
          const PhysicalDataInstance* obj = *iter;
          Dimension3Vector overlap = GetOverlapSize(obj->region(), region);
          if (!HasOverlap(overlap)) {continue;}

          dbg(DBG_TRANSLATE, "Saving FaceArray into physical object %lu.\n", obj->id());
          PhysBAMData* data = static_cast<PhysBAMData*>(obj->data());
          scalar_t* buffer = reinterpret_cast<scalar_t*>(data->buffer());

          Dimension3Vector src = GetOffset(region, obj->region());
          Dimension3Vector dest = GetOffset(obj->region(), region);

          //  x, y and z values are stored separately due to the
          // difference in number of x, y and z values in face arrays
          for (int dim = X_COORD; dim <= Z_COORD; dim++) {
            int mult_x = 1;
            int mult_y = obj->region()->dx();
            int mult_z = obj->region()->dy() * obj->region()->dx();
            int range_x = overlap(X_COORD);
            int range_y = overlap(Y_COORD);
            int range_z = overlap(Z_COORD);
            int dst_offset = 0;
            switch (dim) {
                case X_COORD:
                  range_x += 1;
                  mult_y  += 1;
                  mult_z  += obj->region()->dy();
                  break;
                case Y_COORD:
                  range_y += 1;
                  mult_z  += obj->region()->dx();
                  dst_offset += (obj->region()->dx() + 1) *
                                (obj->region()->dy()) *
                                (obj->region()->dz());
                  break;
                case Z_COORD:
                  range_z += 1;
                  dst_offset += ((obj->region()->dx()) *
                                (obj->region()->dy()+1) *
                                (obj->region()->dz())) +
                                ((obj->region()->dx() + 1) *
                                (obj->region()->dy()) *
                                (obj->region()->dz()));
                  break;
              }
            for (int z = 0; z < range_z; z++) {
              for (int y = 0; y < range_y; y++) {
                for (int x = 0; x < range_x; x++) {
                  int dest_x = x + dest(X_COORD);
                  int dest_y = y + dest(Y_COORD);
                  int dest_z = z + dest(Z_COORD);

                  int destination_index = dest_x * mult_x +
                                          dest_y * mult_y +
                                          dest_z * mult_z;
                  destination_index += dst_offset;

                  int source_x = x + src(X_COORD) + region->x();
                  int source_y = y + src(Y_COORD) + region->y();
                  int source_z = z + src(Z_COORD) + region->z();

                  typename PhysBAM::VECTOR<int, 3>
                    sourceIndex(source_x, source_y, source_z);

                  // The PhysBAM FACE_ARRAY object abstracts away whether
                  // the data is stored in struct of array or array of struct
                  // form (in practice, usually struct of arrays
                  assert(destination_index < data->size() && destination_index >= 0);
                  buffer[destination_index] = (*fa)(dim, sourceIndex);
                }
              }
            }
          }
        }
      }
      return true;
    }

    /* Read the particles from the PhysicalDataInstances specified
     * by instances, limited by the GeometricRegion specified by region,
     * into the PhysBAM::PARTICLE_LEVELSET_UNIFORM specified by dest.
     * This will clear out any existing data in particles first if "merge" flag
     * is not set. */
    bool ReadParticles(const GeometricRegion* region,
                       const PdiVector* instances,
                       ParticleContainer& particle_container,
                       bool positive,
                       bool merge = false) {
      ParticleArray* particles;
      if (positive) {
        particles = &particle_container.positive_particles;
      } else {
        particles = &particle_container.negative_particles;
      }

      // Checks whether the geometric region in the particle array is valid, and
      // clears corresponding buckets inside the geometric region if necessary.
      for (int z = region->z(); z < region->z() + region->dz(); z++)
        for (int y = region->y(); y < region->y() + region->dy(); y++)
          for (int x = region->x(); x < region->x() + region->dx(); x++) {
            TV_INT bucket_index(x, y, z);
            if (!particles->Valid_Index(bucket_index)) {
              dbg(DBG_WARN, "Bucket index (%d, %d, %d) out of range.\n",
                            x, y, z);
              // Warning: might be too strict.
              return false;
            }
            if (!merge) {
              particle_container.Free_Particle_And_Clear_Pointer(
                  (*particles)(bucket_index));
            }
          }

      if (instances == NULL) {
        dbg(DBG_WARN, "Physical data instances are empty.\n");
        return false;
      }

      PdiVector::const_iterator iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        ParticleInternal* buffer =
            reinterpret_cast<ParticleInternal*>(data->buffer());
        ParticleInternal* buffer_end = buffer + static_cast<int>(data->size())
             / static_cast<int>(sizeof(ParticleInternal));

        for (ParticleInternal* p = buffer; p != buffer_end; ++p) {
          int_dimension_t xi = p->index[0];
          int_dimension_t yi = p->index[1];
          int_dimension_t zi = p->index[2];
          if (xi >= region->x() &&
              xi < region->x()+region->dx() &&
              yi >= region->y() &&
              yi < region->y()+region->dy() &&
              zi >= region->z() &&
              zi < region->z()+region->dz()) {
            TV_INT bucket_index(xi, yi, zi);
            assert(particles->Valid_Index(bucket_index));
            // NOTE(By Chinmayee): Please comment out these changes and don't
            // delete them when pushing any updates, till we verify that the
            // code works correctly, and does not give any assertion failure or
            // seg fault on both Linux and Mac.
            if (!(*particles)(bucket_index)) {
              (*particles)(bucket_index) = particle_container.
                  Allocate_Particles(particle_container.template_particles);
            }
            ParticleBucket* particle_bucket =
                (*particles)(bucket_index);

            // Note that Add_Particle traverses a linked list of particle
            // buckets, so it's O(N^2) time. Blech.
            // I think things might not be that bad, because the number of
            // bucket is expected to be constant.  -- quhang
            int index = particle_container.Add_Particle(particle_bucket);
            particle_bucket->X(index) = VECTOR_TYPE(p->delta[0],
                                                    p->delta[1],
                                                    p->delta[2]);
            particle_bucket->radius(index) = p->radius;
            particle_bucket->quantized_collision_distance(index) =
              p->quantized_collision_distance;
            if (particle_container.store_unique_particle_id) {
              PhysBAM::ARRAY_VIEW<int>* id = particle_bucket->array_collection->
                  template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
              (*id)(index) = p->id;
            }
          }
        }  // End the loop for buffer.
      }
      return true;
    }


    /* Write the Particle data in particles into the
     * PhysicalDataInstances specified by instances, limited by the
     * GeometricRegion region. */
    bool WriteParticles(const GeometricRegion* region,
                        PdiVector* instances,
                        ParticleContainer& particle_container,
                        bool positive) {
      PdiVector::iterator iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        data->ClearTempBuffer();
      }

      ParticleArray* particles;
      if (positive) {
        particles = &particle_container.positive_particles;
      } else {
        particles = &particle_container.negative_particles;
      }

      // Loop through each particle bucket in the specified region.
      for (int z = region->z(); z < region->z() + region->dz(); z++)
        for (int y = region->y(); y < region->y() + region->dy(); y++)
          for (int x = region->x(); x < region->x() + region->dx(); x++) {
            TV_INT bucket_index(x, y, z);
            if (!particles->Valid_Index(bucket_index)) {
              dbg(DBG_WARN, "Bucket index (%d, %d, %d) out of range.\n",
                            x, y, z);
              return false;
            }
            ParticleBucket* particle_bucket = (*particles)(bucket_index);
            while (particle_bucket) {
              for (int i = 1; i <= particle_bucket->array_collection->Size();
                   i++) {
                PdiVector::iterator iter = instances->begin();
                // Iterate across instances, checking each one.
                for (; iter != instances->end(); ++iter) {
                  const PhysicalDataInstance* instance = *iter;
                  GeometricRegion* instanceRegion = instance->region();
                  // If it's inside the region of the physical data instance.
                  if (x >= instanceRegion->x() &&
                      x < (instanceRegion->x() + instanceRegion->dx()) &&
                      y >= instanceRegion->y() &&
                      y < (instanceRegion->y() + instanceRegion->dy()) &&
                      z >= instanceRegion->z() &&
                      z < (instanceRegion->z() + instanceRegion->dz())) {
                    ParticleInternal particle_buffer;
                    particle_buffer.index[0] = x;
                    particle_buffer.index[1] = y;
                    particle_buffer.index[2] = z;
                    VECTOR_TYPE particle_delta = particle_bucket->X(i);
                    particle_buffer.delta[0] = particle_delta.x;
                    particle_buffer.delta[1] = particle_delta.y;
                    particle_buffer.delta[2] = particle_delta.z;
                    particle_buffer.radius = particle_bucket->radius(i);
                    particle_buffer.quantized_collision_distance =
                        particle_bucket->quantized_collision_distance(i);
                    if (particle_container.store_unique_particle_id) {
                      PhysBAM::ARRAY_VIEW<int>* id =
                          particle_bucket->array_collection->
                          template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                      particle_buffer.id = (*id)(i);
                    }
                    PhysBAMData* data =
                        static_cast<PhysBAMData*>(instance->data());
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
      iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        data->CommitTempBuffer();
      }
      return true;
    }

    /* Read the removed particles from the PhysicalDataInstances specified
     * by instances, limited by the GeometricRegion specified by region,
     * into the PhysBAM::PARTICLE_LEVELSET_UNIFORM specified by dest.
     * This will clear out any existing data in particles first if "merge" flag
     * is not set. */
    bool ReadRemovedParticles(const GeometricRegion* region,
                              const PdiVector* instances,
                              ParticleContainer& particle_container,
                              bool positive,
                              bool merge = false) {
      RemovedParticleArray* particles;
      if (positive) {
        particles = &particle_container.removed_positive_particles;
      } else {
        particles = &particle_container.removed_negative_particles;
      }

      // Checks whether the geometric region in the particle array is valid, and
      // clears corresponding buckets inside the geometric region if necessary.
      for (int z = region->z(); z < region->z() + region->dz(); z++)
        for (int y = region->y(); y < region->y() + region->dy(); y++)
          for (int x = region->x(); x < region->x() + region->dx(); x++) {
            TV_INT bucket_index(x, y, z);
            if (!particles->Valid_Index(bucket_index)) {
              dbg(DBG_WARN, "Bucket index (%d, %d, %d) out of range.\n",
                            x, y, z);
              // Warning: might be too strict.
              return false;
            }
            if (!merge) {
              if ((*particles)(bucket_index)) {
                delete (*particles)(bucket_index);
                (*particles)(bucket_index) = NULL;
              }
            }
          }

      if (instances == NULL) {
        dbg(DBG_WARN, "Physical data instances are empty.\n");
        return false;
      }

      PdiVector::const_iterator iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        RemovedParticleInternal* buffer =
            reinterpret_cast<RemovedParticleInternal*>(data->buffer());
        RemovedParticleInternal* buffer_end = buffer
            + static_cast<int>(data->size())
            / static_cast<int>(sizeof(RemovedParticleInternal));

        for (RemovedParticleInternal* p = buffer; p != buffer_end; ++p) {
          int_dimension_t xi = p->index[0];
          int_dimension_t yi = p->index[1];
          int_dimension_t zi = p->index[2];
          if (xi >= region->x() &&
              xi < region->x()+region->dx() &&
              yi >= region->y() &&
              yi < region->y()+region->dy() &&
              zi >= region->z() &&
              zi < region->z()+region->dz()) {
            TV_INT bucket_index(xi, yi, zi);
            assert(particles->Valid_Index(bucket_index));
            // NOTE(By Chinmayee): Please comment out these changes and don't
            // delete them when pushing any updates, till we verify that the
            // code works correctly, and does not give any assertion failure or
            // seg fault on both Linux and Mac.
            if (!(*particles)(bucket_index)) {
              (*particles)(bucket_index) =
                  particle_container.Allocate_Particles(
                  particle_container.template_removed_particles);
            }
            RemovedParticleBucket* particle_bucket =
                (*particles)(bucket_index);

            // Note that Add_Particle traverses a linked list of particle
            // buckets, so it's O(N^2) time. Blech.
            // I think things might not be that bad, because the number of
            // bucket is expected to be constant.  -- quhang
            int index = particle_container.Add_Particle(particle_bucket);
            particle_bucket->X(index) = VECTOR_TYPE(p->delta[0],
                                                    p->delta[1],
                                                    p->delta[2]);
            particle_bucket->radius(index) = p->radius;
            particle_bucket->quantized_collision_distance(index) =
              p->quantized_collision_distance;
            if (particle_container.store_unique_particle_id) {
              PhysBAM::ARRAY_VIEW<int>* id = particle_bucket->array_collection->
                  template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
              (*id)(index) = p->id;
            }
            particle_bucket->V(index) = VECTOR_TYPE(p->v[0],
                                                    p->v[1],
                                                    p->v[2]);
          }
        }  // End the loop for buffer.
      }
      return true;
    }

    /* Write the Particle data in removed particles into the
     * PhysicalDataInstances specified by instances, limited by the
     * GeometricRegion region. */
    bool WriteRemovedParticles(const GeometricRegion* region,
                        PdiVector* instances,
                        ParticleContainer& particle_container,
                        bool positive) {
      PdiVector::iterator iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        data->ClearTempBuffer();
      }

      RemovedParticleArray* particles;
      if (positive) {
        particles = &particle_container.removed_positive_particles;
      } else {
        particles = &particle_container.removed_negative_particles;
      }

      // Loop through each particle bucket in the specified region.
      for (int z = region->z(); z < region->z() + region->dz(); z++)
        for (int y = region->y(); y < region->y() + region->dy(); y++)
          for (int x = region->x(); x < region->x() + region->dx(); x++) {
            TV_INT bucket_index(x, y, z);
            if (!particles->Valid_Index(bucket_index)) {
              dbg(DBG_WARN, "Bucket index (%d, %d, %d) out of range.\n",
                            x, y, z);
              return false;
            }
            RemovedParticleBucket* particle_bucket = (*particles)(bucket_index);
            // Note: the outer while loop might not be necessary for removed
            // particle bucket.
            while (particle_bucket) {
              for (int i = 1; i <= particle_bucket->array_collection->Size();
                   i++) {
                PdiVector::iterator iter = instances->begin();
                // Iterate across instances, checking each one.
                for (; iter != instances->end(); ++iter) {
                  const PhysicalDataInstance* instance = *iter;
                  GeometricRegion* instanceRegion = instance->region();
                  // If it's inside the region of the physical data instance.
                  if (x >= instanceRegion->x() &&
                      x < (instanceRegion->x() + instanceRegion->dx()) &&
                      y >= instanceRegion->y() &&
                      y < (instanceRegion->y() + instanceRegion->dy()) &&
                      z >= instanceRegion->z() &&
                      z < (instanceRegion->z() + instanceRegion->dz())) {
                    RemovedParticleInternal particle_buffer;
                    particle_buffer.index[0] = x;
                    particle_buffer.index[1] = y;
                    particle_buffer.index[2] = z;
                    VECTOR_TYPE particle_delta = particle_bucket->X(i);
                    particle_buffer.delta[0] = particle_delta.x;
                    particle_buffer.delta[1] = particle_delta.y;
                    particle_buffer.delta[2] = particle_delta.z;
                    particle_buffer.radius = particle_bucket->radius(i);
                    particle_buffer.quantized_collision_distance =
                        particle_bucket->quantized_collision_distance(i);
                    if (particle_container.store_unique_particle_id) {
                      PhysBAM::ARRAY_VIEW<int>* id =
                          particle_bucket->array_collection->
                          template Get_Array<int>(PhysBAM::ATTRIBUTE_ID_ID);
                      particle_buffer.id = (*id)(i);
                    }
                    particle_buffer.v[0] = particle_bucket->V(i).x;
                    particle_buffer.v[1] = particle_bucket->V(i).y;
                    particle_buffer.v[2] = particle_bucket->V(i).z;
                    PhysBAMData* data =
                        static_cast<PhysBAMData*>(instance->data());
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
      iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        data->CommitTempBuffer();
      }
      return true;
    }

    /* Read scalar array from PhysicalDataInstances specified by instances,
     * limited by the GeometricRegion specified by region, into the
     * ScalarArray specified by dest. This allocates a new scalar array. */
    virtual ScalarArray* ReadScalarArray(const GeometricRegion* region,
                                         const PdiVector* instances,
                                         ScalarArray *sa) {
        if (instances != NULL) {
            PdiVector::const_iterator iter = instances->begin();
            for (; iter != instances->end(); ++iter) {
                const PhysicalDataInstance* inst = *iter;
                Dimension3Vector overlap = GetOverlapSize(inst->region(), region);

                if (HasOverlap(overlap)) {
                    dbg(DBG_TRANSLATE, "Incorporating physical object %lu into scalar array.\n",
                            inst->id());
                    PhysBAMData* data = static_cast<PhysBAMData*>(inst->data());
                    scalar_t* buffer  = reinterpret_cast<scalar_t*>(data->buffer());

                    Dimension3Vector src  = GetOffset(inst->region(), region);
                    Dimension3Vector dest = GetOffset(region, inst->region());

                    for (int z = 0; z < overlap(Z_COORD); z++) {
                        for (int y = 0; y < overlap(Y_COORD); y++) {
                            for (int x = 0; x < overlap(X_COORD); x++) {
                                int source_x = x + src(X_COORD);
                                int source_y = y + src(Y_COORD);
                                int source_z = z + src(Z_COORD);
                                int source_index =
                                    (source_z * (inst->region()->dy() * inst->region()->dx())) +
                                    (source_y * (inst->region()->dx())) +
                                    source_x;
                                int dest_x = x + dest(X_COORD) + region->x();
                                int dest_y = y + dest(Y_COORD) + region->y();
                                int dest_z = z + dest(Z_COORD) + region->z();
                                Int3Vector destination_index(dest_x, dest_y, dest_z);
                                assert(source_index < data->size() && source_index >= 0);
                                (*sa)(destination_index) = buffer[source_index];
                            }
                        }
                    }
                }
            }
        }
        return sa;
    }

    /* Write scalar array data into PhysicalDataInstances specified by instances,
     * limited by the GeometricRegion region. This frees the physbam scalar array. */
    virtual bool WriteScalarArray(const GeometricRegion* region,
                                  PdiVector* instances,
                                  ScalarArray* sa) {
        if (sa->counts != Int3Vector(region->dx(), region->dy(), region->dz())) {
            dbg(DBG_WARN, "WARN: writing to a scalar array of a different size\n");
            Int3Vector cp = sa->counts;
            dbg(DBG_WARN, "WARN: physbam array has size %i, %i, %i\n", cp.x, cp.y, cp.z);
            dbg(DBG_WARN, "WARN: original array has size %llu, %llu, %llu\n", region->dx(), region->dy(), region->dz()); // NOLINT
        }

        if (instances != NULL) {
            PdiVector::iterator iter = instances->begin();

            for (; iter != instances->end(); ++iter) {
                const PhysicalDataInstance* inst = *iter;
                GeometricRegion *temp = inst->region();
                Dimension3Vector overlap = GetOverlapSize(temp, region);

                if (HasOverlap(overlap)) {
                    dbg(DBG_TRANSLATE, "Saving scalar array into physical object %lu.\n",
                                       inst->id());
                    PhysBAMData* data = static_cast<PhysBAMData*>(inst->data());
                    scalar_t* buffer  = reinterpret_cast<scalar_t*>(data->buffer());

                    Dimension3Vector src  = GetOffset(region, inst->region());
                    Dimension3Vector dest = GetOffset(inst->region(), region);

                    for (int z = 0; z < overlap(Z_COORD); z++) {
                        for (int y = 0; y < overlap(Y_COORD); y++) {
                            for (int x = 0; x < overlap(X_COORD); x++) {
                                int dest_x = x + dest(X_COORD);
                                int dest_y = y + dest(Y_COORD);
                                int dest_z = z + dest(Z_COORD);
                                int destination_index =
                                    (dest_z * (inst->region()->dy() * inst->region()->dx())) +
                                    (dest_y * (inst->region()->dx())) +
                                    dest_x;
                                int source_x = x + src(X_COORD) + region->x();
                                int source_y = y + src(Y_COORD) + region->y();
                                int source_z = z + src(Z_COORD) + region->z();
                                Int3Vector source_index(source_x, source_y, source_z);
                                assert(destination_index < data->size() && destination_index >= 0);
                                buffer[destination_index] = (*sa)(source_index);
                            }
                        }
                    }
                }
            }
        }
        return true;
    }


  private:
    /* Return a vector describing what the offset of b
       within a, such that a.x + offset = b.x. If
       offset is negative, return 0. */
    virtual Dimension3Vector GetOffset(const GeometricRegion* a,
                                       const GeometricRegion* b) {
      Dimension3Vector result;

      // If source is > than dest, its offset is zero (it's contained),
      // otherwise the offset is the difference between the values.
      int_dimension_t x = b->x() - a->x();
      int_dimension_t y = b->y() - a->y();
      int_dimension_t z = b->z() - a->z();
      result(X_COORD) = (x >= 0)? x:0;
      result(Y_COORD) = (y >= 0)? y:0;
      result(Z_COORD) = (z >= 0)? z:0;

      return result;
    }

    virtual Dimension3Vector GetOverlapSize(const GeometricRegion* src,
                                            const GeometricRegion* dest) {
      Dimension3Vector result;

      int_dimension_t x_start = std::max(src->x(), dest->x());
      int_dimension_t x_end   = std::min(src->x() + src->dx(),
                                         dest->x() + dest->dx());
      int_dimension_t x_size = x_end - x_start;

      int_dimension_t y_start = std::max(src->y(), dest->y());
      int_dimension_t y_end   = std::min(src->y() + src->dy(),
                                         dest->y() + dest->dy());
      int_dimension_t y_size = y_end - y_start;

      int_dimension_t z_start = std::max(src->y(), dest->z());
      int_dimension_t z_end   = std::min(src->z() + src->dz(),
                                         dest->z() + dest->dz());
      int_dimension_t z_size = z_end - z_start;

      result(X_COORD) = (x_size >= 0)? x_size:0;
      result(Y_COORD) = (y_size >= 0)? y_size:0;
      result(Z_COORD) = (z_size >= 0)? z_size:0;

      return result;
    }

    virtual bool HasOverlap(Dimension3Vector overlapSize) {
      return (overlapSize(X_COORD) > 0 &&
              overlapSize(Y_COORD) > 0 &&
              overlapSize(Z_COORD) > 0);
    }
};

}  // namespace nimbus

#endif  // NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_H_
