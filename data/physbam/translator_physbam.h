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
    typedef PhysBAM::GRID<TV> Grid;
    typedef typename TV::SCALAR scalar_t;
    typedef PhysBAM::VECTOR<int_dimension_t, 3> Dimension3Vector;
    typedef PhysBAM::VECTOR<int, 3> Int3Vector;
    typedef typename PhysBAM::FACE_INDEX<TV::dimension> FaceIndex;
    typedef typename PhysBAM::ARRAY<scalar_t, FaceIndex> FaceArray;
    typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> Particles;

    typedef typename PhysBAM::ARRAY<Particles*, PhysBAM::VECTOR<int, 3> > ParticlesArray;

    //  typedef typename PhysBAM::ARRAY<Particles*, int> ParticlesArray;

    typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleLevelset;
    typedef typename PhysBAM::ARRAY<scalar_t, Int3Vector> ScalarArray;

    enum {
      X_COORD = 1,
      Y_COORD = 2,
      Z_COORD = 3
    };

    explicit TranslatorPhysBAM() {}
    virtual ~TranslatorPhysBAM() {}

    virtual void ReadFaceArray(const GeometricRegion* region,
                               const PdiVector* objects,
                               FaceArray* fa) {
      Dimension3Vector vec;
      vec(X_COORD) = region->dx();
      vec(Y_COORD) = region->dy();
      vec(Z_COORD) = region->dz();

      // Create a FACE_ARRAY of the right size.
      PhysBAM::RANGE<PhysBAM::VECTOR<int, 3> > range(1, region->dx(),
                                                     1, region->dy(),
                                                     1, region->dz());

      fa->Resize(range);

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

            Dimension3Vector dest = GetOffset(region, obj->region());
            Dimension3Vector src = GetOffset(obj->region(), region);

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

                    int dest_x = x + dest(X_COORD) + 1;
                    int dest_y = y + dest(Y_COORD) + 1;
                    int dest_z = z + dest(Z_COORD) + 1;

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

                  int source_x = x + src(X_COORD) + 1;
                  int source_y = y + src(Y_COORD) + 1;
                  int source_z = z + src(Z_COORD) + 1;

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
     * into the PhysBAM::PARTICLE_LEVELSET_PARTICLES specified by dest.
     * This will clear out any existing data in particles first. */
    virtual bool ReadParticles(GeometricRegion* region,
                               CPdiVector* instances,
                               ParticleLevelset& particleLevelset,
                               bool positive) {
      for (int z = 1; z < region->dz(); z++) {
        for (int y = 1; y < region->dy(); y++) {
          for (int x = 1; x < region->dx(); x++) {
            Particles* particleBucket;
            // printf("Offset of positive particles is %i\n", offsetof(ParticleLevelset, positive_particles)); // NOLINT
            if (positive) {
              // ParticlesArray& array = particleLevelset.get_positive_particles();
              ParticlesArray* arrayPtr =
                reinterpret_cast<ParticlesArray*>((char*)&particleLevelset.positive_particles + 88); // NOLINT
              ParticlesArray& array = *arrayPtr;
              // printf("Pointer to positive particles is %p\n", arrayPtr);
              particleBucket = array(Int3Vector(x, y, z));
            } else {
              particleBucket = particleLevelset.negative_particles(Int3Vector(x, y, z));
            }

            while (particleBucket != NULL) {
              particleBucket->array_collection->Delete_All_Elements();
              particleBucket = particleBucket->next;
            }
          }
        }
      }

      if (instances == NULL) {
        dbg(DBG_WARN, "Tried to read particles from a NULL vector of PhysicalDataInstances\n");
        return false;
      }

      CPdiVector::iterator iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        scalar_t* buffer = reinterpret_cast<scalar_t*>(data->buffer());

        for (int i = 0; i < data->size(); i+= 5) {
          // This instance may have particles inside the region. So
          // we have to iterate over them and insert them accordingly.
          // Particles are stored as (x,y,z) triples in data array
          scalar_t x = buffer[i];
          scalar_t y = buffer[i + 1];
          scalar_t z = buffer[i + 2];
          scalar_t radius = buffer[i + 3];

          uint16_t collision_distance = (uint16_t)buffer[i + 4];  // NOLINT

          VECTOR_TYPE position;
          position.x = x;
          position.y = y;
          position.z = z;
          int_dimension_t xi = (int_dimension_t)x;
          int_dimension_t yi = (int_dimension_t)y;
          int_dimension_t zi = (int_dimension_t)z;

          // If particle is within region, copy it to particles
          if (xi >= region->x() &&
              xi <= (region->x() + region->dx()) &&
              yi >= region->y() &&
              yi <= (region->y() + region->dy()) &&
              zi >= region->z() &&
              zi <= (region->z() + region->dz())) {
            // I have no idea if this works -- indexing by a 3D vector?
            Particles* cellParticles;
            if (positive) {
              ParticlesArray* arrayPtr =
                reinterpret_cast<ParticlesArray*>((char*)&particleLevelset.positive_particles + 88); // NOLINT
              ParticlesArray& array = *arrayPtr;
              printf("Pointer to positive particles is %p\n", arrayPtr);
              cellParticles = array(Int3Vector(xi - region->x(),
                                               yi - region->y(),
                                               zi - region->z()));
              // cellParticles = particleLevelset.positive_particles(Int3Vector(x, y, z));
            } else {
              cellParticles = particleLevelset.negative_particles(Int3Vector(x, y, z));
            }

            // Note that Add_Particle traverses a linked list of particle
            // buckets, so it's O(N^2) time. Blech.
            int index = particleLevelset.Add_Particle(cellParticles);
            cellParticles->quantized_collision_distance(index) =
              collision_distance;
            cellParticles->X(index) = position;
            cellParticles->radius(index) = radius;
          }
        }
      }
      return true;
    }

    /* Write the Particle data in particles into the
     * PhysicalDataInstances specified by instances, limited by the
     * GeometricRegion region. */

    virtual bool WriteParticles(GeometricRegion* region,
                                CPdiVector* instances,
                                ParticleLevelset& particleLevelset,
                                bool positive) {
      CPdiVector::iterator iter = instances->begin();
      for (; iter != instances->end(); ++iter) {
        const PhysicalDataInstance* instance = *iter;
        PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
        data->ClearTempBuffer();
      }

      for (int z = 0; z < region->dz(); z++) {
        for (int y = 0; y < region->dy(); y++) {
          for (int x = 0; x < region->dx(); x++) {
            ParticlesArray* arrayPtr =
              reinterpret_cast<ParticlesArray*>((char*)&particleLevelset.positive_particles + 88); // NOLINT
            ParticlesArray& array = *arrayPtr;
            printf("Pointer to positive particles is %p\n", arrayPtr);
            Particles* particles = array(Int3Vector(x, y, z));
            while (particles) {
              for (int i = 1; i <= particles->array_collection->Size(); i++) {
                VECTOR_TYPE particle = particles->X(i);
                scalar_t x = particle.x;
                scalar_t y = particle.y;
                scalar_t z = particle.z;
                int_dimension_t xi = (int_dimension_t)x;
                int_dimension_t yi = (int_dimension_t)y;
                int_dimension_t zi = (int_dimension_t)z;
                // If it's inside the region,
                if (xi >= region->x() &&
                    xi <= (region->x() + region->dx()) &&
                    yi >= region->y() &&
                    yi <= (region->y() + region->dy()) &&
                    zi >= region->z() &&
                    zi <= (region->z() + region->dz())) {
                  CPdiVector::iterator iter = instances->begin();

                  // Iterate across instances, checking each one
                  for (; iter != instances->end(); ++iter) {
                    const PhysicalDataInstance* instance = *iter;
                    GeometricRegion* instanceRegion = instance->region();
                    PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());

                    // If it's inside the region of the physical data instance
                    if (xi >= instanceRegion->x() &&
                        xi <= (instanceRegion->x() + instanceRegion->dx()) &&
                        yi >= instanceRegion->y() &&
                        yi <= (instanceRegion->y() + instanceRegion->dy()) &&
                        zi >= instanceRegion->z() &&
                        zi <= (instanceRegion->z() + instanceRegion->dz())) {
                      scalar_t particleBuffer[5];
                      particleBuffer[0] = particle.x;
                      particleBuffer[1] = particle.y;
                      particleBuffer[2] = particle.z;
                      particleBuffer[3] = particles->radius(i);
                      particleBuffer[4] = particles->quantized_collision_distance(i);
                      data->AddToTempBuffer(reinterpret_cast<char*>(particleBuffer),
                                            sizeof(particleBuffer));
                    }
                  }
                }
              }
              particles = particles->next;
            }
          }
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

        Int3Vector min = (sa->domain).Minimum_Corner();
        Int3Vector delta = (sa->domain).Edge_Lengths() + 1;
        GeometricRegion dest_region(min.x, min.y, min.z, delta.x, delta.y, delta.z);

        if (instances != NULL) {
            PdiVector::const_iterator iter = instances->begin();
            for (; iter != instances->end(); ++iter) {
                const PhysicalDataInstance* inst = *iter;
                Dimension3Vector overlap = GetOverlapSize(inst->region(), region);

                if (HasOverlap(overlap)) {
                    dbg(DBG_TRANSLATE, "Incorporating physical object %lu into FaceArray.\n",
                            inst->id());
                    PhysBAMData* data = static_cast<PhysBAMData*>(inst->data());
                    scalar_t* buffer  = reinterpret_cast<scalar_t*>(data->buffer());

                    Dimension3Vector src  = GetOffset(inst->region(), region);
                    Dimension3Vector dest = GetOffset(&dest_region, inst->region());

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
                                int dest_x = x + dest(X_COORD) + 1;
                                int dest_y = y + dest(Y_COORD) + 1;
                                int dest_z = z + dest(Z_COORD) + 1;
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
        }

        Int3Vector min = (sa->domain).Minimum_Corner();
        Int3Vector delta = (sa->domain).Edge_Lengths() + 1;
        GeometricRegion src_region(min.x, min.y, min.z, delta.x, delta.y, delta.z);

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

                    Dimension3Vector src  = GetOffset(&src_region, inst->region());
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
                                int source_x = x + src(X_COORD) + 1;
                                int source_y = y + src(Y_COORD) + 1;
                                int source_z = z + src(Z_COORD) + 1;
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
