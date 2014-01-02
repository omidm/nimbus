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
  *
  * Work on particle-related functionality.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_TEST_H_
#define NIMBUS_DATA_PHYSBAM_TRANSLATOR_PHYSBAM_TEST_H_

#include <algorithm>
#include <cmath>

#include "data/physbam/physbam_include.h"
#include "data/physbam/physbam_data.h"

#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "worker/physical_data_instance.h"
#include "worker/worker.h"

namespace nimbus {

template <class VECTOR_TYPE> class TranslatorPhysBAMTest {
 public:
  typedef VECTOR_TYPE TV;
  typedef typename TV::SCALAR scalar_t;
  typedef typename TV::template REBIND<int>::TYPE TV_INT;
  typedef PhysBAM::GRID<TV> Grid;

  // TODO(quhang) Document particles-related classes well and give them better
  // names.
  typedef typename PhysBAM::PARTICLE_LEVELSET<Grid> ParticleProfile;
  typedef typename PhysBAM::PARTICLE_LEVELSET_UNIFORM<Grid> ParticleContainer;
  typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> ParticlesUnit;
  typedef typename PhysBAM::ARRAY<ParticlesUnit*, TV_INT> ParticlesArray;

  bool ReadParticles(GeometricRegion* region,
                     CPdiVector* instances,
                     ParticleContainer& particle_container,
                     bool positive) {
    // TODO(quhang) Check whether particle_container has valid particle data
    // before moving on.
    ParticlesArray* particles;
    if (positive) {
      particles = &particle_container.positive_particles;
    } else {
      particles = &particle_container.negative_particles;
    }

    // Allocates buckets.
    for (int z = 1; z <= region->dz(); z++)
      for (int y = 1; y <= region->dy(); y++)
        for (int x = 1; x <= region->dx(); x++) {
          TV_INT block_index(x, y, z);
          particle_container.Free_Particle_And_Clear_Pointer(
              (*particles)(block_index));
          if (!(*particles)(block_index)) {
            (*particles)(block_index) = particle_container.Allocate_Particles(
                particle_container.template_particles);
          }
        }
    dbg(DBG_WARN, "Successfully clean particle data structures!\n");

    // TODO(quhang) warning rather than assert.
    assert(instances != NULL);

    CPdiVector::iterator iter = instances->begin();
    for (; iter != instances->end(); ++iter) {
      const PhysicalDataInstance* instance = *iter;
      PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
      scalar_t* buffer = reinterpret_cast<scalar_t*>(data->buffer());

      // @Phil. The data size is returned in the unit of bytes.  --quhang
      // TODO(anyone) anyone knows how to avoid the ugly casting?
      for (int i = 0;
           i < static_cast<int>(data->size())
           / static_cast<int>(sizeof(float));  // NOLINT
           i+= 5) {
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
        // TODO(quhang): Check whether the cast is safe.
        int_dimension_t xi = (int_dimension_t)floor(x - region->x() + 1);
        int_dimension_t yi = (int_dimension_t)floor(y - region->y() + 1);
        int_dimension_t zi = (int_dimension_t)floor(z - region->z() + 1);

        // TODO(quhang): The condition is not accurate.
        // If particle is within region, copy it to particles
        if (xi >= 1 &&
            xi <= region->dx() &&
            yi >= 1 &&
            yi <= region->dy() &&
            zi >= 1 &&
            zi <= region->dz()) {
          ParticlesUnit* cellParticles = (*particles)(TV_INT(xi, yi, zi));

          // Note that Add_Particle traverses a linked list of particle
          // buckets, so it's O(N^2) time. Blech.
          int index = particle_container.Add_Particle(cellParticles);
          cellParticles->quantized_collision_distance(index) =
            collision_distance;
          cellParticles->X(index) = position;
          cellParticles->radius(index) = radius;
        }
      }
    }
    dbg(DBG_WARN, "Successfully finished ReadParticles!\n");
    return true;
  }

  bool WriteParticles(GeometricRegion* region,
                      CPdiVector* instances,
                      ParticleContainer& particle_container,
                      bool positive) {
    CPdiVector::iterator iter = instances->begin();
    for (; iter != instances->end(); ++iter) {
      const PhysicalDataInstance* instance = *iter;
      PhysBAMData* data = static_cast<PhysBAMData*>(instance->data());
      data->ClearTempBuffer();
    }

    for (int z = 1; z <= region->dz(); z++) {
      for (int y = 1; y <= region->dy(); y++) {
        for (int x = 1; x <= region->dx(); x++) {
          ParticlesArray* arrayPtr;
          if (positive) {
            arrayPtr = &particle_container.positive_particles;
          } else {
            arrayPtr = &particle_container.negative_particles;
          }
          ParticlesUnit* particles = (*arrayPtr)(TV_INT(x, y, z));
          while (particles) {
            for (int i = 1; i <= particles->array_collection->Size(); i++) {
              VECTOR_TYPE particle = particles->X(i);
              scalar_t x = particle.x;
              scalar_t y = particle.y;
              scalar_t z = particle.z;
              double xi = x;
              double yi = y;
              double zi = z;
              // TODO(quhang): I am almost 100% sure this is wrong. The
              // condition is not accurate. But needs time to figure out.
              // If it's inside the region,
              if (xi >= region->x() &&
                  xi < (region->x() + region->dx()) &&
                  yi >= region->y() &&
                  yi < (region->y() + region->dy()) &&
                  zi >= region->z() &&
                  zi < (region->z() + region->dz())) {
                CPdiVector::iterator iter = instances->begin();
                // Iterate across instances, checking each one
                for (; iter != instances->end(); ++iter) {
                  const PhysicalDataInstance* instance = *iter;
                  GeometricRegion* instanceRegion = instance->region();
                  PhysBAMData* data =
                      static_cast<PhysBAMData*>(instance->data());

                  // If it's inside the region of the physical data instance
                  if (xi >= instanceRegion->x() &&
                      xi < (instanceRegion->x() + instanceRegion->dx()) &&
                      yi >= instanceRegion->y() &&
                      yi < (instanceRegion->y() + instanceRegion->dy()) &&
                      zi >= instanceRegion->z() &&
                      zi < (instanceRegion->z() + instanceRegion->dz())) {
                    scalar_t particleBuffer[5];
                    particleBuffer[0] = particle.x;
                    particleBuffer[1] = particle.y;
                    particleBuffer[2] = particle.z;
                    particleBuffer[3] = particles->radius(i);
                    particleBuffer[4] =
                        particles->quantized_collision_distance(i);
                    data->AddToTempBuffer(
                        reinterpret_cast<char*>(particleBuffer),
                        sizeof(scalar_t)*5);
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
};

}  // namespace nimbus

#endif
