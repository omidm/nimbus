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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_multiple/data_particle_array.h"
#include "data/physbam/physbam_data.h"
#include "shared/nimbus.h"
#include "string.h"

namespace application {

DataParticleArray::DataParticleArray(std::string name) {
    set_name(name);
}

nimbus::Data* DataParticleArray::Clone() {
    return (new DataParticleArray(name()));
}

void DataParticleArray::Create() {
    set_size(0);
    nimbus::PhysBAMData::Create();
}

void DataParticleArray::MergeParticles(const std::vector<nimbus::Data *> &scratch) {
    nimbus::int_dimension_t new_size = 0;
    for (size_t i = 0; i < scratch.size(); i++) {
        DataParticleArray *s = static_cast<DataParticleArray *>(scratch[i]);
        new_size += s->size();
    }
    char *new_buffer = new char[new_size];
    char *buf_merged = new_buffer;
    for (size_t i = 0; i < scratch.size(); i++) {
        DataParticleArray *s = static_cast<DataParticleArray *>(scratch[i]);
        if (s->size() == 0)
            continue;
        memcpy(buf_merged, s->buffer(), s->size());
        buf_merged += s->size();
    }
    delete buffer();
    set_buffer(new_buffer, new_size);
}

} // namespace application
