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
 * Author: Andrew Lim <alim16@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_SMOKE_SMOKE_APP_H_
#define NIMBUS_APPLICATION_SMOKE_SMOKE_APP_H_

#include "src/shared/log.h"
#include "src/shared/nimbus.h"

namespace application {

    /* Smoke application lunached by a nimbus worker. */
    class SmokeApp : public nimbus::Application {
        public:
            Log *translator_log;
            SmokeApp();
            virtual void Load();

            void set_scale(uint64_t scale);
            void set_part_num_x(uint64_t part_num_x);
            void set_part_num_y(uint64_t part_num_y);
            void set_part_num_z(uint64_t part_num_z);
            void set_projection_part_num_x(uint64_t projection_part_num_x);
            void set_projection_part_num_y(uint64_t projection_part_num_y);
            void set_projection_part_num_z(uint64_t projection_part_num_z);
            void set_last_frame(uint64_t last_frame);
            void set_max_iterations(uint64_t max_iterations);
            void set_global_write(bool flag);

            uint64_t projection_part_num_x() {return projection_part_num_x_;}
            uint64_t projection_part_num_y() {return projection_part_num_y_;}
            uint64_t projection_part_num_z() {return projection_part_num_z_;}

        private:
            uint64_t scale_;
            uint64_t part_num_x_;
            uint64_t part_num_y_;
            uint64_t part_num_z_;
            uint64_t projection_part_num_x_;
            uint64_t projection_part_num_y_;
            uint64_t projection_part_num_z_;
            uint64_t last_frame_;
            uint64_t max_iterations_;
            bool global_write_;
    };

} // namespace application

#endif  // NIMBUS_APPLICATION_SMOKE_SMOKE_APP_H_
