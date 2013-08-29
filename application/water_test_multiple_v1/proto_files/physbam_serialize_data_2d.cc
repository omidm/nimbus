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

#include "assert.h"
#include "physbam_serialize_data_2d.h"

namespace physbam_pb {

    void make_pb_object(VI2 *phys_vec,
            ::communication::PhysbamVectorInt2d *pb_vec) {
        assert(pb_vec);
        if (phys_vec != NULL) {
            pb_vec->set_x((*phys_vec)[0]);
            pb_vec->set_y((*phys_vec)[1]);
        }
    }

    void make_pb_object(VF2 *phys_vec,
            ::communication::PhysbamVectorFloat2d *pb_vec) {
        assert(pb_vec);
        if (phys_vec != NULL) {
            pb_vec->set_x((*phys_vec)[0]);
            pb_vec->set_y((*phys_vec)[1]);
        }
    }

    void make_pb_object(RangeI2 *phys_range,
            ::communication::PhysbamRangeInt2d *pb_range) {
        assert(pb_range);
        if (phys_range != NULL) {
            VI2 phys_range_min = phys_range->Minimum_Corner();
            VI2 phys_range_max = phys_range->Maximum_Corner();
            ::communication::PhysbamVectorInt2d *pb_range_min = 
                pb_range->mutable_min_corner();
            ::communication::PhysbamVectorInt2d *pb_range_max = 
                pb_range->mutable_max_corner();
            make_pb_object(&phys_range_min, pb_range_min);
            make_pb_object(&phys_range_max, pb_range_max);
        }
    }

    void make_pb_object(RangeF2 *phys_range,
            ::communication::PhysbamRangeFloat2d *pb_range) {
        assert(pb_range);
        if (phys_range != NULL) {
            VF2 phys_range_min = phys_range->Minimum_Corner();
            VF2 phys_range_max = phys_range->Maximum_Corner();
            ::communication::PhysbamVectorFloat2d *pb_range_min = 
                pb_range->mutable_min_corner();
            ::communication::PhysbamVectorFloat2d *pb_range_max = 
                pb_range->mutable_max_corner();
            make_pb_object(&phys_range_min, pb_range_min);
            make_pb_object(&phys_range_max, pb_range_max);
        }
    }

    void make_pb_object(Grid2 *phys_grid,
            ::communication::PhysbamGrid2d *pb_grid) {

        assert(pb_grid);

        if (phys_grid != NULL) {

            VI2 phys_grid_counts = phys_grid->counts;
            ::communication::PhysbamVectorInt2d *pb_grid_counts =
                pb_grid->mutable_counts();
            make_pb_object(&phys_grid_counts, pb_grid_counts);

            RangeF2 phys_grid_domain = phys_grid->domain;
            ::communication::PhysbamRangeFloat2d *pb_grid_domain =
                pb_grid->mutable_domain();
            make_pb_object(&phys_grid_domain, pb_grid_domain);

            VF2 phys_grid_dx = phys_grid->dX;
            ::communication::PhysbamVectorFloat2d *pb_grid_dx =
                pb_grid->mutable_dx();
            make_pb_object(&phys_grid_dx, pb_grid_dx);

            VF2 phys_grid_one_over_dx = phys_grid->one_over_dX;
            ::communication::PhysbamVectorFloat2d *pb_grid_one_over_dx =
                pb_grid->mutable_one_over_dx();
            make_pb_object(&phys_grid_one_over_dx, pb_grid_one_over_dx);

            VI2 phys_grid_num_cells = phys_grid->numbers_of_cells;
            ::communication::PhysbamVectorInt2d *pb_grid_num_cells =
                pb_grid->mutable_numbers_of_cells();
            make_pb_object(&phys_grid_num_cells, pb_grid_num_cells);

            pb_grid->set_min_dx(phys_grid->min_dX);
            pb_grid->set_mac_offset(phys_grid->MAC_offset);
        }
    }

} // namespace physbam_pb
