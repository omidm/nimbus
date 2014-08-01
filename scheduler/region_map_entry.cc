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
  * Each element in the region map which defines an unstructured region.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "scheduler/region_map_entry.h"

namespace nimbus {

RegionMapEntry::RegionMapEntry() {
}

RegionMapEntry::RegionMapEntry(const RegionMapEntry& other) {
}

RegionMapEntry::~RegionMapEntry() {
}

RegionMapEntry& RegionMapEntry::operator= (const RegionMapEntry& right) {
  return *this;
}

void RegionMapEntry::AddRegion(const GeometricRegion *region) {
  RegionList regions_to_add;
  regions_to_add.push_back(*region);

  RegionListIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    RegionList result;
    RegionListIter it = regions_to_add.begin();
    for (; it != regions_to_add.end(); ++it) {
      RemoveIntersect(&(*it), &(*iter), &result, true);
    }
    regions_to_add = result;
  }

  iter = regions_to_add.begin();
  for (; iter != regions_to_add.end(); ++iter) {
    region_list_.push_back(*iter);
  }
}

void RegionMapEntry::RemoveRegion(const GeometricRegion *region) {
}

double RegionMapEntry::CommonSurface(const GeometricRegion *region) {
  return 0;
}


void RegionMapEntry::RemoveIntersect(const GeometricRegion *o,
                                     const GeometricRegion *r,
                                     RegionList *result,
                                     bool append) {
  if (!append) {
    result->clear();
  }

  if (!(o->Intersects(r))) {
    result->push_back(*o);
  } else {
    if (o->x() < r->x()) {
      GeometricRegion a(o->x(),
                        o->y(),
                        o->z(),
                        r->x() - o->x(),
                        o->dy(),
                        o->dz());
      result->push_back(a);

      GeometricRegion n(r->x(),
                        o->y(),
                        o->z(),
                        o->x() + o->dx() - r->x(),
                        o->dy(),
                        o->dz());
      RemoveIntersect(&n, r, result, true);
      return;
    } else if ((r->x() + r->dx()) < (o->x() + o->dx())) {
      GeometricRegion a(r->x() + r->dx(),
                        o->y(),
                        o->z(),
                        o->x() + o->dx() - r->x() - r->dx(),
                        o->dy(),
                        o->dz());
      result->push_back(a);

      GeometricRegion n(o->x(),
                        o->y(),
                        o->z(),
                        r->x() + r->dx() - o->x(),
                        o->dy(),
                        o->dz());
      RemoveIntersect(&n, r, result, true);
      return;
    }

    if (o->y() < r->y()) {
      GeometricRegion a(o->x(),
                        o->y(),
                        o->z(),
                        o->dx(),
                        r->y() - o->y(),
                        o->dz());
      result->push_back(a);

      GeometricRegion n(r->x(),
                        o->y(),
                        o->z(),
                        o->dx(),
                        o->y() + o->dy() - r->y(),
                        o->dz());
      RemoveIntersect(&n, r, result, true);
      return;
    } else if ((r->y() + r->dy()) < (o->y() + o->dy())) {
      GeometricRegion a(o->x(),
                        r->y() + r->dy(),
                        o->z(),
                        o->dx(),
                        o->y() + o->dy() - r->y() - r->dy(),
                        o->dz());
      result->push_back(a);

      GeometricRegion n(o->x(),
                        o->y(),
                        o->z(),
                        o->dx(),
                        r->y() + r->dy() - o->y(),
                        o->dz());
      RemoveIntersect(&n, r, result, true);
      return;
    }

    if (o->z() < r->z()) {
      GeometricRegion a(o->x(),
                        o->y(),
                        o->z(),
                        o->dx(),
                        o->dy(),
                        r->z() - o->z());
      result->push_back(a);

      GeometricRegion n(r->x(),
                        o->y(),
                        o->z(),
                        o->dx(),
                        o->dy(),
                        o->z() + o->dz() - r->z());
      RemoveIntersect(&n, r, result, true);
      return;
    } else if ((r->y() + r->dy()) < (o->y() + o->dy())) {
      GeometricRegion a(o->x(),
                        o->y(),
                        r->z() + r->dz(),
                        o->dx(),
                        o->dy(),
                        o->z() + o->dz() - r->z() - r->dz());
      result->push_back(a);

      GeometricRegion n(o->x(),
                        o->y(),
                        o->z(),
                        o->dx(),
                        o->dy(),
                        r->z() + r->dz() - o->z());
      RemoveIntersect(&n, r, result, true);
      return;
    }
  }
}

}  // namespace nimbus
