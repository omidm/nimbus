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
  * Unstructured region which is a set of non-overlapping geometric regions.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "scheduler/unstructured_region.h"

namespace nimbus {

UnstructuredRegion::UnstructuredRegion() {
}

UnstructuredRegion::UnstructuredRegion(const UnstructuredRegion& other) {
  region_list_ = other.region_list_;
}

UnstructuredRegion::~UnstructuredRegion() {
}

UnstructuredRegion& UnstructuredRegion::operator= (const UnstructuredRegion& right) {
  region_list_ = right.region_list_;
  return *this;
}

void UnstructuredRegion::AddRegion(const GeometricRegion *region) {
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

void UnstructuredRegion::RemoveRegion(const GeometricRegion *region) {
  RegionList result;

  RegionListIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
      RemoveIntersect(&(*iter), region, &result, true);
  }

  region_list_ = result;
}


void UnstructuredRegion::Grow(const UnstructuredRegion *u_region) {
  RegionListConstIter iter = u_region->region_list_.begin();
  for (; iter != u_region->region_list_.end(); ++iter) {
    AddRegion(&(*iter));
  }
}


int_dimension_t UnstructuredRegion::CommonSurface(const GeometricRegion *region) {
  int_dimension_t result = 0;
  RegionListIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    result += GeometricRegion::GetIntersection(*iter, *region).GetSurfaceArea();
  }

  return result;
}

std::string UnstructuredRegion::Print() {
  std::string str;
  RegionListIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    str += iter->ToNetworkData();
    str += "\n";
  }
  return str;
}

void UnstructuredRegion::RemoveIntersect(const GeometricRegion *o,
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

      GeometricRegion n(o->x(),
                        r->y(),
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

      GeometricRegion n(o->x(),
                        o->y(),
                        r->z(),
                        o->dx(),
                        o->dy(),
                        o->z() + o->dz() - r->z());
      RemoveIntersect(&n, r, result, true);
      return;
    } else if ((r->z() + r->dz()) < (o->z() + o->dz())) {
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


bool UnstructuredRegion::Adjacent(const UnstructuredRegion *u_region) const {
  return (AdjacentOrIntersects(u_region) && !Intersects(u_region));
}

bool UnstructuredRegion::Intersects(const UnstructuredRegion *u_region) const {
  RegionListConstIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    if (u_region->Intersects(&(*iter))) {
      return true;
    }
  }

  return false;
}

bool UnstructuredRegion::AdjacentOrIntersects(const UnstructuredRegion *u_region) const {
  RegionListConstIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    if (u_region->AdjacentOrIntersects(&(*iter))) {
      return true;
    }
  }

  return false;
}

bool UnstructuredRegion::Adjacent(const GeometricRegion *region) const {
  return (AdjacentOrIntersects(region) && !Intersects(region));
}

bool UnstructuredRegion::Intersects(const GeometricRegion *region) const {
  RegionListConstIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    if (region->Intersects(&(*iter))) {
      return true;
    }
  }

  return false;
}

bool UnstructuredRegion::AdjacentOrIntersects(const GeometricRegion *region) const {
  RegionListConstIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    if (region->AdjacentOrIntersects(&(*iter))) {
      return true;
    }
  }

  return false;
}

int_dimension_t UnstructuredRegion::GetDistance(const UnstructuredRegion *u_region) const {
  if (region_list_.size() == 0) {
    return 0;
  }

  std::vector<int_dimension_t> dist;

  RegionListConstIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    dist.push_back(u_region->GetDistance(&(*iter)));
  }

  std::sort(dist.begin(), dist.end());
  return *dist.begin();
}

int_dimension_t UnstructuredRegion::GetDistance(const GeometricRegion *region) const {
  if (region_list_.size() == 0) {
    return 0;
  }

  std::vector<int_dimension_t> dist;

  RegionListConstIter iter = region_list_.begin();
  for (; iter != region_list_.end(); ++iter) {
    dist.push_back(GeometricRegion::GetDistance(&(*iter), region));
  }

  std::sort(dist.begin(), dist.end());
  return *dist.begin();
}

}  // namespace nimbus
