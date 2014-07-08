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
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <vector>
#include "shared/dbg.h"
#include "shared/ldo_index_rtree.h"

namespace nimbus {

LdoIndexRtree::LdoIndexRtree() {}

LdoIndexRtree::~LdoIndexRtree() {
  LdoVariableIndexRtree::iterator it = index_rtree_.begin();
  for (; it != index_rtree_.end(); ++it) {
    Rtree* rtree = (*it).second;
    dbg(DBG_MEMORY, "Deleting Rtree 0x%x\n", rtree);
    delete rtree;
  }
  exists_.clear();
  index_rtree_.clear();
}

LdoIndexRtree::Box LdoIndexRtree::RegionToBox(const GeometricRegion& region) {
  return Box(Point(region.x(), region.y(), region.z()),
             Point(region.x() + region.dx(),
                   region.y() + region.dy(),
                   region.z() + region.dx()));
}

bool LdoIndexRtree::AddObject(LogicalDataObject *object) {
  // If the object already exists, reject it. Otherwise, insert
  // it into a list, creating a list if needed.

  // TODO(quhang): cannot handle zero-area case.
  assert(object->region()->NoneZeroArea());

  dbg(DBG_DATA_OBJECTS, "LdoIndexRtree: adding %llu as %s.\n",
      object->id(), object->variable().c_str());

  if (exists_.find(object->id()) != exists_.end()) {
    dbg(DBG_DATA_OBJECTS|DBG_ERROR,
        "  - FAIL LdoIndexRtree: object %llu already in index.\n", object->id());
    return false;
  } else {
    std::string var = object->variable();
    Rtree* rtree;
    if (index_rtree_.find(var) == index_rtree_.end()) {
      rtree = new Rtree;
      dbg(DBG_MEMORY, "Allocating Rtree 0x%x\n", rtree);
      index_rtree_[var] = rtree;
    } else {
      rtree = index_rtree_[var];
    }
    rtree->insert(Value(RegionToBox(*object->region()), object->id()));
    exists_[object->id()] = object;
    return true;
  }
}

bool LdoIndexRtree::HasObject(logical_data_id_t id) {
  return (exists_.find(id) != exists_.end());
}

bool LdoIndexRtree::RemoveObject(logical_data_id_t id) {
  dbg(DBG_TEMP, "Trying to remove object %llu.\n", id);

  if (HasObject(id)) {
    LogicalDataObject* obj = exists_[id];
    assert(obj->id() == id);
    std::string variable = obj->variable();
    Rtree* rtree = index_rtree_[variable];
    rtree->remove(Value(RegionToBox(*obj->region()), id));
    int cnt = exists_.erase(id);
    dbg(DBG_TEMP, "Removing object %llu, removed %i elements from exists_.\n", id, cnt);
    return true;
  } else {
    return false;
  }
}

bool LdoIndexRtree::RemoveObject(LogicalDataObject* object) {
  return RemoveObject(object->id());
}

LogicalDataObject* LdoIndexRtree::SpecificObject(logical_data_id_t id) {
  if (!HasObject(id)) {
    return NULL;
  } else {
    return exists_[id];
  }
}

int LdoIndexRtree::AllObjects(CLdoVector* dest) {
  dest->clear();
  int count = 0;
  LdoIdIndex::iterator it = exists_.begin();
  for (; it != exists_.end(); ++it) {
    LogicalDataObject* obj = (*it).second;
    dest->push_back(obj);
    count++;
  }
  return count;
}

void LdoIndexRtree::Transform(const std::vector<Value>& result,
                              CLdoVector* dest) {
  for (std::vector<Value>::const_iterator iter = result.begin();
       iter != result.end();
       ++iter) {
    assert(exists_.find(iter->second) != exists_.end());
    dest->push_back(exists_[iter->second]);
  }
}

int LdoIndexRtree::AllObjects(const std::string& variable,
                              CLdoVector* dest) {
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  std::vector<Value> result;
  rtree->query(boost::geometry::index::satisfies(AlwaysTrue),
               std::back_inserter(result));
  Transform(result, dest);
  return dest->size();
}

int LdoIndexRtree::IntersectingObjects(const std::string& variable,
                                       const GeometricRegion* region,
                                       CLdoVector* dest) {
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  std::vector<Value> result;
  rtree->query(boost::geometry::index::intersects(RegionToBox(*region)),
               std::back_inserter(result));
  Transform(result, dest);
  return dest->size();
}

int LdoIndexRtree::CoveredObjects(const std::string& variable,
                                  const GeometricRegion* region,
                                  CLdoVector* dest) {
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  std::vector<Value> result;
  rtree->query(boost::geometry::index::covered_by(RegionToBox(*region)),
               std::back_inserter(result));
  Transform(result, dest);
  return dest->size();
}

int LdoIndexRtree::AdjacentObjects(const std::string& variable,
                                   const GeometricRegion* region,
                                   CLdoVector* dest) {
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  std::vector<Value> result;
  // TODO(quhang): not sure if the two sets have overlap.
  rtree->query(boost::geometry::index::covered_by(RegionToBox(*region)),
               std::back_inserter(result));
  rtree->query(boost::geometry::index::overlaps(RegionToBox(*region)),
               std::back_inserter(result));
  Transform(result, dest);
  return dest->size();
}

}  // namespace nimbus
