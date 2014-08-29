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

#include <algorithm>
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
  exists_rtree_.clear();
  index_rtree_.clear();
}

LdoIndexRtree::Box LdoIndexRtree::RegionToBox(const GeometricRegion& region) {
  return Box(Point(region.x(), region.y(), region.z()),
             Point(region.x() + region.dx(),
                   region.y() + region.dy(),
                   region.z() + region.dz()));
}

bool LdoIndexRtree::AddObject(LogicalDataObject *object) {
#ifdef LDO_REFER
  refer.AddObject(new LogicalDataObject(*object));
#endif  // LDO_REFER
  // If the object already exists, reject it. Otherwise, insert
  // it into a list, creating a list if needed.

  // TODO(quhang): cannot handle zero-area case.
  assert(object->region()->NoneZeroArea());

  dbg(DBG_DATA_OBJECTS, "LdoIndexRtree: adding %llu as %s.\n",
      object->id(), object->variable().c_str());

  if (exists_rtree_.find(object->id()) != exists_rtree_.end()) {
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
    exists_rtree_[object->id()] = object;
    return true;
  }
}

bool LdoIndexRtree::HasObject(logical_data_id_t id) {
  return (exists_rtree_.find(id) != exists_rtree_.end());
}

bool LdoIndexRtree::RemoveObject(logical_data_id_t id) {
#ifdef LDO_REFER
  refer.RemoveObject(id);
#endif  // LDO_REFER
  dbg(DBG_TEMP, "Trying to remove object %llu.\n", id);

  if (HasObject(id)) {
    LogicalDataObject* obj = exists_rtree_[id];
    assert(obj->id() == id);
    std::string variable = obj->variable();
    Rtree* rtree = index_rtree_[variable];
    rtree->remove(Value(RegionToBox(*obj->region()), id));
    int cnt = exists_rtree_.erase(id);
    dbg(DBG_TEMP, "Removing object %llu, removed %i elements from exists_.\n", id, cnt);
    return true;
  } else {
    return false;
  }
}

bool LdoIndexRtree::RemoveObject(LogicalDataObject* object) {
#ifdef LDO_REFER
  refer.RemoveObject(object);
#endif  // LDO_REFER
  return RemoveObject(object->id());
}

LogicalDataObject* LdoIndexRtree::SpecificObject(logical_data_id_t id) {
  if (!HasObject(id)) {
    return NULL;
  } else {
    return exists_rtree_[id];
  }
}

int LdoIndexRtree::AllObjects(CLdoVector* dest) {
  dest->clear();
  int count = 0;
  LdoIdIndex::iterator it = exists_rtree_.begin();
  for (; it != exists_rtree_.end(); ++it) {
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
    assert(exists_rtree_.find(iter->second) != exists_rtree_.end());
    dest->push_back(exists_rtree_[iter->second]);
  }
}

void LdoIndexRtree::Sort(CLdoVector* dest) {
  assert(dest != NULL);
  for (uint32_t i = 0; i < dest->size(); ++i) {
    for (uint32_t j = i + 1; j < dest->size(); ++j) {
      if ((*dest)[i]->id() > (*dest)[j]->id()) {
        std::swap((*dest)[i], (*dest)[j]);
      }
    }
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
  assert(dest != NULL);
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  int count = 0;
  for (Rtree::const_query_iterator iter =
       rtree->qbegin(boost::geometry::index::intersects(RegionToBox(*region)));
       iter != rtree->qend();
       ++iter) {
    LogicalDataObject* object = exists_rtree_[iter->second];
    assert(object != NULL);
    if (region->Intersects(object->region())) {
      dest->push_back(object);
      ++count;
    }
  }
  // Sort(dest);
#ifdef LDO_REFER
  CLdoVector refer_dest;
  assert(count == refer.IntersectingObjects(variable, region, &refer_dest));
#endif  // LDO_REFER
  return count;
}

int LdoIndexRtree::CoveredObjects(const std::string& variable,
                                  const GeometricRegion* region,
                                  CLdoVector* dest) {
  assert(dest != NULL);
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  int count = 0;
  for (Rtree::const_query_iterator iter =
       rtree->qbegin(boost::geometry::index::covered_by(RegionToBox(*region)));
       iter != rtree->qend();
       ++iter) {
    LogicalDataObject* object = exists_rtree_[iter->second];
    assert(object != NULL);
    assert(region->Covers(object->region()));
    dest->push_back(object);
    ++count;
  }
  // Sort(dest);
#ifdef LDO_REFER
  CLdoVector refer_dest;
  assert(count == refer.CoveredObjects(variable, region, &refer_dest));
#endif  // LDO_REFER
  return count;
}

int LdoIndexRtree::AdjacentObjects(const std::string& variable,
                                   const GeometricRegion* region,
                                   CLdoVector* dest) {
  assert(dest != NULL);
  dest->clear();
  if (index_rtree_.find(variable) == index_rtree_.end()) {  // No such variable
    return 0;
  }

  Rtree* rtree = index_rtree_[variable];
  int count = 0;
  for (Rtree::const_query_iterator iter =
       rtree->qbegin(boost::geometry::index::intersects(RegionToBox(*region)));
       iter != rtree->qend();
       ++iter) {
    LogicalDataObject* object = exists_rtree_[iter->second];
    assert(object != NULL);
    assert(region->Adjacent(object->region()));
    dest->push_back(object);
    ++count;
  }
  // Sort(dest);
#ifdef LDO_REFER
  CLdoVector refer_dest;
  assert(count == refer.AdjacentObjects(variable, region, &refer_dest));
#endif  // LDO_REFER
  return count;
}

}  // namespace nimbus
