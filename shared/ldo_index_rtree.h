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
  * An implementation of ldoindex using r-tree.
  * If macro LDO_REFER is defined, both the list implementation and the r-tree
  * implementation will be used and compared against each other.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_SHARED_LDO_INDEX_RTREE_H_
#define NIMBUS_SHARED_LDO_INDEX_RTREE_H_

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "shared/ldo_index.h"
#include "shared/logical_data_object.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"

#ifndef LDO_REFER
// #define LDO_REFER
#endif  // LDO_REFER

namespace nimbus {

class LdoIndexRtree : public LdoIndex {
 public:
  LdoIndexRtree();
  virtual ~LdoIndexRtree();

  virtual bool AddObject(LogicalDataObject* object);
  virtual bool HasObject(logical_data_id_t id);
  virtual bool RemoveObject(logical_data_id_t id);
  virtual bool RemoveObject(LogicalDataObject* object);

  virtual LogicalDataObject* SpecificObject(logical_data_id_t id);
  virtual int AllObjects(CLdoVector* dest);
  virtual int AllObjects(const std::string& variable,
                         CLdoVector* dest);
  virtual int IntersectingObjects(const std::string& variable,
                                  const GeometricRegion* region,
                                  CLdoVector* dest);
  virtual int CoveredObjects(const std::string& variable,
                             const GeometricRegion* region,
                             CLdoVector* dest);
  virtual int AdjacentObjects(const std::string& variable,
                              const GeometricRegion* region,
                              CLdoVector* dest);

 private:
  typedef boost::geometry::model::
      point<int_dimension_t, 3, boost::geometry::cs::cartesian> Point;
  typedef boost::geometry::model::box<Point> Box;
  typedef std::pair<Box, logical_data_id_t> Value;
  typedef boost::geometry::index::
      rtree<Value, boost::geometry::index::linear<16> > Rtree;
  typedef std::map<std::string, Rtree*> LdoVariableIndexRtree;
  LdoVariableIndexRtree index_rtree_;
  Box RegionToBox(const GeometricRegion& region);
  static bool AlwaysTrue(Value const& v) {
    return true;
  }
  void Transform(const std::vector<Value>& result, CLdoVector* dest);
  void Sort(CLdoVector* dest);
  // For existance query.
  // typedef std::map<logical_data_id_t, LogicalDataObject*> LdoIdIndex;
  LdoIdIndex exists_rtree_;

#ifdef LDO_REFER
  LdoIndex refer;
#endif  // LDO_REFER
};
}  // namespace nimbus

#endif  // NIMBUS_SHARED_LDO_INDEX_RTREE_H_
