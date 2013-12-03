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
  * The pair elements in the list are copied around whenever this object
  * is passed around to guarantee correctness. Very low efficiency.
  *
  * The plus operator might be a bad idea, needs more thinking.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include <algorithm>
#include "test/experimental_quh/variable_region_set.h"

namespace nimbus {

VariableRegionSet::VariableRegionSet() {
  Init();
}

VariableRegionSet::VariableRegionSet(
    const std::string& variable, const GeometricRegion& region) {
  Init();
  _var_region_list.push_back(VarRegion(variable, region));
}

VariableRegionSet::VariableRegionSet(const VariableRegionSet& vrs) {
  Init();
  std::copy(vrs._var_region_list.begin(), vrs._var_region_list.end(),
      std::back_inserter(_var_region_list));
}

VariableRegionSet& VariableRegionSet::CopyFrom(const VariableRegionSet& vrs) {
  Init();
  std::copy(vrs._var_region_list.begin(), vrs._var_region_list.end(),
      std::back_inserter(_var_region_list));
  return *this;
}

VariableRegionSet& VariableRegionSet::operator+=(
    const VariableRegionSet right_vrs) {
  std::copy(
      right_vrs._var_region_list.begin(),
      right_vrs._var_region_list.end(),
      std::back_inserter(_var_region_list));
  return *this;
}

void VariableRegionSet::DebugPrint() {
  for (VarRegionList::iterator index = _var_region_list.begin();
       index != _var_region_list.end();
       index++) {
    // TODO(quhang) Use correct log functionality.
    std::cout << index->first << ':' << index->second.toString() << std::endl;
  }
}

void VariableRegionSet::Init() {
  _var_region_list.clear();
}

bool VariableRegionSet::IntersectsAndDelete(const VariableRegionSet& vrs) {
  return IntersectsImpl(vrs, true);
}

bool VariableRegionSet::IntersectsTest(const VariableRegionSet& vrs) {
  return IntersectsImpl(vrs, false);
}

void VariableRegionSet::AddRegionHelper(
    const std::string& variable,
    const GeometricRegion& origin,
    const GeometricRegion& target) {
  // TODO(quhang) Not the best place for impl here.
  std::vector<int_dimension_t> x_seg, y_seg, z_seg;
  PushSegXHelper(origin, &x_seg);
  PushSegXHelper(target, &x_seg);
  PushSegYHelper(origin, &y_seg);
  PushSegYHelper(target, &y_seg);
  PushSegZHelper(origin, &z_seg);
  PushSegZHelper(target, &z_seg);
  SortAndDeduplication(&x_seg);
  SortAndDeduplication(&y_seg);
  SortAndDeduplication(&z_seg);
  for (std::size_t i = 0; i < x_seg.size() - 1; ++i)
    for (std::size_t j = 0; j < y_seg.size() - 1; ++j)
      for (std::size_t k = 0; k < z_seg.size() - 1; ++k) {
        GeometricRegion temp(x_seg[i], y_seg[j], z_seg[k],
            x_seg[i+1]-x_seg[i], y_seg[j+1]-y_seg[j], z_seg[k+1]-z_seg[k]);
        if ((const_cast<GeometricRegion*> (&origin))->Covers(&temp)
            && !(const_cast<GeometricRegion*> (&target))->Covers(&temp)) {
          _var_region_list.push_back(VarRegion(variable, temp));
        }
      }
}

void VariableRegionSet::PushSegXHelper(
    const GeometricRegion& gr,
    std::vector<int_dimension_t>* vec) {
  vec->push_back(gr.x());
  vec->push_back(gr.x()+gr.dx());
}

void VariableRegionSet::PushSegYHelper(
    const GeometricRegion& gr,
    std::vector<int_dimension_t>* vec) {
  vec->push_back(gr.y());
  vec->push_back(gr.y()+gr.dy());
}

void VariableRegionSet::PushSegZHelper(
    const GeometricRegion& gr,
    std::vector<int_dimension_t>* vec) {
  vec->push_back(gr.z());
  vec->push_back(gr.z()+gr.dz());
}

void VariableRegionSet::SortAndDeduplication(
    std::vector<int_dimension_t>* vec) {
  std::sort(vec->begin(), vec->end());
  vec->erase(std::unique(vec->begin(), vec->end()), vec->end());
}

bool VariableRegionSet::IntersectsImpl(
    const VariableRegionSet& vrs, bool remove_mode) {
  bool result = false;
  VarRegionList::iterator origin = _var_region_list.begin();
  while (origin != _var_region_list.end()) {
    for (VarRegionList::iterator target
            = (const_cast<VariableRegionSet*> (&vrs))->_var_region_list.begin();
         target != vrs._var_region_list.end();
         target++)
      if (origin->first == target->first
          && origin->second.Intersects(&target->second)) {
        result = true;
        if (!remove_mode) return result;
        if (!target->second.Covers(&origin->second)) {
          AddRegionHelper(origin->first, origin->second, target->second);
        }
        VarRegionList::iterator temp = origin;
        origin++;
        _var_region_list.erase(temp);
        continue;
      }  // end if
    origin++;
  }  // end while
  return result;
}

VariableRegionSet operator+(
    VariableRegionSet left_vrs,
    const VariableRegionSet& right_vrs) {
  std::copy(
      right_vrs._var_region_list.begin(),
      right_vrs._var_region_list.end(),
      std::back_inserter(left_vrs._var_region_list));
  return left_vrs;
}
}  // namespace nimbus
