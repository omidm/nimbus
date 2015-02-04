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
  * The data structure "unordered_map_view" is used to accelerate table lookup
  * which is the most common operation in a Nimbus scheduler.
  * The idea is to create a static view of the hash table which is optimized for
  * lookup.
  *
  * Author: Hang Qu <quhang@stanford.edu>
  */

#ifndef NIMBUS_SHARED_UNORDERED_MAP_VIEW_H_
#define NIMBUS_SHARED_UNORDERED_MAP_VIEW_H_

#include <boost/unordered_map.hpp>

#include <algorithm>
#include <cassert>
#include <vector>

namespace boost {

template <class K, class T>
class unordered_map_view {
 public:
  typedef K key_type;
  typedef T mapped_type;
  typedef T* iterator;

  typedef std::size_t size_type;

  explicit unordered_map_view() {
    is_initialized_ = false;
  }
  ~unordered_map_view() {}

  template<class H, class P, class A>
  bool snapshot(const unordered_map<K, T, H, P, A>& unordered_map) {
    typedef typename boost::unordered_map<K, T, H, P, A>::const_iterator
        const_iterator;
    if (unordered_map.empty()) {
      is_initialized_ = false;
      return false;
    }
    const_iterator iter = unordered_map.cbegin();
    down_limit_ = iter->first;
    up_limit_ = iter->first;
    ++iter;
    for (; iter != unordered_map.cend(); ++iter) {
      if (iter->first > up_limit_) up_limit_ = iter->first;
      if (iter->first < down_limit_) down_limit_ = iter->first;
    }
    size_ = unordered_map.size();
    is_continuous_ =
        (static_cast<key_type>(size_) == (up_limit_ - down_limit_ + 1));
    if (up_limit_ - down_limit_ + 1 >
        max_factor * static_cast<key_type>(size_)) {
      is_initialized_ = false;
      return false;
    }
    table_.clear();
    table_.resize(up_limit_ - down_limit_ + 1);
    if (!is_continuous_) {
      exist_.resize(up_limit_ - down_limit_ + 1);
      std::fill(exist_.begin(), exist_.end(), false);
    }
    iter = unordered_map.cbegin();
    for (; iter != unordered_map.cend(); ++iter) {
      table_[iter->first - down_limit_] = iter->second;
      if (!is_continuous_) {
        exist_[iter->first - down_limit_] = iter->second;
      }
    }
    is_initialized_ = true;
    return true;
  }

  bool empty() const {
    assert(is_initialized_);
    return size_ == 0;
  }

  size_type size() const {
    assert(is_initialized_);
    return size_;
  }

  mapped_type const& operator[](const key_type& key) const {
    assert(is_initialized_ && key >= down_limit_ && key <= up_limit_);
    if (!is_continuous_) {
      assert(exist_[key - down_limit_]);
    }
    return table_[key - down_limit_];
  }
  mapped_type const& at(const key_type& key) const {
    assert(is_initialized_ && key >= down_limit_ && key <= up_limit_);
    if (!is_continuous_) {
      assert(exist_[key - down_limit_]);
    }
    return table_[key - down_limit_];
  }
  iterator find(const key_type& key) {
    assert(is_initialized_);
    if (key < down_limit_ && key > up_limit_) {
      return end();
    }
    if ((!is_continuous_) && (!exist_[key - down_limit_])) {
      return end();
    }
    return &(table_[key - down_limit_]);
  }
  iterator end() const {
    return NULL;
  }

  void clear() {
    table_.clear();
    exist_.clear();
    is_initialized_ = false;
  }

  bool has_element(const key_type& key) {
    assert(is_initialized_);
    if (key < down_limit_ || key > up_limit_) {
      return false;
    }
    if (!is_continuous_) {
      return exist_[key - down_limit_];
    } else {
      return true;
    }
  }

 private:
  static const int max_factor = 100;
  std::vector<T> table_;
  std::vector<bool> exist_;
  key_type down_limit_;
  key_type up_limit_;
  std::size_t size_;
  bool is_continuous_;
  bool is_initialized_;
};

}  // namespace boost

#endif  // NIMBUS_SHARED_UNORDERED_MAP_VIEW_H_
