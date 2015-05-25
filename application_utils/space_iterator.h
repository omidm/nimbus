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
  * Object used to traverse the space based on a bandwidth size.
  * Used to create core, face, vertex, and edge partitions.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_APPLICATION_UTILS_SPACE_ITERATOR_H_
#define NIMBUS_APPLICATION_UTILS_SPACE_ITERATOR_H_

#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <vector>
#include <map>
#include <set>
#include "shared/nimbus_types.h"



namespace nimbus {

class SpaceIterator {
 public:
  SpaceIterator(int size,
                int bw,
                size_t part_num,
                bool with_global_bw);

  ~SpaceIterator();

  struct Cursor {
    int_dimension_t point_;
    int_dimension_t delta_;
    enum Type {
      GHOST = 0,
      CORE  = 1
    };
    Type type_;
  };

  bool Initialize();

  bool Advance();

  void ReadCursor(Cursor *cursor);

 private:
  int bw_;
  int size_;
  size_t part_num_;
  bool wgbw_;
  bool initialized_;

  enum State {
    NEXT_GHOST,
    NEXT_CORE
  };

  State state_;
  Cursor cursor_;
};

}  // namespace nimbus

#endif  // NIMBUS_APPLICATION_UTILS_SPACE_ITERATOR_H_
