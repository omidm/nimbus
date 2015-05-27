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

#include <assert.h>
#include "application_utils/space_iterator.h"

using namespace nimbus; // NOLINT


SpaceIterator::SpaceIterator(int size,
                             int bw,
                             size_t part_num,
                             bool with_global_bw) {
  bw_ = bw;
  size_ = size;
  part_num_ = part_num;
  wgbw_ = with_global_bw;
  initialized_ = false;

  assert(bw_ >= 0);
  assert(size_ > 0);
  assert(part_num_ > 0);
  assert((size % part_num_) == 0);
  assert((size_ / part_num_) > (2 * bw_));
}

SpaceIterator::~SpaceIterator() {
}

bool SpaceIterator::Initialize() {
  if (bw_ > 0) {
    if (wgbw_) {
      cursor_.point_ = START_INDEX - bw_;
      cursor_.delta_ = bw_;
      cursor_.type_ = Cursor::GHOST;
      state_ = NEXT_GHOST;
    } else {
      cursor_.point_ = START_INDEX;
      cursor_.delta_ = bw_;
      cursor_.type_ = Cursor::GHOST;
      state_ = NEXT_CORE;
    }
  } else {
    assert(bw_ == 0);
    cursor_.point_ = START_INDEX;
    cursor_.delta_ = size_ / part_num_;
    cursor_.type_ = Cursor::CORE;
    state_ = NEXT_CORE;
  }

  initialized_ = true;
  return true;
}

bool SpaceIterator::Advance() {
  assert(initialized_);

  if (wgbw_) {
    if ((cursor_.point_ + cursor_.delta_) == (START_INDEX + size_ + bw_)) {
      return false;
    }
  } else {
    if ((cursor_.point_ + cursor_.delta_) == (START_INDEX + size_)) {
      return false;
    }
  }

  if (bw_ > 0) {
    switch (state_) {
      case NEXT_CORE:
        state_ = NEXT_GHOST;
        cursor_.point_ = cursor_.point_ + cursor_.delta_;
        cursor_.delta_ = (size_ / part_num_) - (2 * bw_);
        cursor_.type_ = Cursor::CORE;
        break;
      case NEXT_GHOST:
        if (cursor_.type_ == Cursor::GHOST) {
          state_ = NEXT_CORE;
        } else {
          state_ = NEXT_GHOST;
        }
        cursor_.point_ = cursor_.point_ + cursor_.delta_;
        cursor_.delta_ = bw_;
        cursor_.type_ = Cursor::GHOST;
        break;
      default:
        assert(false);
    }
  } else {
    assert(bw_ == 0);
    assert(state_ == NEXT_CORE);
    cursor_.point_ = cursor_.point_ + cursor_.delta_;
    cursor_.type_ = Cursor::CORE;
  }

  return true;
}

void SpaceIterator::ReadCursor(Cursor *cursor) {
  assert(initialized_);
  cursor->point_ = cursor_.point_;
  cursor->delta_ = cursor_.delta_;
  cursor->type_ = cursor_.type_;
  return;
}


