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
  * The graph graph abstraction class, build as template class so that one can
  * choose the data structure to keep the data at each vertex dynamically.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_PRIORITY_QUEUE_H_
#define NIMBUS_SHARED_PRIORITY_QUEUE_H_

#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <list>
#include <set>
#include "shared/dbg.h"

namespace nimbus {

enum PRIORITY_TYPE {
  LOW_PRIORITY,
  HIGH_PRIORITY
};

template<typename T>
class PriorityQueue {
  public:
    typedef typename std::list<T> Queue;
    typedef typename std::list<T>::reference reference;
    typedef typename std::list<T>::const_reference const_reference;

    PriorityQueue();
    virtual ~PriorityQueue();

    bool empty() const;
    reference front();
    const_reference front() const;
    void pop_front();
    void push_back(const T& val);
    void push_back(const T& val, PRIORITY_TYPE priority);

  private:
    Queue l_queue_;
    Queue h_queue_;
};


template<typename T>
PriorityQueue<T>::PriorityQueue() {
}

template<typename T>
PriorityQueue<T>::~PriorityQueue() {
}

template<typename T>
bool PriorityQueue<T>::empty() const {
  return (l_queue_.empty() && h_queue_.empty());
}

template<typename T>
typename PriorityQueue<T>::reference PriorityQueue<T>::front() {
  if (!h_queue_.empty()) {
    return h_queue_.front();
  } else {
    return l_queue_.front();
  }
}

template<typename T>
typename PriorityQueue<T>::const_reference PriorityQueue<T>::front() const {
  if (!h_queue_.empty()) {
    return h_queue_.front();
  } else {
    return l_queue_.front();
  }
}

template<typename T>
void PriorityQueue<T>::pop_front() {
  if (!h_queue_.empty()) {
    h_queue_.pop_front();
  } else {
    l_queue_.pop_front();
  }
}

template<typename T>
void PriorityQueue<T>::push_back(const T& val) {
    h_queue_.push_back(val);
}

template<typename T>
void PriorityQueue<T>::push_back(const T& val, PRIORITY_TYPE priority) {
  switch (priority) {
    case HIGH_PRIORITY:
      h_queue_.push_back(val);
      break;
    case LOW_PRIORITY:
      l_queue_.push_back(val);
      break;
    default:
      assert(false);
  }
}

}  // namespace nimbus

#endif  // NIMBUS_SHARED_PRIORITY_QUEUE_H_
