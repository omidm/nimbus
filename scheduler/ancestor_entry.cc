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
  * Job Ancestor Entry.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/ancestor_entry.h"

using namespace nimbus; // NOLINT

AncestorEntry::AncestorEntry(job_id_t id,
    boost::shared_ptr<VersionMap> version_map) {
  id_ = id;
  version_map_ = version_map;
}

AncestorEntry::AncestorEntry(const AncestorEntry& other) {
  id_ = other.id_;
  version_map_ = other.version_map_;
}

AncestorEntry::~AncestorEntry() {
}

job_id_t AncestorEntry::id() const {
  return id_;
}

boost::shared_ptr<VersionMap> AncestorEntry::version_map() const {
  return version_map_;
}

AncestorEntry& AncestorEntry::operator=(const AncestorEntry& right) {
  id_ = right.id_;
  version_map_ = right.version_map_;
  return (*this);
}


