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
  * A logical data object, which may have multiple physical copies throughout
  * the system. Provides methods for finding its physical copies and
  * its geometric regions. Every logical object has a variable associated
  * with it. So, for example, a logical object represents the storage of
  * velocity field of some subset of the simulation domain. The GeometricRegion
  * object describes a conservative geometric region that encloses the
  * associated data.
  *
  * Logical data objects can have one of three
  * basic structures: partitions, hierarchies, and overlaps. Within a
  * simulation, partition objects represent elements from a disjoint
  * partitioning of the simulation domain. Hierarchy objects represent
  * one node in a hierarchy of different scales, so they have parents
  * and potentially children. Overlap objects are free-moving data objects
  * that may overlap with one another.
  */

#ifndef NIMBUS_SHARED_LOGICAL_DATA_OBJECT_H_
#define NIMBUS_SHARED_LOGICAL_DATA_OBJECT_H_

#include <list>
#include <set>
#include <string>
#include <vector>
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"

namespace nimbus {

  class LogicalDataObject {
  public:
    LogicalDataObject();
    LogicalDataObject(logical_data_id_t id,
                      std::string variable,
                      GeometricRegion* region);

    LogicalDataObject(logical_data_id_t id,
                      std::string variable,
                      GeometricRegion* region,
                      partition_id_t partition);

    virtual ~LogicalDataObject();

    virtual logical_data_id_t id() const;
    virtual std::string variable() const;
    virtual GeometricRegion* region() const;
    virtual partition_id_t partition() const;

    virtual void FillInMessage(LdoMessage* mg);

    virtual bool Parse(std::istream* is);
    virtual bool Parse(const std::string& data);
    virtual bool Serialize(std::ostream* os);
    virtual bool SerializeToString(std::string* output);

  private:
    logical_data_id_t id_;
    GeometricRegion* region_;
    std::string variable_;
    partition_id_t partition_;
  };

  typedef std::set<LogicalDataObject*> LdoSet;
  typedef std::list<LogicalDataObject*> LdoList;
  typedef std::vector<LogicalDataObject*> LdoVector;
  typedef std::vector<const LogicalDataObject*> CLdoVector;

}  // namespace nimbus

#endif  // NIMBUS_SHARED_LOGICAL_DATA_OBJECT_H_
