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
  * Maintains an index of all of the logical data objects (LDOs) in Nimbus.
  * Allows for quering which logical data objects cover what geometric
  * regions. Note that an LdoIndex does not manage memory for
  * LogicalDataObject parameters. A LogicalDataObject passed in
  * AddObject is assumed to be allocated and deleted externally, although
  * the LdoIndex maintains a reference. The DataManager class provides
  * for full-fledged memory management for a scheduler.
  */

#ifndef NIMBUS_SHARED_LDO_INDEX_H_
#define NIMBUS_SHARED_LDO_INDEX_H_

#include <unordered_map>
#include <map>
#include <string>
#include "shared/logical_data_object.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"

namespace nimbus {

  typedef std::unordered_map<std::string, LdoList*> LdoVariableIndex;
  typedef std::unordered_map<logical_data_id_t, LogicalDataObject*> LdoIdIndex;

  class LdoIndex {
  public:
    LdoIndex();
    virtual ~LdoIndex();

    virtual bool AddObject(LogicalDataObject* object);
    virtual bool HasObject(logical_data_id_t id);
    virtual bool RemoveObject(logical_data_id_t id);
    virtual bool RemoveObject(LogicalDataObject* object);

    virtual LogicalDataObject* SpecificObject(logical_data_id_t id);
    virtual int AllObjects(CLdoVector* dest);
    virtual int AllObjects(std::string variable,
                           CLdoVector* dest);
    virtual int IntersectingObjects(std::string variable,
                                    GeometricRegion* region,
                                    CLdoVector* dest);
    virtual int CoveredObjects(std::string variable,
                               GeometricRegion* region,
                               CLdoVector* dest);
    virtual int AdjacentObjects(std::string variable,
                                GeometricRegion* region,
                                CLdoVector* dest);

  private:
    LdoVariableIndex index_;
    LdoIdIndex exists_;
  };
}  // namespace nimbus

#endif  // NIMBUS_SHARED_LDO_INDEX_H_
