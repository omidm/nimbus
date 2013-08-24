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
  * Nimbus abstraction of data. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_WORKER_DATA_H_
#define NIMBUS_WORKER_DATA_H_

#include <set>
#include <map>
#include <string>
#include "shared/cluster.h"
#include "shared/idset.h"
#include "shared/nimbus_types.h"

namespace nimbus {
class Data;
typedef std::set<Data*> Neighbors;
typedef std::map<data_id_t, Data*> DataMap;
typedef std::map<std::string, Data*> DataTable;

class Data {
 public:
  Data();
  virtual ~Data() {}

  virtual Data* clone();
  virtual void create() {}
  virtual void destroy(data_id_t id) {}
  virtual void LocalCopy(data_id_t id_source, data_id_t id_destination) {}
  // TODO(omidm): Computer needs to be changed to worker, will be changed soon.
  virtual void RemoteCopy(Computer source, data_id_t id_source,
      Computer destination, data_id_t id_destination) {}


  virtual void duplicate(Computer source, Computer destination) {}
  virtual void migrate(Computer sourcer, Computer destination) {}
  virtual void split(Data *, Data *) {}
  virtual void merge(Data *, Data *) {}

  virtual int get_debug_info();

  data_id_t id();
  void set_id(data_id_t id);

 private:
  data_id_t id_;
  partition_t partition_id_;
  bool advanceData_;
  Hosts hosts_;

  // Set of data ids that could be involved in SYNC jobs with this data.
  IDSet<data_id_t> neighbors_;

  // Set of partition ids neighbor to this partition.
  IDSet<data_id_t> neighbor_partitions_;
};

}  // namespace nimbus

#endif  // NIMBUS_WORKER_DATA_H_
