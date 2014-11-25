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
  * Distributed DB class implements a service that is used to save the state
  * required for checkpoint creation in order to provide fault tolerance.
  * 
  * It is basically a key-value store. Upon saving a value with a specific key,
  * the service returns a handle that can be used later to retrieve the value.
  * The handle is overloaded version of the key such that any node in the
  * cluster could find the primary node that saved the data and retrieve it
  * form that node. This way even if a worker is down, other nodes could access
  * the saved state from the non-volatile memory that worker save the data on.
  *
  * Currently each node of the distributed db has a LevelDB module to save the
  * value on local disk. The handle interface makes it easy to change the
  * implementation anytime in the future to use other options like AWS S3.
  * For LevelDB implementation handle has the following format:
  * 
  * handle: <node_ip_address>-<leveldb_root>-<key>
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_DISTRIBUTED_DB_H_
#define NIMBUS_SHARED_DISTRIBUTED_DB_H_

#include <boost/thread.hpp>
#include <stdlib.h>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <map>
#include <vector>
#include <string>
#include "shared/nimbus_types.h"
#include "leveldb/db.h"

namespace nimbus {


class DistributedDB {
  public:
    DistributedDB();
    ~DistributedDB();

    void Initialize(const std::string& ip_address,
                    const worker_id_t& worker_id);

    bool Put(const std::string& key,
             const std::string& value,
             const checkpoint_id_t& checkpoint_id,
             std::string *handle);

    bool Get(const std::string& handle,
             std::string *value);

    bool RemoveCheckpoint(checkpoint_id_t checkpoint_id);


  private:
    boost::mutex mutex_;
    bool initialized_;
    std::string ip_address_;
    worker_id_t worker_id_;

    typedef std::map<std::string, leveldb::DB*> Map;
    Map db_map_;

    leveldb::DB* GetDB(const std::string& ip_address,
                       const std::string& leveldb_root);

    bool DBExistsLocally(std::string leveldb_root);

    bool RetrieveDBFromOtherNode(const std::string& ip_address,
                                 const std::string& leveldb_root);
};



}  // namespace nimbus


#endif  // NIMBUS_SHARED_DISTRIBUTED_DB_H_

