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

/***********************************************************************
 * AUTHOR: Philip Levis <pal>
 *   FILE: .//data_map.cc
 *   DATE: Fri Nov  1 12:44:48 2013
 *  DESCR:
 ***********************************************************************/
#include "scheduler/data_map.h"

namespace nimbus {
/**
 * \fn nimbus::DataMap::DataMap()
 * \brief Brief description.
 * \return
*/
nimbus::DataMap::DataMap() {}


/**
 * \fn nimbus::DataMap::~DataMap()
 * \brief Brief description.
 * \return
*/
nimbus::DataMap::~DataMap() {
  DataMapType::iterator it = data_map_.begin();
  while (it != data_map_.end()) {
    std::pair<data_id_t, PhysicalDataVector*> pair = *it;
    delete pair.second;
  }
}


/**
 * \fn void nimbus::DataMap::AddLogicalObject(LogicalDataObject *object)
 * \brief Brief description.
 * \param object
 * \return
*/
void nimbus::DataMap::AddLogicalObject(LogicalDataObject *object) {
  data_id_t id = object->id();

  // Does not exist, insert
  if (data_map_.find(id) == data_map_.end()) {
    PhysicalDataVector* v = new PhysicalDataVector();
    data_map_[id] = v;
  } else {   // Exists, error
  }
}


/**
 * \fn void nimbus::DataMap::AddPhysicalInstance(PhysicalData *instance)
 * \brief Brief description.
 * \param instance
 * \return
*/
bool nimbus::DataMap::AddPhysicalInstance(LogicalDataObject* obj,
                                          PhysicalData instance) {
  if (data_map_.find(obj->id()) == data_map_.end()) {
    return false;
  } else {
    PhysicalDataVector* v = data_map_[obj->id()];
    v->push_back(instance);
    return true;
  }
}


/**
 * \fn void nimbus::DataMap::RemovePhysicalInstance(PhysicalData *instance)
 * \brief Brief description.
 * \param instance
 * \return
*/
bool nimbus::DataMap::RemovePhysicalInstance(LogicalDataObject* obj,
                                             PhysicalData instance) {
  if (data_map_.find(obj->id()) == data_map_.end()) {
    return false;
  } else {
    PhysicalDataVector* v = data_map_[obj->id()];
    PhysicalDataVector::iterator it = v->begin();
    for (; it != v->end(); ++it) {
      PhysicalData pd = *it;
      if (pd.id() == instance.id()) {
        v->erase(it);
        return true;
      }
    }
  }
  return false;
}


/**
 * \fn const PhysicalDataVector * nimbus::DataMap::AllInstances(LogicalDataObject *object)
 * \brief Brief description.
 * \param object
 * \return
*/
const PhysicalDataVector * nimbus::DataMap::AllInstances(LogicalDataObject *object) {
  if (data_map_.find(object->id()) == data_map_.end()) {
    return NULL;
  } else {
    PhysicalDataVector* v = data_map_[object->id()];
    return v;
  }
}


/**
 * \fn int nimbus::DataMap::AllInstances(LogicalDataObject *object,
                              PhysicalDataVector *dest)
 * \brief Brief description.
 * \param object
 * \param dest
 * \return
*/
int nimbus::DataMap::AllInstances(LogicalDataObject *object,
                                  PhysicalDataVector *dest) {
  if (data_map_.find(object->id()) == data_map_.end()) {
    return 0;
  } else {
    PhysicalDataVector* v = data_map_[object->id()];
    int len = v->size();
    for (int i = 0; i < len; ++i) {
      dest[i] = v[i];
    }
    return len;
  }
}


/**
 * \fn int nimbus::DataMap::InstancesByWorker(LogicalDataObject *object,
                                   worker_id_t worker,
                                   PhysicalDataVector *dest)
 * \brief Brief description.
 * \param object
 * \param worker
 * \param dest
 * \return
*/
int nimbus::DataMap::InstancesByWorker(LogicalDataObject *object,
                                   worker_id_t worker,
                                   PhysicalDataVector *dest) {
  if (data_map_.find(object->id()) == data_map_.end()) {
    return 0;
  } else {
    PhysicalDataVector* v = data_map_[object->id()];
    PhysicalDataVector::iterator it = v->begin();
    int count = 0;

    for (; it != v->end(); ++it) {
      PhysicalData pd = *it;
      if (pd.worker() == worker) {
        (*dest)[count] = pd;
        count++;
      }
    }

    return count;
  }
}


/**
 * \fn int nimbus::DataMap::InstancesByVersion(LogicalDataObject *object,
                                    data_version_t version,
                                    PhysicalDataVector *dest)
 * \brief Brief description.
 * \param object
 * \param version
 * \param dest
 * \return
*/
int nimbus::DataMap::InstancesByVersion(LogicalDataObject *object,
                                        data_version_t version,
                                        PhysicalDataVector *dest) {
  if (data_map_.find(object->id()) == data_map_.end()) {
    return 0;
  } else {
    PhysicalDataVector* v = data_map_[object->id()];
    PhysicalDataVector::iterator it = v->begin();
    int count = 0;
    for (; it != v->end(); ++it) {
      PhysicalData pd = *it;
      if (pd.version() == version) {
        (*dest)[count] = pd;
        count++;
      }
    }
    return count;
  }
}

}  // namespace nimbus
