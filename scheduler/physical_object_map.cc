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
 * AUTHOR: Philip Levis <pal@cs.stanford.edu>
 * AUTHOR: Omid Mashayekhi <omidm@stanford.edu>
 *   FILE: .//physical_object_map.cc
 *   DATE: Fri Nov  1 12:44:48 2013
 *  DESCR:
 ***********************************************************************/
#include "scheduler/physical_object_map.h"
#include "shared/dbg.h"

namespace nimbus {
/**
 * \fn nimbus::PhysicalObjectMap::PhysicalObjectMap()
 * \brief Brief description.
 * \return
*/
nimbus::PhysicalObjectMap::PhysicalObjectMap() {
  is_static_ = false;
}


/**
 * \fn nimbus::PhysicalObjectMap::~PhysicalObjectMap()
 * \brief Brief description.
 * \return
*/
nimbus::PhysicalObjectMap::~PhysicalObjectMap() {
  PhysicalObjectMapType::iterator it = data_map_.begin();
  for (; it != data_map_.end(); ++it) {
    std::pair<physical_data_id_t, PhysicalDataList*> pair = *it;
    dbg(DBG_DATA_OBJECTS, "Freeing physical vector 0x%llx\n", pair.second);
    delete pair.second;
  }
}


/**
 * \fn void nimbus::PhysicalObjectMap::BuildDataMapViewIfNeeded()
 * \brief Brief description.
 * \param object
 * \return
*/
bool nimbus::PhysicalObjectMap::BuildDataMapViewIfNeeded() {
  if (is_static_) {
    return false;
  }
  data_map_view_.snapshot(data_map_);
  is_static_ = true;
  return true;
}

/**
 * \fn void nimbus::PhysicalObjectMap::AddLogicalObject(LogicalDataObject *object)
 * \brief Brief description.
 * \param object
 * \return
*/
bool nimbus::PhysicalObjectMap::AddLogicalObject(LogicalDataObject *object) {
  assert(!is_static_);
  logical_data_id_t id = object->id();

  // Does not exist, insert
  if (data_map_.find(id) == data_map_.end()) {
    PhysicalDataList* v = new PhysicalDataList();
    dbg(DBG_MEMORY, "Allocating physical vector 0x%llx\n", v);
    data_map_[id] = v;
    return true;
  } else {   // Exists, error
    return false;
  }
}

bool nimbus::PhysicalObjectMap::RemoveLogicalObject(logical_data_id_t id) {
  assert(!is_static_);
  if (data_map_.find(id) == data_map_.end()) {  // Exists
    return false;
  } else {
    PhysicalDataList* v = data_map_[id];
    dbg(DBG_MEMORY, "Freeing physical vector 0x%llx\n", v);
    delete v;
    data_map_.erase(id);
    return true;
  }
}
/**
 * \fn void nimbus::PhysicalObjectMap::AddPhysicalInstance(PhysicalData *instance)
 * \brief Brief description.
 * \param instance
 * \return
*/
bool nimbus::PhysicalObjectMap::AddPhysicalInstance(LogicalDataObject* obj,
                                          const PhysicalData& instance) {
  BuildDataMapViewIfNeeded();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(obj->id());
  if (iter == data_map_view_.end()) {
    return false;
  } else {
    PhysicalDataList* v = *iter;
    v->push_back(instance);
    return true;
  }
}


/**
 * \fn void nimbus::PhysicalObjectMap::RemovePhysicalInstance(PhysicalData *instance)
 * \brief Brief description.
 * \param instance
 * \return
*/
bool nimbus::PhysicalObjectMap::RemovePhysicalInstance(LogicalDataObject* obj,
                                             const PhysicalData& instance) {
  BuildDataMapViewIfNeeded();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(obj->id());
  if (iter == data_map_view_.end()) {
    return false;
  } else {
    PhysicalDataList* v = *iter;
    PhysicalDataList::iterator it = v->begin();
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

size_t nimbus::PhysicalObjectMap::RemoveAllInstanceByWorker(worker_id_t worker_id) {
  // Not optimized yet.
  dbg(DBG_WARN, "Lookup not optimized!\n");
  size_t count = 0;
  PhysicalObjectMapType::iterator iter = data_map_.begin();
  for (; iter != data_map_.end(); ++iter) {
    PhysicalDataList *pdv = iter->second;
    PhysicalDataList::iterator it = pdv->begin();
    for (; it != pdv->end();) {
      PhysicalData pd = *it;
      if (pd.worker() == worker_id) {
        pdv->erase(it++);
        ++count;
      } else {
        ++it;
      }
    }
  }

  return count;
}

size_t nimbus::PhysicalObjectMap::ResetAllInstances() {
  // Not optimized yet.
  dbg(DBG_WARN, "Lookup not optimized!\n");
  size_t count = 0;
  PhysicalObjectMapType::iterator iter = data_map_.begin();
  for (; iter != data_map_.end(); ++iter) {
    PhysicalDataList *pdv = iter->second;
    PhysicalDataList::iterator it = pdv->begin();
    for (; it != pdv->end(); ++it) {
      it->set_version(NIMBUS_INIT_DATA_VERSION);
      IDSet<job_id_t> empty;
      it->set_list_job_read(empty);
      it->set_last_job_write(NIMBUS_KERNEL_JOB_ID);
      ++count;
    }
  }

  return count;
}


/**
 * \fn void nimbus::PhysicalObjectMap::UpdatePhysicalInstance(PhysicalData *instance)
 * \brief Brief description.
 * \param instance
 * \return
*/
bool nimbus::PhysicalObjectMap::UpdatePhysicalInstance(LogicalDataObject* obj,
                                                       const PhysicalData& old_instance,
                                                       const PhysicalData& new_instance) {
  BuildDataMapViewIfNeeded();
  assert(old_instance.id() == new_instance.id());

  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(obj->id());
  if (iter == data_map_view_.end()) {
    return false;
  } else {
    PhysicalDataList* v = *iter;
    PhysicalDataList::iterator it = v->begin();
    for (; it != v->end(); ++it) {
      if (it->id() == old_instance.id()) {
        it->set_worker(new_instance.worker());
        it->set_version(new_instance.version());
        it->set_list_job_read(new_instance.list_job_read());
        it->set_last_job_write(new_instance.last_job_write());
        return true;
      }
    }
  }
  return false;
}

/**
 * \fn const PhysicalDataList * nimbus::PhysicalObjectMap::AllInstances(LogicalDataObject *object)
 * \brief Brief description.
 * \param object
 * \return
*/
const PhysicalDataList * nimbus::PhysicalObjectMap::AllInstances(LogicalDataObject *object) {
  BuildDataMapViewIfNeeded();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(object->id());
  if (iter == data_map_view_.end()) {
    return NULL;
  } else {
    PhysicalDataList* v = *iter;
    return v;
  }
}


/**
 * \fn int nimbus::PhysicalObjectMap::AllInstances(LogicalDataObject *object,
                              PhysicalDataList *dest)
 * \brief Brief description.
 * \param object
 * \param dest
 * \return
*/
int nimbus::PhysicalObjectMap::AllInstances(LogicalDataObject *object,
                                  PhysicalDataList *dest) {
  BuildDataMapViewIfNeeded();
  dest->clear();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(object->id());
  if (iter == data_map_view_.end()) {
    return 0;
  } else {
    PhysicalDataList* v = *iter;
    PhysicalDataList::iterator it = v->begin();
    for (; it != v->end(); ++it) {
      dest->push_back(*it);
    }
    return v->size();
  }
}


/**
 * \fn int nimbus::PhysicalObjectMap::InstancesByWorker(LogicalDataObject *object,
                                   worker_id_t worker,
                                   PhysicalDataList *dest)
 * \brief Brief description.
 * \param object
 * \param worker
 * \param dest
 * \return
*/
int nimbus::PhysicalObjectMap::InstancesByWorker(LogicalDataObject *object,
                                   worker_id_t worker,
                                   PhysicalDataList *dest) {
  BuildDataMapViewIfNeeded();
  dest->clear();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(object->id());
  if (iter == data_map_view_.end()) {
    return 0;
  } else {
    PhysicalDataList* v = *iter;
    PhysicalDataList::iterator it = v->begin();
    int count = 0;

    for (; it != v->end(); ++it) {
      if (it->worker() == worker) {
        dest->push_back(*it);
        count++;
      }
    }

    return count;
  }
}


/**
 * \fn int nimbus::PhysicalObjectMap::InstancesByVersion(LogicalDataObject *object,
                                    data_version_t version,
                                    PhysicalDataList *dest)
 * \brief Brief description.
 * \param object
 * \param version
 * \param dest
 * \return
*/
int nimbus::PhysicalObjectMap::InstancesByVersion(LogicalDataObject *object,
                                        data_version_t version,
                                        PhysicalDataList *dest) {
  BuildDataMapViewIfNeeded();
  dest->clear();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(object->id());
  if (iter == data_map_view_.end()) {
    return 0;
  } else {
    PhysicalDataList* v = *iter;
    PhysicalDataList::iterator it = v->begin();
    int count = 0;
    for (; it != v->end(); ++it) {
      if (it->version() == version) {
        dest->push_back(*it);
        count++;
      }
    }
    return count;
  }
}

/**
 * \fn int nimbus::PhysicalObjectMap::InstancesByWorkerAndVersion(LogicalDataObject *object,
                                   worker_id_t worker,
                                   data_version_t version,
                                   PhysicalDataList *dest)
 * \brief Brief description.
 * \param object
 * \param worker
 * \param version
 * \param dest
 * \return
*/
int nimbus::PhysicalObjectMap::InstancesByWorkerAndVersion(LogicalDataObject *object,
                                   worker_id_t worker,
                                   data_version_t version,
                                   PhysicalDataList *dest) {
  BuildDataMapViewIfNeeded();
  dest->clear();
  PhysicalObjectMapViewType::iterator iter = data_map_view_.find(object->id());
  if (iter == data_map_view_.end()) {
    return 0;
  } else {
    PhysicalDataList* v = *iter;
    PhysicalDataList::iterator it = v->begin();
    int count = 0;

    for (; it != v->end(); ++it) {
      if ((it->worker() == worker) && (it->version() == version)) {
        dest->push_back(*it);
        count++;
      }
    }

    return count;
  }
}







}  // namespace nimbus
