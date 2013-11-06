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
 *   FILE: .//data_manager.cc
 *   DATE: Sun Nov  3 10:34:57 2013
 *  DESCR:
 ***********************************************************************/
#include "scheduler/data_manager.h"
#include "shared/dbg.h"

namespace nimbus {
/**
 * \fn nimbus::DataManager::DataManager()
 * \brief Brief description.
 * \return
*/
nimbus::DataManager::DataManager() {}


/**
 * \fn nimbus::DataManager::~DataManager()
 * \brief Because the DataManager allocates the LogicalDataObjects,
 *        this destructor must free them.
 * \return
*/
nimbus::DataManager::~DataManager() {
  std::map<data_id_t, LogicalDataObject*>::iterator it = ldo_map_.begin();
  for (; it != ldo_map_.end(); ++it) {
    LogicalDataObject* o = (*it).second;
    dbg(DBG_DATA_OBJECTS|DBG_MEMORY, "Invoking delete on object %llu.\n", o->id());
    delete o;
  }
  ldo_map_.clear();
}


/**
 * \fn bool nimbus::DataManager::AddLogicalObject(data_id_t id,
                                      std::string variable,
                                      GeometricRegion *region)
 * \brief Brief description.
 * \param id
 * \param variable
 * \param region
 * \return
*/
bool nimbus::DataManager::AddLogicalObject(data_id_t id,
                                           std::string variable,
                                           GeometricRegion r) {
  dbg(DBG_DATA_OBJECTS, "Adding %llu as type %s.\n", id, variable.c_str());
  if (ldo_index_.HasObject(id)) {
    dbg(DBG_DATA_OBJECTS|DBG_ERROR, "  - FAIL DataManager: tried adding existing object %llu.\n", id); // NOLINT
    return false;
  } else {
    GeometricRegion* region = new GeometricRegion(r);
    dbg(DBG_MEMORY, "Allocated geo region 0x%x\n", region);
    LogicalDataObject* ldo = new LogicalDataObject(id, variable, region);
    dbg(DBG_MEMORY, "Allocated ldo 0x%x\n", ldo);
    ldo_map_[id] = ldo;
    data_map_.AddLogicalObject(ldo);
    ldo_index_.AddObject(ldo);
    return true;
  }
}


/**
 * \fn bool nimbus::DataManager::RemoveLogicalObject(data_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
bool nimbus::DataManager::RemoveLogicalObject(data_id_t id) {
  if (!ldo_index_.HasObject(id)) {
    dbg(DBG_ERROR, "  - FAIL DataManager: tried deleting non-existent object %llu.\n", id); // NOLINT
    return false;
  } else {
    dbg(DBG_DATA_OBJECTS, "Removing logical object %llu from data manager.\n", id);
    LogicalDataObject* obj = ldo_map_[id];
    ldo_map_.erase(id);
    data_map_.RemoveLogicalObject(id);
    ldo_index_.RemoveObject(id);
    delete obj;
    return true;
  }
}


/**
 * \fn const LogicalDataObject * nimbus::DataManager::FindLogicalObject(data_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
const LogicalDataObject * nimbus::DataManager::FindLogicalObject(data_id_t id) {
  if (ldo_map_.find(id) == ldo_map_.end()) {
    return NULL;
  } else {
    return ldo_map_[id];
  }
}


/**
 * \fn int nimbus::DataManager::FindLogicalObjects(std::string variable,
                                        CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param dest
 * \return
*/
int nimbus::DataManager::FindLogicalObjects(std::string variable,
                                            CLdoVector* dest) {
  return ldo_index_.AllObjects(variable, dest);
}


/**
 * \fn int nimbus::DataManager::FindIntersectingLogicalObjects(std::string variable,
                                                    GeometricRegion *region,
                                                    CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param region
 * \param dest
 * \return
*/
int nimbus::DataManager::FindIntersectingLogicalObjects(std::string variable,
                                                    GeometricRegion *region,
                                                        CLdoVector *dest) {
  return ldo_index_.IntersectingObjects(variable, region, dest);
}


/**
 * \fn int nimbus::DataManager::FindCoveredLogicalObjects(std::string variable,
                                               GeometricRegion *region,
                                               CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param region
 * \param dest
 * \return
*/
int nimbus::DataManager::FindCoveredLogicalObjects(std::string variable,
                                                   GeometricRegion *region,
                                                   CLdoVector *dest) {
  return ldo_index_.CoveredObjects(variable, region, dest);
}


/**
 * \fn int nimbus::DataManager::FindAdjacentLogicalObjects(std::string variable,
                                                GeometricRegion *region,
                                                CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param region
 * \param dest
 * \return
*/
int nimbus::DataManager::FindAdjacentLogicalObjects(std::string variable,
                                                    GeometricRegion *region,
                                                    CLdoVector *dest) {
  return ldo_index_.AdjacentObjects(variable, region, dest);
}


/**
 * \fn bool nimbus::DataManager::AddPhysicalInstance(LogicalDataObject *object,
                                         PhysicalData instance)
 * \brief Brief description.
 * \param object
 * \param instance
 * \return
*/
bool nimbus::DataManager::AddPhysicalInstance(LogicalDataObject *object,
                                              PhysicalData instance) {
  return data_map_.AddPhysicalInstance(object, instance);
}


/**
 * \fn bool nimbus::DataManager::RemovePhysicalInstance(LogicalDataObject *object,
                                            PhysicalData instance)
 * \brief Brief description.
 * \param object
 * \param instance
 * \return
*/
bool nimbus::DataManager::RemovePhysicalInstance(LogicalDataObject *object,
                                                 PhysicalData instance) {
  return data_map_.RemovePhysicalInstance(object, instance);
}


/**
 * \fn const PhysicalDataVector * nimbus::DataManager::AllInstances(LogicalDataObject *object)
 * \brief Brief description.
 * \param object
 * \return
*/
const PhysicalDataVector * nimbus::DataManager::AllInstances(LogicalDataObject *object) {
  return data_map_.AllInstances(object);
}


/**
 * \fn int nimbus::DataManager::AllInstances(LogicalDataObject *object,
                                  PhysicalDataVector *dest)
 * \brief Brief description.
 * \param object
 * \param dest
 * \return
*/
int nimbus::DataManager::AllInstances(LogicalDataObject *object,
                                  PhysicalDataVector *dest) {
  return data_map_.AllInstances(object, dest);
}


/**
 * \fn int nimbus::DataManager::AllInstancesByWorker(LogicalDataObject *object,
                                          worker_id_t worker,
                                          PhysicalDataVector *dest)
 * \brief Brief description.
 * \param object
 * \param worker
 * \param dest
 * \return
*/
int nimbus::DataManager::InstancesByWorker(LogicalDataObject *object,
                                          worker_id_t worker,
                                          PhysicalDataVector *dest) {
  return data_map_.InstancesByWorker(object, worker, dest);
}
/**
 * \fn int nimbus::DataManager::AllInstancesByVersion(LogicalDataObject *object,
                                           data_version_t version,
                                           PhysicalDataVector *dest)
 * \brief Brief description.
 * \param object
 * \param version
 * \param dest
 * \return
*/
int nimbus::DataManager::InstancesByVersion(LogicalDataObject *object,
                                           data_version_t version,
                                           PhysicalDataVector *dest) {
  return data_map_.InstancesByVersion(object, version, dest);
}


}  // namespace nimbus
