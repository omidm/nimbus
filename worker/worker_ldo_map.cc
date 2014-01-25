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
 *   FILE: .//worker_ldo_map.cc
 *   DATE: Mon Nov 11 13:27:54 2013
 *  DESCR:
 ***********************************************************************/
#include "worker/worker_ldo_map.h"

namespace nimbus {

/**
 * \fn nimbus::WorkerLdoMap::WorkerLdoMap()
 * \brief Brief description.
 * \return
*/
nimbus::WorkerLdoMap::WorkerLdoMap() {}


/**
 * \fn nimbus::WorkerLdoMap::~WorkerLdoMap()
 * \brief Brief description.
 * \return
*/
nimbus::WorkerLdoMap::~WorkerLdoMap() {
  CLdoVector all;
  ldo_index_.AllObjects(&all);
  CLdoVector::iterator it = all.begin();
  for (; it != all.end(); ++it) {
    const LogicalDataObject* obj = *it;
    delete obj;
  }
}


/**
 * \fn bool nimbus::WorkerLdoMap::AddPartition(partition_id_t id,
                                              GeometricRegion region)
 * \brief Brief description.
 * \param id
 * \param region
 * \return
*/
bool nimbus::WorkerLdoMap::AddPartition(partition_id_t id,
                                       GeometricRegion r) {
  dbg(DBG_DATA_OBJECTS, "Adding %llu partition.\n", id);
  if (HasPartition(id)) {
    GeometricRegion r_p;
    FindPartition(id, &r_p);
    if (r_p.IsEqual(&r)) {
      return true;
    } else {
      dbg(DBG_DATA_OBJECTS|DBG_ERROR, "  - FAIL WorkerLdoMap: tried adding existing partition %llu.\n", id); // NOLINT
      return false;
    }
  } else {
    partition_map_.insert(std::pair<partition_id_t, GeometricRegion>(id, r));
    return true;
  }
}

/**
 * \fn bool nimbus::WorkerLdoMap::RemovePartition(partition_id_t id)
 *
 * \brief Brief description.
 * \param id
 * \return
*/
bool nimbus::WorkerLdoMap::RemovePartition(partition_id_t id) {
  dbg(DBG_DATA_OBJECTS, "Removing partition %llu.\n", id);
  if (!HasPartition(id)) {
    dbg(DBG_DATA_OBJECTS|DBG_ERROR, "  - FAIL WorkerLdoMap: tried removing non-existent partition %llu.\n", id); // NOLINT
    return false;
  } else {
    partition_map_.erase(id);
    return true;
  }
}

/**
 * \fn bool nimbus::WorkerLdoMap::HasPartition(partition_id_t id)
 *
 * \brief Brief description.
 * \param id
 * \return
*/
bool nimbus::WorkerLdoMap::HasPartition(partition_id_t id) {
  return partition_map_.find(id) != partition_map_.end();
}

/**
 * \fn bool nimbus::WorkerLdoMap::FindPartition(partition_id_t id, GeometricRegion* r)
 *
 * \brief Brief description.
 * \param id
 * \return
*/
bool nimbus::WorkerLdoMap::FindPartition(partition_id_t id, GeometricRegion* r) {
  if (HasPartition(id)) {
    *r =  partition_map_[id];
    return true;
  } else {
    dbg(DBG_ERROR, "  - FAIL WorkerLdoMap: partition id %llu does not exist.\n", id); // NOLINT
    return false;
  }
}

/**
 * \fn bool nimbus::WorkerLdoMap::AddLogicalObject(logical_data_id_t id,
                                       std::string variable,
                                       GeometricRegion region)
 * \brief Brief description.
 * \param id
 * \param variable
 * \param region
 * \return
*/
bool nimbus::WorkerLdoMap::AddLogicalObject(logical_data_id_t id,
                                            std::string variable,
                                            GeometricRegion region) {
  if (ldo_index_.HasObject(id)) {
    LogicalDataObject* ldo = ldo_index_.SpecificObject(id);
    if (ldo->variable() == variable && region.IsEqual(ldo->region())) {
      return true;
    } else {
      dbg(DBG_DATA_OBJECTS|DBG_ERROR, "  - FAIL WorkerLdoMap: tried adding existing object %llu.\n", id); // NOLINT
      return false;
    }
  } else {
    GeometricRegion* r = new GeometricRegion(region);
    LogicalDataObject* obj = new LogicalDataObject(id, variable, r);
    return ldo_index_.AddObject(obj);
  }
}

/**
 * \fn bool nimbus::WorkerLdoMap::AddLogicalObject(logical_data_id_t id,
                                      std::string variable,
                                      partition_id_t partition)
 * \brief Brief description.
 * \param id
 * \param variable
 * \param partition
 * \return
*/
bool nimbus::WorkerLdoMap::AddLogicalObject(logical_data_id_t id,
                                           std::string variable,
                                           partition_id_t partition) {
  dbg(DBG_DATA_OBJECTS, "Adding %llu as type %s.\n", id, variable.c_str());
  if (!HasPartition(partition)) {
    dbg(DBG_DATA_OBJECTS|DBG_ERROR, "  - FAIL WorkerLdoMap: tried adding object %llu with invalid partition %llu.\n", id, partition); // NOLINT
    return false;
  } else {
    GeometricRegion r;
    FindPartition(partition, &r);
    return AddLogicalObject(id, variable, r);
  }
}

/**
 * \fn bool nimbus::WorkerLdoMap::RemoveLogicalObject(logical_data_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
bool nimbus::WorkerLdoMap::RemoveLogicalObject(logical_data_id_t id) {
  if (!ldo_index_.HasObject(id)) {
    return false;
  } else {
    LogicalDataObject* obj = ldo_index_.SpecificObject(id);
    ldo_index_.RemoveObject(id);
    delete obj;
    return true;
  }
}


/**
 * \fn const LogicalDataObject * nimbus::WorkerLdoMap::FindLogicalObject(logical_data_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
const LogicalDataObject * nimbus::WorkerLdoMap::FindLogicalObject(logical_data_id_t id) {
  return ldo_index_.SpecificObject(id);
}


/**
 * \fn int nimbus::WorkerLdoMap::FindLogicalObjects(std::string variable,
                                         CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param dest
 * \return
*/
int nimbus::WorkerLdoMap::FindLogicalObjects(const std::string& variable,
                                             CLdoVector *dest) {
  return ldo_index_.AllObjects(variable, dest);
}


/**
 * \fn int nimbus::WorkerLdoMap::FindIntersectingLogicalObjects(std::string variable,
                                                     GeometricRegion *region,
                                                     CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param region
 * \param dest
 * \return
*/
int nimbus::WorkerLdoMap::FindIntersectingLogicalObjects(
    const std::string& variable,
    const GeometricRegion* region,
    CLdoVector *dest) {
  return ldo_index_.IntersectingObjects(variable, region, dest);
}


/**
 * \fn int nimbus::WorkerLdoMap::FindCoveredLogicalObjects(std::string variable,
                                                GeometricRegion *region,
                                                CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param region
 * \param dest
 * \return
*/
int nimbus::WorkerLdoMap::FindCoveredLogicalObjects(
    const std::string& variable,
    const GeometricRegion* region,
    CLdoVector *dest) {
  return ldo_index_.CoveredObjects(variable, region, dest);
}


/**
 * \fn int nimbus::WorkerLdoMap::FindAdjacentLogicalObjects(std::string variable,
                                                 GeometricRegion *region,
                                                 CLdoVector *dest)
 * \brief Brief description.
 * \param variable
 * \param region
 * \param dest
 * \return
*/
int nimbus::WorkerLdoMap::FindAdjacentLogicalObjects(
    const std::string& variable,
    const GeometricRegion* region,
    CLdoVector *dest) {
  return ldo_index_.AdjacentObjects(variable, region, dest);
}

}  // namespace nimbus
