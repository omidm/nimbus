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
 *   FILE: .//ldo_index.cc
 *   DATE: Fri Oct 25 10:20:13 2013
 *  DESCR:
 ***********************************************************************/
#include "shared/ldo_index.h"

namespace nimbus {
/**
 * \fn nimbus::LdoIndex::LdoIndex()
 * \brief Brief description.
 * \return
*/
LdoIndex::LdoIndex() {}


/**
 * \fn LdoIndex::~LdoIndex()
 * \brief Brief description.
 * \return
*/
LdoIndex::~LdoIndex() {}


/**
 * \fn void LdoIndex::addObject(LogicalDataObject *object)
 * \brief Brief description.
 * \param object
 * \return
*/
void LdoIndex::addObject(LogicalDataObject *object) {
  std::string var = object->variable();
  LdoList* list;
  if (index_.find(var) == index_.end()) {
    list = new LdoList();
    index_[var] = list;
  } else {
    list = index_[var];
  }

  list->push_back(object);
}


/**
 * \fn LdoList * LdoIndex::intersectingObjects(std::string variable,
                                                       GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
LdoVector * LdoIndex::intersectingObjects(std::string variable,
                                                  GeometricRegion *region) {
  LdoVector* output = new LdoVector();
  if (index_.find(variable) == index_.end()) {  // No such variable
    return output;
  }

  LdoList* list = index_[variable];
  LdoList::iterator iter = list->begin();
  for (; iter != list->end(); ++iter) {
    LogicalDataObject* object = *iter;
    if (region->intersects(object->region())) {
      output->push_back(object);
    }
  }
  return output;
}


/**
 * \fn LdoList * LdoIndex::coveredObjects(std::string variable,
                                                  GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
LdoVector * LdoIndex::coveredObjects(std::string variable,
                                             GeometricRegion *region) {
  LdoVector* output = new LdoVector();
  if (index_.find(variable) == index_.end()) {  // No such variable
    return output;
  }

  LdoList* list = index_[variable];
  LdoList::iterator iter = list->begin();
  for (; iter != list->end(); ++iter) {
    LogicalDataObject* object = *iter;
    if (region->covers(object->region())) {
      output->push_back(object);
    }
  }
  return output;
}


/**
 * \fn LdoList * LdoIndex::adjacentObjects(std::string variable,
                                                   GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
LdoVector * LdoIndex::adjacentObjects(std::string variable,
                                              GeometricRegion *region) {
  LdoVector* output = new LdoVector();
  if (index_.find(variable) == index_.end()) {  // No such variable
    return output;
  }

  LdoList* list = index_[variable];
  LdoList::iterator iter = list->begin();
  for (; iter != list->end(); ++iter) {
    LogicalDataObject* object = *iter;
    if (region->adjacent(object->region())) {
      output->push_back(object);
    }
  }
  return output;
}


}  // namespace nimbus
