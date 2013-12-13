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
 *   FILE: .//physical_data_object.cc
 *   DATE: Tue Dec 10 10:12:56 2013
 *  DESCR:
 ***********************************************************************/
#include "worker/physical_data_instance.h"

namespace nimbus {

/**
 * \fn nimbus::PhysicalDataInstance::PhysicalDataInstance()
 * \brief Brief description.
 * \return
*/
PhysicalDataInstance::PhysicalDataInstance() {
}


/**
 * \fn PhysicalDataInstance::PhysicalDataInstance(physical_data_id_t id,
					       LogicalDataObject *lobj,
					       Data *data,
					       data_version_t version)
 * \brief Brief description.
 * \param id
 * \param lobj
 * \param data
 * \param version
 * \return
*/
PhysicalDataInstance::PhysicalDataInstance(physical_data_id_t i,
                                           LogicalDataObject *lobj,
                                           Data *d,
                                           data_version_t v) {
  id_ = i;
  logical_object_ = lobj;
  data_ = d;
  version_ = v;
}


/**
 * \fn PhysicalDataInstance::~PhysicalDataInstance()
 * \brief Brief description.
 * \return
*/
PhysicalDataInstance::~PhysicalDataInstance() {
}


/**
 * \fn physical_data_id_t PhysicalDataInstance::id()
 * \brief Brief description.
 * \return
*/
physical_data_id_t PhysicalDataInstance::id() const {
  return id_;
}


/**
 * \fn std::string PhysicalDataInstance::variable()
 * \brief Brief description.
 * \return
*/
std::string PhysicalDataInstance::variable() const {
  return logical_object()->variable();
}


/**
 * \fn GeometricRegion * PhysicalDataInstance::region()
 * \brief Brief description.
 * \return
*/
GeometricRegion * PhysicalDataInstance::region() const {
  return logical_object()->region();
}


/**
 * \fn partition_id_t PhysicalDataInstance::partition()
 * \brief Brief description.
 * \return
*/
partition_id_t PhysicalDataInstance::partition() const {
  return logical_object()->partition();
}


/**
 * \fn LogicalDataObject * PhysicalDataInstance::logical_object()
 * \brief Brief description.
 * \return
*/
LogicalDataObject * PhysicalDataInstance::logical_object() const {
  return logical_object_;
}


/**
 * \fn data_version_t PhysicalDataInstance::version()
 * \brief Brief description.
 * \return
*/
data_version_t PhysicalDataInstance::version() const {
  return version_;
}


/**
 * \fn void PhysicalDataInstance::set_version(data_version_t version)
 * \brief Brief description.
 * \param version
 * \return
*/
void PhysicalDataInstance::set_version(data_version_t v) {
  version_ = v;
}


/**
 * \fn data_version_t PhysicalDataInstance::IncrementVersion()
 * \brief Brief description.
 * \return
*/
data_version_t PhysicalDataInstance::IncrementVersion() {
  version_++;
  return version_;
}


/**
 * \fn Data * PhysicalDataInstance::data()
 * \brief Brief description.
 * \return
*/
Data * PhysicalDataInstance::data() const {
  return data_;
}

}  // namespace nimbus


