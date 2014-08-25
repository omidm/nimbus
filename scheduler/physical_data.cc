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
 *   FILE: .//physical_data.cc
 *   DATE: Fri Nov  1 10:41:39 2013
 *  DESCR:
 ***********************************************************************/
#include "scheduler/physical_data.h"

namespace nimbus {

PhysicalData::PhysicalData() {}

PhysicalData::PhysicalData(const physical_data_id_t& id, const worker_id_t& worker) {
  id_ = id;
  worker_ = worker;
  version_ = NIMBUS_INIT_DATA_VERSION;
  last_job_write_ = NIMBUS_KERNEL_JOB_ID;
}


PhysicalData::PhysicalData(const physical_data_id_t& id, const worker_id_t& worker,
    const data_version_t& version) {
  id_ = id;
  worker_ = worker;
  version_ = version;
  last_job_write_ = NIMBUS_KERNEL_JOB_ID;
}


PhysicalData::PhysicalData(const physical_data_id_t& id, const worker_id_t& worker,
    const data_version_t& version,
    const IDSet<job_id_t>& list_job_read, const job_id_t& last_job_write) {
  id_ = id;
  worker_ = worker;
  version_ = version;
  list_job_read_ = list_job_read;
  last_job_write_ = last_job_write;
}

PhysicalData::PhysicalData(const PhysicalData& other) {
  id_ = other.id_;
  worker_ = other.worker_;
  version_ = other.version_;
  list_job_read_ = other.list_job_read_;
  last_job_write_ = other.last_job_write_;
}

PhysicalData& nimbus::PhysicalData::operator= (const PhysicalData& right) {
  id_ = right.id_;
  worker_ = right.worker_;
  version_ = right.version_;
  list_job_read_ = right.list_job_read_;
  last_job_write_ = right.last_job_write_;
  return *this;
}

/**
 * \fn nimbus::PhysicalData::~PhysicalData()
 * \brief Brief description.
 * \return
*/
nimbus::PhysicalData::~PhysicalData() {}


/**
 * \fn physical_data_id_t nimbus::PhysicalData::id()
 * \brief Brief description.
 * \return
*/
physical_data_id_t nimbus::PhysicalData::id() const {
  return id_;
}


/**
 * \fn worker_id_t nimbus::PhysicalData::worker()
 * \brief Brief description.
 * \return
*/
worker_id_t nimbus::PhysicalData::worker() const {
  return worker_;
}


/**
 * \fn data_version_t nimbus::PhysicalData::version()
 * \brief Brief description.
 * \return
*/
data_version_t nimbus::PhysicalData::version() const {
  return version_;
}


IDSet<job_id_t> nimbus::PhysicalData::list_job_read() const {
  return list_job_read_;
}


const IDSet<job_id_t>* nimbus::PhysicalData::list_job_read_p() const {
  return &list_job_read_;
}

IDSet<job_id_t>* nimbus::PhysicalData::list_job_read_p() {
  return &list_job_read_;
}


/**
 * \fn job_id_t nimbus::PhysicalData::last_job_write()
 * \brief Brief description.
 * \return
*/
job_id_t nimbus::PhysicalData::last_job_write() const {
  return last_job_write_;
}

void nimbus::PhysicalData::set_id(physical_data_id_t id) {
  id_ = id;
}

void nimbus::PhysicalData::set_worker(worker_id_t worker) {
  worker_ = worker;
}

/**
 * \fn void nimbus::PhysicalData::set_version(data_version_t v)
 * \brief Brief description.
 * \param v
 * \return
*/
void nimbus::PhysicalData::set_version(data_version_t v) {
  version_ = v;
}

void nimbus::PhysicalData::set_list_job_read(IDSet<job_id_t> list_job_read) {
  list_job_read_ = list_job_read;
}

/**
 * \fn void nimbus::PhysicalData::set_last_job_write(job_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
void nimbus::PhysicalData::set_last_job_write(job_id_t id) {
  last_job_write_ = id;
}

void nimbus::PhysicalData::clear_list_job_read() {
  list_job_read_.clear();
}

void nimbus::PhysicalData::add_to_list_job_read(job_id_t job_id) {
  list_job_read_.insert(job_id);
}

}  // namespace nimbus
