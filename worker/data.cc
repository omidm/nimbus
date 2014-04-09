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

#include "worker/data.h"

using namespace nimbus; // NOLINT

Data::Data() {}

Data* Data::Clone() {
  std::cout << "cloning the base data\n";
  Data* d = new Data();
  return d;
}

logical_data_id_t Data::logical_id() {
  return logical_id_;
}

physical_data_id_t Data::physical_id() {
  return physical_id_;
}

std::string Data::name() {
    return name_;
}

GeometricRegion Data::region() {
    return region_;
}

data_version_t Data::version() {
  return version_;
}

void Data::set_logical_id(logical_data_id_t logical_id) {
  logical_id_ = logical_id;
}

void Data::set_physical_id(physical_data_id_t physical_id) {
  physical_id_ = physical_id;
}
int Data::get_debug_info() {
    return -1;
}

void Data::set_name(std::string name) {
    name_ = name;
}

void Data::set_region(const GeometricRegion& region) {
    region_ = region;
}

void Data::set_version(data_version_t version) {
  version_ = version;
}



