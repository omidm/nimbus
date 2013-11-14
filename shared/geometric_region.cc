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
 *   FILE: .//geometric_region.cc
 *   DATE: Thu Oct 24 16:13:23 2013
 *  DESCR:
 ***********************************************************************/
#include "shared/geometric_region.h"

namespace nimbus {

  GeometricRegion::GeometricRegion() {
    x_ = y_ = z_ = dx_ = dy_ = dz_ = -1;
  }

/**
 * \fn GeometricRegion::GeometricRegion(int_dimension_t x,
                                         int_dimension_t y,
                                         int_dimension_t z,
                                         int_dimension_t dx,
                                         int_dimension_t dy,
                                         int_dimension_t dz)
 * \brief Brief description.
 * \param x
 * \param y
 * \param z
 * \param dx
 * \param dy
 * \param dz
 * \return
*/
GeometricRegion::GeometricRegion(int_dimension_t x,
                                 int_dimension_t y,
                                 int_dimension_t z,
                                 int_dimension_t dx,
                                 int_dimension_t dy,
                                 int_dimension_t dz) {
  x_ = x;
  y_ = y;
  z_ = z;
  dx_ = dx;
  dy_ = dy;
  dz_ = dz;
}

GeometricRegion::GeometricRegion(const GeometricRegion& r) {
  x_ = r.x();
  y_ = r.y();
  z_ = r.z();
  dx_ = r.dx();
  dy_ = r.dy();
  dz_ = r.dz();
}

GeometricRegion::GeometricRegion(const int_dimension_t* values) {
  fillInValues(values);
}

GeometricRegion::GeometricRegion(const GeometricRegionMessage* msg) {
  fillInValues(msg);
}

GeometricRegion::GeometricRegion(std::istream* is) {
  GeometricRegionMessage msg;
  msg.ParseFromIstream(is);
  fillInValues(&msg);
}

GeometricRegion::GeometricRegion(const std::string& data) {
  GeometricRegionMessage msg;
  msg.ParseFromString(data);
  fillInValues(&msg);
}

void GeometricRegion::Rebuild(int_dimension_t x,
                                 int_dimension_t y,
                                 int_dimension_t z,
                                 int_dimension_t dx,
                                 int_dimension_t dy,
                                 int_dimension_t dz) {
  x_ = x;
  y_ = y;
  z_ = z;
  dx_ = dx;
  dy_ = dy;
  dz_ = dz;
}

void GeometricRegion::FillInMessage(GeometricRegionMessage* msg) {
  msg->set_x(x_);
  msg->set_y(y_);
  msg->set_z(z_);
  msg->set_dx(dx_);
  msg->set_dy(dy_);
  msg->set_dz(dz_);
}

/**
 * \fn bool GeometricRegion::intersects(GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
bool GeometricRegion::AdjacentOrIntersects(GeometricRegion *region) {
  return !((x() + dx() < region->x()) ||
           (x() > region->x() + region->dx()) ||
           (y() + dy() < region->y()) ||
           (y() > region->y() + region->dy()) ||
           (z() + dz() < region->z()) ||
           (z() > region->z() + region->dz()));
}

/**
 * \fn bool GeometricRegion::intersects(GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
bool GeometricRegion::Intersects(GeometricRegion *region) {
  return !((x() + dx() <= region->x()) ||
           (x() >= region->x() + region->dx()) ||
           (y() + dy() <= region->y()) ||
           (y() >= region->y() + region->dy()) ||
           (z() + dz() <= region->z()) ||
           (z() >= region->z() + region->dz()));
}

/**
 * \fn bool GeometricRegion::Adjacent(GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return Whether the two regions share face surfaces.
*/
bool GeometricRegion::Adjacent(GeometricRegion *region) {
  return (AdjacentOrIntersects(region) && !Intersects(region));
}


/**
 * \fn bool GeometricRegion::Covers(GeometricRegion *region)
 * \brief Returns whether this region (encompasses) covers the region
          described by the argument.
 * \param region
 * \return
*/
bool GeometricRegion::Covers(GeometricRegion *region) {
  return ((x() <= region->x()) &&
          (x() + dx() >= region->x() + region->dx()) &&
          (y() <= region->y()) &&
          (y() + dy() >= region->y() + region->dy()) &&
          (z() <= region->z()) &&
          (z() + dz() >= region->z() + region->dz()));
}


/**
 * \fn int_dimension_t GeometricRegion::x()
 * \brief Brief description.
 * \return
*/
int_dimension_t GeometricRegion::x() const {
  return x_;
}


/**
 * \fn int_dimension_t GeometricRegion::y()
 * \brief Brief description.
 * \return
*/
int_dimension_t GeometricRegion::y() const {
  return y_;
}


/**
 * \fn int_dimension_t GeometricRegion::z()
 * \brief Brief description.
 * \return
*/
int_dimension_t GeometricRegion::z() const {
  return z_;
}


/**
 * \fn int_dimension_t GeometricRegion::dx()
 * \brief Brief description.
 * \return
*/
int_dimension_t GeometricRegion::dx() const {
  return dx_;
}


/**
 * \fn int_dimension_t GeometricRegion::dy()
 * \brief Brief description.
 * \return
*/
int_dimension_t GeometricRegion::dy() const {
  return dy_;
}


/**
 * \fn int_dimension_t GeometricRegion::dz()
 * \brief Brief description.
 * \return
*/
int_dimension_t GeometricRegion::dz() const {
  return dz_;
}

/**
 * \fn int_dimension_t GeometricRegion::fillInValues()
 * \brief Fills in six parameters from array.
 * \return
*/
void GeometricRegion::fillInValues(const int_dimension_t* values) {
  x_  = values[0];
  y_  = values[1];
  z_  = values[2];
  dx_ = values[3];
  dy_ = values[4];
  dz_ = values[5];
}

void GeometricRegion::fillInValues(const GeometricRegionMessage* msg) {
  x_ = msg->x();
  y_ = msg->y();
  z_ = msg->z();
  dx_ = msg->dx();
  dy_ = msg->dy();
  dz_ = msg->dz();
}

std::string GeometricRegion::toString() {
  std::string str;
  char buf[2048];
  snprintf(buf, sizeof(buf), "bbox:%ld,%ld,%ld,%ld,%ld,%ld",
           x_, y_, z_, dx_, dy_, dz_);
  str += buf;
  return str;
}

}  // namespace nimbus
