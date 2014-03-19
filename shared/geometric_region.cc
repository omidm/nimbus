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

using boost::tokenizer;
using boost::char_separator;

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

GeometricRegion::GeometricRegion(const Coord &min, const Coord &delta) {
  x_ = min.x;
  y_ = min.y;
  z_ = min.z;
  dx_ = delta.x;
  dy_ = delta.y;
  dz_ = delta.z;
}

GeometricRegion GeometricRegion::GeometricRegionFromRange(const Coord &min,
                                                          const Coord &max) {
  GeometricRegion result(min.x,
                         min.y,
                         min.z,
                         max.x - min.x + 1,
                         max.y - min.y + 1,
                         max.z - min.z + 1);
  return result;
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
bool GeometricRegion::AdjacentOrIntersects(GeometricRegion *region) const {
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
bool GeometricRegion::Intersects(GeometricRegion *region) const {
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
bool GeometricRegion::Adjacent(GeometricRegion *region) const {
  return (AdjacentOrIntersects(region) && !Intersects(region));
}


/**
 * \fn bool GeometricRegion::Covers(GeometricRegion *region)
 * \brief Returns whether this region (encompasses) covers the region
          described by the argument.
 * \param region
 * \return
*/
bool GeometricRegion::Covers(GeometricRegion *region) const {
  return ((x() <= region->x()) &&
          (x() + dx() >= region->x() + region->dx()) &&
          (y() <= region->y()) &&
          (y() + dy() >= region->y() + region->dy()) &&
          (z() <= region->z()) &&
          (z() + dz() >= region->z() + region->dz()));
}

/**
 * \fn bool GeometricRegion::IsEqual(GeometricRegion *region)
 * \brief Returns whether this region is equal to the region
          described by the argument.
 * \param region
 * \return
*/
bool GeometricRegion::IsEqual(GeometricRegion *region) const {
  return ((x() == region->x()) &&
          (y() == region->y()) &&
          (z() == region->z()) &&
          (dx() == region->dx()) &&
          (dy() == region->dy()) &&
          (dz() == region->dz()));
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

Coord GeometricRegion::MinCorner() const {
  Coord result(x_, y_, z_);
  return result;
}

Coord GeometricRegion::MaxCorner() const {
  Coord result(x_ + dx_ - 1, y_ + dy_ - 1, z_ + dz_ - 1);
  return result;
}

Coord GeometricRegion::Delta() const {
  Coord result(dx_, dy_, dz_);
  return result;
}

void GeometricRegion::Enlarge(const int_dimension_t delta) {
  x_ -= delta;
  y_ -= delta;
  z_ -= delta;
  dx_ += 2*delta;
  dy_ += 2*delta;
  dz_ += 2*delta;
}

void GeometricRegion::Enlarge(const Coord delta) {
  x_ -= delta.x;
  y_ -= delta.y;
  z_ -= delta.z;
  dx_ += 2*delta.x;
  dy_ += 2*delta.y;
  dz_ += 2*delta.z;
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

std::string GeometricRegion::toString() const {
  std::string str;
  char buf[2048];
  snprintf(buf, sizeof(buf), "bbox:%lld,%lld,%lld,%lld,%lld,%lld",
           x_, y_, z_, dx_, dy_, dz_);
  str += buf;
  return str;
}

bool GeometricRegion::Parse(const std::string& input) {
  std::string header = "bbox:";
  if (input.substr(0, header.size()) != header) {
    std::cout << "ERROR: wrong format for GeometricRegion" << std::endl;
    return false;
  }

  int num = 6;
  std::string str = input.substr(header.size());
  char_separator<char> separator(",");
  tokenizer<char_separator<char> > tokens(str, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();

  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: GeometricRegion has only " << i <<
        " fields (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: GeometricRegion has more than "<<
      num << " fields." << std::endl;
    return false;
  }

  int_dimension_t dim;

  iter = tokens.begin();
  std::stringstream ss_x(*iter);
  ss_x >> dim;
  if (ss_x.fail()) {
    std::cout << "ERROR: wrong element in bbox." << std::endl;
    return false;
  }
  x_ = dim;

  iter++;
  std::stringstream ss_y(*iter);
  ss_y >> dim;
  if (ss_y.fail()) {
    std::cout << "ERROR: wrong element in bbox." << std::endl;
    return false;
  }
  y_ = dim;

  iter++;
  std::stringstream ss_z(*iter);
  ss_z >> dim;
  if (ss_z.fail()) {
    std::cout << "ERROR: wrong element in bbox." << std::endl;
    return false;
  }
  z_ = dim;

  iter++;
  std::stringstream ss_dx(*iter);
  ss_dx >> dim;
  if (ss_dx.fail()) {
    std::cout << "ERROR: wrong element in bbox." << std::endl;
    return false;
  }
  dx_ = dim;

  iter++;
  std::stringstream ss_dy(*iter);
  ss_dy >> dim;
  if (ss_dy.fail()) {
    std::cout << "ERROR: wrong element in bbox." << std::endl;
    return false;
  }
  dy_ = dim;

  iter++;
  std::stringstream ss_dz(*iter);
  ss_dz >> dim;
  if (ss_dz.fail()) {
    std::cout << "ERROR: wrong element in bbox." << std::endl;
    return false;
  }
  dz_ = dim;

  return true;
}

GeometricRegion& GeometricRegion::operator= (const GeometricRegion& right) {
  x_ = right.x_;
  y_ = right.y_;
  z_ = right.z_;
  dx_ = right.dx_;
  dy_ = right.dy_;
  dz_ = right.dz_;
  return *this;
}

Coord::Coord() {
  x = 0;
  y = 0;
  z = 0;
}

Coord::Coord(int_dimension_t xe, int_dimension_t ye, int_dimension_t ze) {
  x = xe;
  y = ye;
  z = ze;
}

Coord ElementWiseMin(Coord a, Coord b) {
  Coord r;
  r.x = (a.x < b.x)? a.x : b.x;
  r.y = (a.y < b.y)? a.y : b.y;
  r.z = (a.z < b.z)? a.z : b.z;
  return r;
}

Coord ElementWiseMax(Coord a, Coord b) {
  Coord r;
  r.x = (a.x > b.x)? a.x : b.x;
  r.y = (a.y > b.y)? a.y : b.y;
  r.z = (a.z > b.z)? a.z : b.z;
  return r;
}
}  // namespace nimbus
