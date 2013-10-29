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
#include "application_utils/geometric_region.h"

namespace nimbus {
/**
 * \fn nimbus::GeometricRegion::GeometricRegion(int_dimension_t x,
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
nimbus::GeometricRegion::GeometricRegion(int_dimension_t x,
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

/**
 * \fn bool nimbus::GeometricRegion::intersects(GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
bool nimbus::GeometricRegion::adjacentOrIntersects(GeometricRegion *region) {
  return !((x() + dx() < region->x()) ||
           (x() > region->x() + region->dx()) ||
           (y() + dy() < region->y()) ||
           (y() > region->y() + region->dy()) ||
           (z() + dz() < region->z()) ||
           (z() > region->z() + region->dz()));
}

/**
 * \fn bool nimbus::GeometricRegion::intersects(GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return
*/
bool nimbus::GeometricRegion::intersects(GeometricRegion *region) {
  return !((x() + dx() <= region->x()) ||
           (x() >= region->x() + region->dx()) ||
           (y() + dy() <= region->y()) ||
           (y() >= region->y() + region->dy()) ||
           (z() + dz() <= region->z()) ||
           (z() >= region->z() + region->dz()));
}

/**
 * \fn bool nimbus::GeometricRegion::adjacent(GeometricRegion *region)
 * \brief Brief description.
 * \param region
 * \return Whether the two regions share face surfaces.
*/
bool nimbus::GeometricRegion::adjacent(GeometricRegion *region) {
  return (adjacentOrIntersects(region) && !intersects(region));
}


/**
 * \fn bool nimbus::GeometricRegion::covers(GeometricRegion *region)
 * \brief Returns whether this region (encompasses) covers the region
          described by the argument.
 * \param region
 * \return
*/
bool nimbus::GeometricRegion::covers(GeometricRegion *region) {
  return ((x() <= region->x()) &&
          (x() + dx() >= region->x() + region->dx()) &&
          (y() <= region->y()) &&
          (y() + dy() >= region->y() + region->dy()) &&
          (z() <= region->z()) &&
          (z() + dz() >= region->z() + region->dz()));
}


/**
 * \fn int_dimension_t nimbus::GeometricRegion::x()
 * \brief Brief description.
 * \return
*/
int_dimension_t nimbus::GeometricRegion::x() {
  return x_;
}


/**
 * \fn int_dimension_t nimbus::GeometricRegion::y()
 * \brief Brief description.
 * \return
*/
int_dimension_t nimbus::GeometricRegion::y() {
  return y_;
}


/**
 * \fn int_dimension_t nimbus::GeometricRegion::z()
 * \brief Brief description.
 * \return
*/
int_dimension_t nimbus::GeometricRegion::z() {
  return z_;
}


/**
 * \fn int_dimension_t nimbus::GeometricRegion::dx()
 * \brief Brief description.
 * \return
*/
int_dimension_t nimbus::GeometricRegion::dx() {
  return dx_;
}


/**
 * \fn int_dimension_t nimbus::GeometricRegion::dy()
 * \brief Brief description.
 * \return
*/
int_dimension_t nimbus::GeometricRegion::dy() {
  return dy_;
}


/**
 * \fn int_dimension_t nimbus::GeometricRegion::dz()
 * \brief Brief description.
 * \return
*/
int_dimension_t nimbus::GeometricRegion::dz() {
  return dz_;
}

}  // namespace nimbus
