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
 * A CacheObject is an application object corresponding to one/ multiple nimbus
 * variables.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_CACHE_CACHE_OBJECT_H_
#define NIMBUS_DATA_CACHE_CACHE_OBJECT_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/cache_defs.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
typedef std::set<Data *> DataSet;
class GeometricRegion;

/**
 * \class CacheObject
 * \details Application object corresponding to one/ multiple nimbus variables.
 * CacheVariable and CacheStruct, that inherit from CacheObject, provide the
 * one variable and multiple variable implementation respectively.
 */
class CacheObject {
    public:
        /**
         * \brief Creates a CacheObject
         * \return Constructed CacheObject instance
         */
        explicit CacheObject();

        /**
         * \brief Creates a CacheObject
         * \param  ob_reg specifies application object region
         * \return Constructed CacheObject instance
         */
        explicit CacheObject(const GeometricRegion &ob_reg);

        /**
         * \brief Makes this instance a prototype. The application writer must
         * make a prototype for every application object he/ she plans to use.
         */
         void MakePrototype();

        /**
         * \brief Flushes data from cache, removes corresponding dirty data
         * mapping
         * \param d is data to flush to
         */
        virtual void PullData(Data *d) = 0;

        /**
         * \brief Unsets mapping between data and CacheObject instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d) = 0;

        /**
         * \brief Unsets dirty data mapping between data and CacheObject
         * instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetDirtyData(Data *d) = 0;

        /**
         * \brief Acquires access to CacheObject instance
         * \param access can be cache::EXCLUSIVE or cache::SHARED
         */
        void AcquireAccess(cache::CacheAccess access);

        /**
         * \brief Releases access to CacheObject instance
         */
        void ReleaseAccessInternal();

        /**
         * \brief Checks if CacheObject instance is available for use in access
         * mode
         * \param access denotes the mode (cache::EXCLUSIVE/ cache::SHARED) that application
         * wants
         * \return Boolean indicating if the instance is available
         */
        bool IsAvailable(cache::CacheAccess access) const;

        /**
         * \brief Accessor for id_ member
         * \return Instance's prototype id, type co_id_t
         */
        cache::co_id_t id() const;

        /**
         * \brief Accessor for object_region_ member
         * \return Instance's object_region_, of type GeometricRegion
         */
        GeometricRegion object_region() const;

        /**
         * \brief Setter for object_region_ member
         * \param object_region is of type GeometricRegion
         */
        void set_object_region(const GeometricRegion &object_region);

        bool pending_flag() {
          return pending_flag_;
        }
        void set_pending_flag() {
          pending_flag_ = true;
        }
        void unset_pending_flag() {
          pending_flag_ = false;
        }

        virtual size_t memory_size() {
          return sizeof(*this);
        }

        uint64_t unique_id() {
          return unique_id_;
        }

        void set_unique_id(const uint64_t unique_id) {
          unique_id_ = unique_id;
        }

        virtual std::string name() {
          return "cache_object";
        }

    private:
        uint64_t unique_id_;
        bool pending_flag_;
        /**
         * \brief Setter for id_ member
         * \param id, of type cache::co_id_t
         */
        void set_id(cache::co_id_t id);

        // prototype information
        static cache::co_id_t ids_allocated_;
        cache::co_id_t id_;

        // access information
        cache::CacheAccess access_;
        int users_;

    protected:
        // read/ write/ object region information
        GeometricRegion object_region_;
        GeometricRegion write_region_;
};  // class CacheObject

typedef std::vector<CacheObject *> CacheObjects;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_CACHE_OBJECT_H_
