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
 * An AppObject is an application object corresponding to one/ multiple nimbus
 * variables.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#ifndef NIMBUS_DATA_APP_DATA_APP_OBJECT_H_
#define NIMBUS_DATA_APP_DATA_APP_OBJECT_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/app_data/app_data_defs.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class Data;
typedef std::vector<Data *> DataArray;
typedef std::set<Data *> DataSet;
class GeometricRegion;

/**
 * \class AppObject
 * \details Application object corresponding to one/ multiple nimbus variables.
 * AppVar and AppStruct, that inherit from AppObject, provide one variable and
 * multiple variable implementation respectively.
 */
class  AppObject {
    // list friend classes here
    friend class CacheManager;
    friend class CacheTable;
    friend class Data;

    public:
        /**
         * \brief Creates an AppObject
         * \return Constructed AppObject instance
         */
        explicit AppObject();

        /**
         * \brief Creates an AppObject
         * \param  ob_reg specifies application object region
         * \return Constructed AppObject instance
         */
        explicit AppObject(const GeometricRegion &ob_reg);

        /**
         * \brief Makes this instance a prototype. The application writer must
         * make a prototype for every application object he/ she plans to use.
         */
         void MakePrototype();

        /**
         * \brief Accessor for id_ member
         * \return Instance's prototype id, type co_id_t
         */
        app_data::ob_id_t id() const;

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

        uint64_t unique_id() {
          return unique_id_;
        }

        void set_unique_id(const uint64_t unique_id) {
          unique_id_ = unique_id;
        }

        std::string name() const {
          return name_;
        }
        void set_name(const std::string name) {
          name_ = name;
        }

        /**
         * \brief Dumps out data to a file, to be implemented by child classes
         * \param file_name is file to dump data to
         */
        virtual void DumpData(std::string file_name) {}

    private:
        std::string name_;
        uint64_t unique_id_;

        // prototype information
        static app_data::ob_id_t ids_allocated_;
        app_data::ob_id_t id_;

        // access information
        app_data::Access access_;
        int users_;

        void set_id(app_data::ob_id_t id);

        void AcquireAccess(app_data::Access access);
        void ReleaseAccessInternal();
        bool IsAvailable(app_data::Access access) const;

    protected:
        // read/ write/ object region information
        GeometricRegion object_region_;
        GeometricRegion write_region_;

        /**
         * \brief Pull data from AppObject to corresponding nimbus physical
         * data object - this is unsafe, should be used with appropriate checks
         * and flags set
         * \param d is data to flush to
         */
        virtual void PullData(Data *d) = 0;

        /**
         * \brief Unsets mapping between data and AppObject instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetData(Data *d) = 0;

        /**
         * \brief Unsets dirty data mapping between data and AppObject instance
         * \param d denotes the data to unmap
         */
        virtual void UnsetDirtyData(Data *d) = 0;

        // method for memory profiling
        virtual size_t memory_size() {
          return sizeof(*this);
        }
};  // class AppObject

typedef std::vector<AppObject *> AppObjects;

}  // namespace nimbus

#endif  // NIMBUS_DATA_APP_DATA_APP_OBJECT_H_
