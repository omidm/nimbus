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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include <limits>
#include <set>
#include <string>
#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/reg_def.h"
#include "data/physbam/physbam_data.h"
#include "data/physbam/protobuf_compiled/water_parameter.pb.h"
#include "shared/logical_data_object.h"
#include "shared/nimbus.h"
#include "worker/physical_data_instance.h"

namespace application {

    bool GetTranslatorData(const nimbus::Job *job,
                           const std::string &name,
                           const nimbus::DataArray& da,
                           nimbus::PdiVector *vec,
                           AccessType access_type) {
        bool success = false;
        if (da.empty()) {
            return success;
        }
        IDSet<nimbus::physical_data_id_t> read_set = job->read_set();
        IDSet<nimbus::physical_data_id_t> write_set = job->write_set();
        std::set<nimbus::Data *> ds;
        for (nimbus::DataArray::const_iterator it = da.begin(); it != da.end(); ++it) {
            nimbus::Data *d = *it;
            bool allowed = false;
            if ((access_type == READ_ACCESS) &&
                read_set.contains(d->physical_id())) {
              allowed = true;
            }
            if ((access_type == WRITE_ACCESS) &&
                write_set.contains(d->physical_id())) {
              allowed = true;
            }
            if (d->name() == name && allowed) {
                if (ds.count(*it) == 0) {
                  nimbus::Data *d = *it;
                  std::string name_str = d->name();
                  const nimbus::LogicalDataObject *ldo =
                      job->GetLogicalObject(d->logical_id());
                  nimbus::PhysicalDataInstance *pdi = new
                      nimbus::PhysicalDataInstance(d->physical_id(),
                                                   ldo, d,
                                                   nimbus::data_version_t(0));
                  vec->push_back(pdi);
                  ds.insert(*it);
                }
                success = true;
            }
        }
        /*
        if (ds.empty()) {
            return success;
        }
        for (std::set<nimbus::Data *>::const_iterator it = ds.begin(); it != ds.end(); ++it) {
            nimbus::Data *d = *it;
            std::string name_str = d->name();
            const nimbus::LogicalDataObject *ldo = job->GetLogicalObject(d->logical_id());
            nimbus::PhysicalDataInstance *pdi = new
                nimbus::PhysicalDataInstance(d->physical_id(),
                                             ldo, d,
                                             nimbus::data_version_t(0));
            vec->push_back(pdi);
            success = true;
        }
        */
        return success;
    }

    void GetReadData(const nimbus::Job *job,
                     const std::string &name,
                     const nimbus::DataArray &da,
                     nimbus::DataArray *read) {
        read->clear();
        if (da.empty())
            return;
        IDSet<nimbus::physical_data_id_t> read_set = job->read_set();
        size_t rs = read_set.size();
        for (size_t i = 0; i < rs; ++i) {
            Data *d = da[i];
            if (d->name() == name &&
                    read_set.contains(d->physical_id())) {
                read->push_back(d);
            }
        }
    }

    void GetWriteData(const nimbus::Job *job,
                      const std::string &name,
                      const nimbus::DataArray &da,
                      nimbus::DataArray *write) {
        write->clear();
        if (da.empty())
            return;
        IDSet<nimbus::physical_data_id_t> read_set = job->read_set();
        IDSet<nimbus::physical_data_id_t> write_set = job->write_set();
        size_t rs = read_set.size();
        size_t ws = write_set.size();
        for (size_t i = rs; i < rs + ws; ++i) {
            Data *d = da[i];
            if (d->name() == name &&
                    write_set.contains(d->physical_id())) {
                write->push_back(d);
            }
        }
    }

    bool GroupSyncData(const nimbus::Job *job,
                       const nimbus::DataArray &da,
                       DataVec *main_copy,
                       DataSetVec *scratch_copies) {
        if (da.size() < 2)
            return false;
        IDSet<nimbus::physical_data_id_t> read_set = job->read_set();
        IDSet<nimbus::physical_data_id_t> write_set = job->write_set();
        for (nimbus::DataArray::const_iterator it = da.begin();
                it != da.end(); it++) {
            nimbus::Data *d = *it;
            if (write_set.contains(d->physical_id()))
                main_copy->push_back(d);
        }
        for (size_t i = 0; i < main_copy->size(); i++) {
            nimbus::GeometricRegion region = main_copy->at(i)->region();
            DataVec *scratch = new DataVec();
            for (nimbus::DataArray::const_iterator it = da.begin();
                    it != da.end(); it++) {
                nimbus::Data *d = *it;
                if (region == d->region() && read_set.contains(d->physical_id()))
                    scratch->push_back(d);
            }
            scratch_copies->push_back(scratch);
        }
        return true;
    }

    bool GetDataSet(const std::string &name,
                    const nimbus::DataArray &da,
                    std::set<nimbus::Data * > *ds) {
        bool success = false;
        if (da.empty()) {
            return success;
        }
        for (nimbus::DataArray::const_iterator it = da.begin(); it != da.end(); ++it) {
            nimbus::Data *d = *it;
            if (d->name() == name) {
                ds->insert(*it);
                success = true;
            }
        }
        return success;
    }

    nimbus::Data* GetFirstData(const std::string &name,
                               const nimbus::DataArray &da) {
        if (da.empty()) {
            return NULL;
        }
        for (nimbus::DataArray::const_iterator it = da.begin(); it != da.end(); ++it) {
            nimbus::Data *d = *it;
            if (d->name() == name) {
                return d;
            }
        }
        return NULL;
    }

    nimbus::Data* GetTheOnlyData(const nimbus::Job *job,
                                 const std::string &name,
                                 const nimbus::DataArray& da,
                                 AccessType access_type) {
      if (da.empty()) {
        return NULL;
      }
      IDSet<nimbus::physical_data_id_t> read_set = job->read_set();
      IDSet<nimbus::physical_data_id_t> write_set = job->write_set();
      nimbus::Data* result = NULL;
      for (nimbus::DataArray::const_iterator it = da.begin();
           it != da.end();
           ++it) {
          nimbus::Data *d = *it;
        bool allowed = false;
        if ((access_type == READ_ACCESS) &&
            read_set.contains(d->physical_id())) {
              allowed = true;
        }
        if ((access_type == WRITE_ACCESS) &&
            write_set.contains(d->physical_id())) {
              allowed = true;
        }
        if (d->name() == name && allowed) {
          if (result == NULL) {
            result = d;
          } else if (result != d) {
            dbg(DBG_ERROR, "More than one physical data instances matches, "
                "but only one is expected.\n");
            // return NULL;
            // [TODO] Variable in read/write set will appear twice. Maybe it
            // will be changed as the asbtraction changes?  --quhang
            return result;
          }
        }
      }
      return result;
    }

    void DestroyTranslatorObjects(nimbus::PdiVector *vec) {
        if (vec->empty())
            return;
        for (nimbus::PdiVector::iterator it = vec->begin(); it != vec->end(); ++it) {
            delete *it;
        }
        vec->clear();
    }

    bool Contains(nimbus::IDSet<nimbus::logical_data_id_t> data_set,
                  nimbus::logical_data_id_t  id) {
        nimbus::IDSet<nimbus::logical_data_id_t>::IDSetIter it;
        for (it = data_set.begin(); it != data_set.end(); ++it) {
            if (*it == id) {
                return true;
            }
        }
        return false;
    }

    void LoadLogicalIdsInSet(nimbus::Job* job,
        nimbus::IDSet<nimbus::logical_data_id_t>* set,
        const nimbus::GeometricRegion& region, ...) {
      nimbus::CLdoVector result;
      va_list vl;
      va_start(vl, region);
      char* arg = va_arg(vl, char*);
      while (arg != NULL) {
        job->GetIntersectingLogicalObjects(&result, arg, &region);
        for (size_t i = 0; i < result.size(); ++i) {
          set->insert(result[i]->id());
        }
        arg = va_arg(vl, char*);
      }
      va_end(vl);
    }

    // TODO(quhang), this is only for temprory usage.
    std::string region_serial_helper(const GeometricRegion& region) {
      std::stringstream ss;
      ss << region.x() << " " << region.y() << " " << region.z() << " " <<
            region.dx() << " "  << region.dy() << " " << region.dz();
      return ss.str();
    }

    bool region_deserial_helper(
        const std::string input,
        GeometricRegion* region) {
      assert(region != NULL);
      std::stringstream ss(input);
      nimbus::int_dimension_t x, y, z, dx, dy, dz;
      ss >> x >> y >> z >> dx >> dy >> dz;
      region->Rebuild(x, y, z, dx, dy, dz);
      return true;
    }

    void SerializeRegionHelper(
        const GeometricRegion& region,
        nimbus_message::WaterParameter::GeometricRegion* output) {
      output->set_x(region.x());
      output->set_y(region.y());
      output->set_z(region.z());
      output->set_dx(region.dx());
      output->set_dy(region.dy());
      output->set_dz(region.dz());
    }

    bool SerializeParameter(
        const int frame,
        const GeometricRegion& global_region,
        std::string* result) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.set_frame(frame);
      SerializeRegionHelper(
          global_region,
          water_parameter.mutable_global_region());
      water_parameter.SerializeToString(result);
      return true;
    }

    bool SerializeParameter(
        const int frame,
        const T time,
        const GeometricRegion& global_region,
        std::string* result) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.set_frame(frame);
      water_parameter.set_time(time);
      SerializeRegionHelper(
          global_region,
          water_parameter.mutable_global_region());
      water_parameter.SerializeToString(result);
      return true;
    }

    bool SerializeParameter(
        const int frame,
        const T time,
        const T dt,
        const GeometricRegion& global_region,
        const GeometricRegion& local_region,
        std::string *result) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.set_frame(frame);
      water_parameter.set_time(time);
      water_parameter.set_dt(dt);
      SerializeRegionHelper(
          global_region,
          water_parameter.mutable_global_region());
      SerializeRegionHelper(
          local_region,
          water_parameter.mutable_local_region());
      water_parameter.SerializeToString(result);
      return true;
    }

    bool SerializeParameter(
        const int frame,
        const T time,
        const T dt,
        const GeometricRegion& global_region,
        const GeometricRegion& local_region,
        const int iteration,
        std::string *result) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.set_frame(frame);
      water_parameter.set_time(time);
      water_parameter.set_dt(dt);
      water_parameter.set_iteration(iteration);
      SerializeRegionHelper(
          global_region,
          water_parameter.mutable_global_region());
      SerializeRegionHelper(
          local_region,
          water_parameter.mutable_local_region());
      water_parameter.SerializeToString(result);
      return true;
    }

    void DeserializeRegionHelper(
        const nimbus_message::WaterParameter::GeometricRegion& region,
        GeometricRegion* output) {
      output->Rebuild(
          region.x(), region.y(), region.z(),
          region.dx(), region.dy(), region.dz());
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        GeometricRegion* global_region) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.ParseFromString(str);
      assert(water_parameter.has_frame());
      *frame = water_parameter.frame();
      assert(water_parameter.has_global_region());
      DeserializeRegionHelper(water_parameter.global_region(),
                              global_region);
      return true;
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        T* time,
        GeometricRegion* global_region) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.ParseFromString(str);
      assert(water_parameter.has_frame());
      *frame = water_parameter.frame();
      assert(water_parameter.has_time());
      *time = water_parameter.time();
      assert(water_parameter.has_global_region());
      DeserializeRegionHelper(water_parameter.global_region(),
                              global_region);
      return true;
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        T* time,
        T* dt,
        GeometricRegion* global_region,
        GeometricRegion* local_region) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.ParseFromString(str);
      assert(water_parameter.has_frame());
      *frame = water_parameter.frame();
      assert(water_parameter.has_time());
      *time = water_parameter.time();
      assert(water_parameter.has_dt());
      *dt = water_parameter.dt();
      assert(water_parameter.has_global_region());
      DeserializeRegionHelper(water_parameter.global_region(),
                              global_region);
      assert(water_parameter.has_local_region());
      DeserializeRegionHelper(water_parameter.local_region(),
                              local_region);
      return true;
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        T* time,
        T* dt,
        GeometricRegion* global_region,
        GeometricRegion* local_region,
        int* iteration) {
      nimbus_message::WaterParameter water_parameter;
      water_parameter.ParseFromString(str);
      assert(water_parameter.has_frame());
      *frame = water_parameter.frame();
      assert(water_parameter.has_time());
      *time = water_parameter.time();
      assert(water_parameter.has_dt());
      *dt = water_parameter.dt();
      assert(water_parameter.has_iteration());
      *iteration = water_parameter.iteration();
      assert(water_parameter.has_global_region());
      DeserializeRegionHelper(water_parameter.global_region(),
                              global_region);
      assert(water_parameter.has_local_region());
      DeserializeRegionHelper(water_parameter.local_region(),
                              local_region);
      return true;
    }

} // namespace application
