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

#include <set>
#include <string>
#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/reg_def.h"
#include "data/physbam/physbam_data.h"
#include "shared/logical_data_object.h"
#include "shared/nimbus.h"
#include "worker/physical_data_instance.h"

namespace application {

    typedef nimbus::Job Job;
    typedef nimbus::Data Data;
    typedef nimbus::DataArray DataArray;

    bool GetTranslatorData(const nimbus::Job *job,
                           const std::string &name,
                           const nimbus::DataArray& da,
                           nimbus::PdiVector *vec,
                           AccessType access_type) {
        bool success = false;
        if (da.empty()) {
            return success;
        }
        IDSet<physical_data_id_t> read_set = job->read_set();
        IDSet<physical_data_id_t> write_set = job->write_set();
        std::set<Data *> ds;
        for (nimbus::DataArray::const_iterator it = da.begin(); it != da.end(); ++it) {
            Data *d = *it;
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
                ds.insert(*it);
            }
        }
        if (ds.empty()) {
            return success;
        }
        for (std::set<Data *>::const_iterator it = ds.begin(); it != ds.end(); ++it) {
            Data *d = *it;
            std::string name_str = d->name();
            const nimbus::LogicalDataObject *ldo = job->GetLogicalObject(d->logical_id());
            nimbus::PhysicalDataInstance *pdi = new
                nimbus::PhysicalDataInstance(d->physical_id(),
                                             ldo, d,
                                             data_version_t(0));
            vec->push_back(pdi);
            success = true;
        }
        return success;
    }

    bool GetDataSet(const std::string &name,
                    const nimbus::DataArray &da,
                    std::set<Data * > *ds) {
        bool success = false;
        if (da.empty()) {
            return success;
        }
        for (nimbus::DataArray::const_iterator it = da.begin(); it != da.end(); ++it) {
            Data *d = *it;
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
            Data *d = *it;
            if (d->name() == name) {
                return d;
            }
        }
        return NULL;
    }

    Data* GetTheOnlyData(const nimbus::Job *job,
                         const std::string &name,
                         const nimbus::DataArray& da,
                         AccessType access_type) {
      if (da.empty()) {
        return NULL;
      }
      IDSet<physical_data_id_t> read_set = job->read_set();
      IDSet<physical_data_id_t> write_set = job->write_set();
      Data* result = NULL;
      for (nimbus::DataArray::const_iterator it = da.begin();
           it != da.end();
           ++it) {
        Data *d = *it;
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
        job->GetCoveredLogicalObjects(&result, arg, &region);
        for (size_t i = 0; i < result.size(); ++i) {
          set->insert(result[i]->id());
          dbg(APP_LOG, "Loaded logical id %d of variable %s to the set.\n", result[i]->id(), arg);
        }
        arg = va_arg(vl, char*);
      }
      va_end(vl);
    }

    // TODO: Get rid of these calls
    void LoadReadWriteSets(nimbus::Job* job,
        nimbus::IDSet<nimbus::logical_data_id_t>* read,
        nimbus::IDSet<nimbus::logical_data_id_t>* write) {
      nimbus::CLdoVector result;

      job->GetCoveredLogicalObjects(&result, APP_FACE_VEL, &kRegW3Central[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_FACE_VEL_GHOST, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_PHI, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_POS_PARTICLES, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_NEG_PARTICLES, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_POS_REM_PARTICLES, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_NEG_REM_PARTICLES, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_LAST_UNIQUE_PARTICLE_ID, &kRegW3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }
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
      int_dimension_t x, y, z, dx, dy, dz;
      ss >> x >> y >> z >> dx >> dy >> dz;
      region->Rebuild(x, y, z, dx, dy, dz);
      return true;
    }

    bool SerializeParameter(
        const int frame,
        const GeometricRegion& global_region,
        std::string* result) {
      std::stringstream ss;
      ss << frame;
      ss << "\n";
      ss << region_serial_helper(global_region);
      ss << "\n";
      *result = ss.str();
      return true;
    }

    bool SerializeParameter(
        const int frame,
        const T time,
        const GeometricRegion& global_region,
        std::string* result) {
      std::stringstream ss;
      ss << frame;
      ss << "\n";
      ss << time;
      ss << "\n";
      ss << region_serial_helper(global_region);
      ss << "\n";
      *result = ss.str();
      return true;
    }

    bool SerializeParameter(
        const int frame,
        const T time,
        const T dt,
        const GeometricRegion& global_region,
        const GeometricRegion& local_region,
        std::string *result) {
      std::stringstream ss;
      ss << frame;
      ss << "\n";
      ss << time;
      ss << "\n";
      ss << dt;
      ss << "\n";
      ss << region_serial_helper(global_region);
      ss << "\n";
      ss << region_serial_helper(local_region);
      ss << "\n";
      *result = ss.str();
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
      std::stringstream ss;
      ss << frame;
      ss << "\n";
      ss << time;
      ss << "\n";
      ss << dt;
      ss << "\n";
      ss << region_serial_helper(global_region);
      ss << "\n";
      ss << region_serial_helper(local_region);
      ss << "\n";
      ss << iteration;
      ss << "\n";
      *result = ss.str();
      return true;
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        GeometricRegion* global_region) {
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      if (global_region == NULL) {
        dbg(DBG_WARN, "Deserialization of job parameter might be wrong.\n");
      } else {
        std::string temp;
        std::getline(ss, temp);
        std::getline(ss, temp);
        region_deserial_helper(temp, global_region);
      }
      return true;
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        T* time,
        GeometricRegion* global_region) {
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      ss >> (*time);
      if (global_region == NULL) {
        dbg(DBG_WARN, "Deserialization of job parameter might be wrong.\n");
      } else {
        std::string temp;
        std::getline(ss, temp);
        std::getline(ss, temp);
        region_deserial_helper(temp, global_region);
      }
      return true;
    }

    bool LoadParameter(
        const std::string str,
        int* frame,
        T* time,
        T* dt,
        GeometricRegion* global_region,
        GeometricRegion* local_region) {
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      ss >> (*time);
      ss >> (*dt);
      if (global_region == NULL) {
        dbg(DBG_WARN, "Deserialization of job parameter might be wrong.\n");
      } else {
        std::string temp;
        std::getline(ss, temp);
        std::getline(ss, temp);
        region_deserial_helper(temp, global_region);
      }
      if (local_region == NULL) {
        dbg(DBG_WARN, "Deserialization of job parameter might be wrong.\n");
      } else {
        std::string temp;
        std::getline(ss, temp);
        region_deserial_helper(temp, local_region);
      }
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
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      ss >> (*time);
      ss >> (*dt);
      if (global_region == NULL) {
        dbg(DBG_WARN, "Deserialization of job parameter might be wrong.\n");
      } else {
        std::string temp;
        std::getline(ss, temp);
        std::getline(ss, temp);
        region_deserial_helper(temp, global_region);
      }
      if (local_region == NULL) {
        dbg(DBG_WARN, "Deserialization of job parameter might be wrong.\n");
      } else {
        std::string temp;
        std::getline(ss, temp);
        region_deserial_helper(temp, local_region);
      }
      ss >> (*iteration);
      return true;
    }

} // namespace application
