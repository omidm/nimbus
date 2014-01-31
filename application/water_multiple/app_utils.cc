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
#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/data_names.h"
#include "application/water_alternate_fine/reg_def.h"
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
                           nimbus::PdiVector *vec) {
        bool success = false;
        if (da.empty()) {
            return success;
        }
        std::set<Data *> ds;
        for (nimbus::DataArray::const_iterator it = da.begin(); it != da.end(); ++it) {
            Data *d = *it;
            if (d->name() == name) {
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

      job->GetCoveredLogicalObjects(&result, APP_FACE_VEL, &kRegGhostw3Inner[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_FACE_VEL_GHOST, &kRegGhostw3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_PHI, &kRegGhostw3Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_PRESSURE, &kRegGhostw1Outer[0]);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_POS_PARTICLES, &kDomainParticles);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_NEG_PARTICLES, &kDomainParticles);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_POS_REM_PARTICLES, &kDomainParticles);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_NEG_REM_PARTICLES, &kDomainParticles);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }

      job->GetCoveredLogicalObjects(&result, APP_LAST_UNIQUE_PARTICLE_ID, &kDomainParticles);
      for (size_t i = 0; i < result.size(); ++i) {
        read->insert(result[i]->id());
        write->insert(result[i]->id());
      }
    }

    bool SerializeParameter(const int frame, std::string* result){
      std::stringstream ss;
      ss << frame;
      *result = ss.str();
      return true;
    }

    bool SerializeParameter(const int frame, const T time, std::string *result) {
      std::stringstream ss;
      ss << frame;
      ss << "\n";
      ss << time;
      *result = ss.str();
      return true;
    }

    bool SerializeParameter(const int frame, const T time, const T dt, std::string *result) {
      std::stringstream ss;
      ss << frame;
      ss << "\n";
      ss << time;
      ss << "\n";
      ss << dt;
      *result = ss.str();
      return true;
    }

    bool LoadParameter(const std::string str, int* frame) {
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      return true;
    }

    bool LoadParameter(const std::string str, int* frame, T* time) {
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      ss >> (*time);
      return true;
    }

    bool LoadParameter(const std::string str, int* frame, T* time, T* dt) {
      std::stringstream ss;
      ss.str(str);
      ss >> (*frame);
      ss >> (*time);
      ss >> (*dt);
      return true;
    }

} // namespace application
