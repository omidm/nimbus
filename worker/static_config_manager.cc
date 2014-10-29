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
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include "shared/dbg.h"
#include "worker/static_config_manager.h"

namespace nimbus {

StaticConfigManager::StaticConfigManager() {
  pthread_mutex_init(&internal_lock_, NULL);
}

StaticConfigManager::~StaticConfigManager() {
  pthread_mutex_destroy(&internal_lock_);
}

// Registration must happen before everything starts.
void StaticConfigManager::RegisterPrototype(
    const std::string& config_name, StaticConfigVariable* config_variable) {
  PrototypeMap::iterator iter = prototype_map_.find(config_name);
  if (iter != prototype_map_.end()) {
    dbg(DBG_ERROR, "Prototype of static config variable %s is already "
        "registered.", config_name.c_str());
    return;
  }
  prototype_map_[config_name] = config_variable;
  static_config_pool_map_[config_name] = new StaticConfigPool;
}

StaticConfigVariable* StaticConfigManager::GetStaticConfigVariable(
      const std::string& config_name, const GeometricRegion& local_region) {
  PrototypeMap::iterator prototype_iter = prototype_map_.find(config_name);
  if (prototype_iter == prototype_map_.end()) {
    dbg(DBG_ERROR, "Prototype of static config variable %s is not found.",
        config_name.c_str());
    return NULL;
  }
  StaticConfigPoolMap::iterator pool_iter =
      static_config_pool_map_.find(config_name);
  if (pool_iter == static_config_pool_map_.end()) {
    dbg(DBG_ERROR, "Config pool of static config variable %s is not found.",
        config_name.c_str());
    return NULL;
  }
  StaticConfigPool* config_pool = pool_iter->second;

  pthread_mutex_lock(&internal_lock_);
  if (config_pool->find(local_region) == config_pool->end()) {
    (*config_pool)[local_region] =
        prototype_iter->second->CreateNew(local_region);
  }
  StaticConfigVariable* result = (*config_pool)[local_region];
  result->ChangeUserCount(1);
  pthread_mutex_unlock(&internal_lock_);

  return result;
}

void StaticConfigManager::ReleasestaticConfigVariable(
    StaticConfigVariable* config_variable) {
  pthread_mutex_lock(&internal_lock_);
  config_variable->ChangeUserCount(-1);
  pthread_mutex_unlock(&internal_lock_);
}

}  // namespace nimbus
