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
  * AfterMap data structure is meant to keep track of workers that needs to be
  * informes about successful finish of a specific job. In other words, if a
  * worker has to run a job say A, then worker has to be informed about all job
  * dones in the before set of A. So, each job done has to be communicated
  * with all the workers that are running it's after set.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_AFTER_MAP_H_
#define NIMBUS_SRC_SCHEDULER_AFTER_MAP_H_

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <list>
#include "src/shared/nimbus_types.h"
#include "src/scheduler/job_entry.h"
#include "src/scheduler/scheduler_worker.h"
#include "src/shared/dbg.h"

namespace nimbus {

class AfterMap {
  public:
    typedef boost::unordered_set<job_id_t> JobDonePool;
    typedef boost::unordered_set<SchedulerWorker*> Pool;
    typedef boost::unordered_map<job_id_t, Pool*> Map;
    typedef Map::iterator Iter;
    typedef Map::const_iterator ConstIter;

    AfterMap();
    virtual ~AfterMap();

    const Map* map() const;
    const Map* late_map() const;

    bool AddEntries(JobEntry* job);

    bool AddEntry(job_id_t job_id, SchedulerWorker *worker);

    bool GetWorkersWaitingOnJob(job_id_t job_id,
                                std::list<SchedulerWorker*> *list);

    bool NotifyJobDone(job_id_t job_id);

    bool PullLateMap(Map*& late_map);

    bool RemoveJobRecords(job_id_t job_id);

    bool RemoveWorkerRecords(SchedulerWorker *worker);

    bool JobIsDone(job_id_t job_id);

    static void DestroyMap(Map *map);

    void Clear();

  private:
    Map *map_;
    Map *late_map_;
    JobDonePool job_done_pool_;
    mutable boost::recursive_mutex mutex_;

    static void AddEntryToMap(Map *map, job_id_t job_id, SchedulerWorker *worker);

    static void RemoveJobRecordFromMap(Map *map, job_id_t job_id);

    static void RemoveWorkerRecordFromMap(Map *map, SchedulerWorker *worker);

    static bool EntryIsInMap(Map *map, job_id_t job_id, SchedulerWorker *worker);

    AfterMap(const AfterMap& other) {}

    AfterMap& operator=(const AfterMap& right) {return (*this);}
};


}  // namespace nimbus
#endif  // NIMBUS_SRC_SCHEDULER_AFTER_MAP_H_
