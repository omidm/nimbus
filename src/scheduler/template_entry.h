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
  * This is TemplateEntry module to hold and instantiate the templates.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_TEMPLATE_ENTRY_H_
#define NIMBUS_SRC_SCHEDULER_TEMPLATE_ENTRY_H_

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <set>
#include "src/shared/nimbus_types.h"
#include "src/shared/dbg.h"
#include "src/shared/log.h"
#include "src/shared/graph.h"
#include "src/scheduler/shadow_job_entry.h"
#include "src/scheduler/complex_job_entry.h"
#include "src/scheduler/template_job_entry.h"
#include "src/scheduler/binding_template.h"

namespace nimbus {

class JobManager;

class TemplateEntry {
  public:
    TemplateEntry(const std::string& template_name,
                  bool worker_template_active_,
                  bool mega_rcr_job_active_);
    ~TemplateEntry();

    bool finalized();
    size_t compute_jobs_num();
    std::string template_name();
    boost::shared_ptr<VersionMap> vmap_base() const;
    logical_data_id_t min_ldid() const;
    logical_data_id_t max_ldid() const;
    size_t ldid_count() const;

    virtual void set_vmap_base(boost::shared_ptr<VersionMap> vmap_base);

    bool Finalize();

    bool CleanPartiallyFilledTemplate();

    bool LoadBeforeSet(IDSet<job_id_t>* before_set,
                       const size_t& index,
                       const std::vector<job_id_t>& inner_job_ids,
                       const std::vector<job_id_t>& outer_job_ids);

    bool Instantiate(JobManager *job_manager,
                     const std::vector<job_id_t>& inner_job_ids,
                     const std::vector<job_id_t>& outer_job_ids,
                     const std::vector<Parameter>& parameters,
                     const job_id_t& parent_job_id);

    bool GetComplexJobEntry(ComplexJobEntry*& complex_job,
                            const job_id_t& job_id,
                            const job_id_t& parent_job_id,
                            const std::vector<job_id_t>& inner_job_ids,
                            const std::vector<job_id_t>& outer_job_ids,
                            const std::vector<Parameter>& parameters);

    TemplateJobEntry* AddComputeJob(const std::string& job_name,
                                    const job_id_t& job_id,
                                    const IDSet<logical_data_id_t>& read_set,
                                    const IDSet<logical_data_id_t>& write_set,
                                    const IDSet<logical_data_id_t>& scratch_set,
                                    const IDSet<logical_data_id_t>& reduce_set,
                                    const IDSet<job_id_t>& before_set,
                                    const IDSet<job_id_t>& after_set,
                                    const job_id_t& parent_job_id,
                                    const job_id_t& future_job_id,
                                    const bool& sterile,
                                    const GeometricRegion& region,
                                    const CombinerMap& combiner_map);

    bool AddExplicitCopyJob();

    size_t GetParentJobIndices(std::list<size_t>* list);

    TemplateJobEntry* GetJobAtIndex(size_t index);

    bool InitializeCursor(ComplexJobEntry::Cursor* cursor);

    bool AdvanceCursorForAssignment(ComplexJobEntry::Cursor* cursor);

    void AddToAccessPattern(const logical_data_id_t& ldid,
                        const data_version_t& diff_version,
                        const size_t& job_index);

    size_t QueryAccessPattern(const logical_data_id_t& ldid,
                              const data_version_t& diff_version,
                              std::list<size_t>* indices);

    bool AddBindingRecord(const ComplexJobEntry *complex_job,
                          BindingTemplate*& binding_template);

    bool QueryBindingRecord(const ComplexJobEntry *complex_job,
                            BindingTemplate*& binding_template);

    void MarkExplicitBinding(const ComplexJobEntry *complex_job);

    bool ExplicitBindingBefore(const ComplexJobEntry *complex_job);

  private:
    typedef std::vector<boost::shared_ptr<job_id_t> > PtrList;
    typedef boost::unordered_set<boost::shared_ptr<job_id_t> > PtrSet;
    typedef boost::unordered_map<job_id_t, boost::shared_ptr<job_id_t> > PtrMap;

    class TemplateComputeJobEntry {
      public:
        TemplateComputeJobEntry(const std::string& job_name,
                                boost::shared_ptr<job_id_t> job_id_ptr,
                                const IDSet<logical_data_id_t>& read_set,
                                const IDSet<logical_data_id_t>& write_set,
                                const IDSet<logical_data_id_t>& scratch_set,
                                const IDSet<logical_data_id_t>& reduce_set,
                                const PtrSet& before_set_ptrs,
                                const PtrSet& after_set_ptrs,
                                boost::shared_ptr<job_id_t> future_job_id_ptr,
                                const bool& sterile,
                                const GeometricRegion& region,
                                const CombinerMap& combiner_map)
          : job_name_(job_name),
            job_id_ptr_(job_id_ptr),
            read_set_(read_set),
            write_set_(write_set),
            scratch_set_(scratch_set),
            reduce_set_(reduce_set),
            before_set_ptrs_(before_set_ptrs),
            after_set_ptrs_(after_set_ptrs),
            future_job_id_ptr_(future_job_id_ptr),
            sterile_(sterile),
            region_(region),
            combiner_map_(combiner_map) {}

        ~TemplateComputeJobEntry() {}

        std::string job_name_;
        boost::shared_ptr<job_id_t> job_id_ptr_;
        IDSet<logical_data_id_t> read_set_;
        IDSet<logical_data_id_t> write_set_;
        IDSet<logical_data_id_t> scratch_set_;
        IDSet<logical_data_id_t> reduce_set_;
        PtrSet before_set_ptrs_;
        PtrSet after_set_ptrs_;
        boost::shared_ptr<job_id_t> future_job_id_ptr_;
        bool sterile_;
        GeometricRegion region_;
        CombinerMap combiner_map_;
    };
    typedef std::vector<TemplateComputeJobEntry> EntryList;
    PtrList job_id_ptrs_;
    PtrMap job_id_ptrs_map_;
    EntryList entry_list_;
    // TODO(omidm): currently we do not support future job id in templates!
    boost::shared_ptr<job_id_t> future_job_id_ptr_;


    typedef std::list<size_t> Bucket;
    typedef boost::unordered_map<data_version_t, Bucket*> VersionIndex;
    typedef boost::unordered_map<logical_data_id_t, VersionIndex*> AccessIndex;

    typedef boost::unordered_set<std::string> BindingSet;
    typedef boost::unordered_map<std::string, BindingTemplate*> BindingMap;

    std::string template_name_;
    bool worker_template_active_;
    bool mega_rcr_job_active_;

    AccessIndex access_pattern_;
    boost::mutex access_pattern_mutex_;

    BindingMap binding_records_;
    BindingSet explicit_binding_history_;

    bool finalized_;
    Graph<TemplateJobEntry, job_id_t> job_graph_;
    std::list<Vertex<TemplateJobEntry, job_id_t>*> traverse_queue_;
    TemplateJobEntryVector compute_jobs_;

    std::list<size_t> parent_job_indices_;
    std::vector<size_t> assign_ordered_indices_;
    boost::unordered_set<size_t> assign_batch_mark_indices_;
    size_t last_assign_index_;

    boost::shared_ptr<VersionMap> vmap_base_;

    logical_data_id_t min_ldid_;
    logical_data_id_t max_ldid_;
    size_t ldid_count_;

    bool AddTemplateJobEntryToJobGraph(TemplateJobEntry *job);

    void CompleteParentJobIndices();

    void CompleteBreadthFirstSearch();

    void CompleteLdidInfo();

    std::string ProduceBindingRecordName(const ComplexJobEntry *job);
};

}  // namespace nimbus
#endif  // NIMBUS_SRC_SCHEDULER_TEMPLATE_ENTRY_H_
