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
  * Job Ancestor Chain.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_ANCESTOR_CHAIN_H_
#define NIMBUS_SCHEDULER_ANCESTOR_CHAIN_H_

#include <boost/shared_ptr.hpp>
#include <vector>
#include <set>
#include <list>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "scheduler/ancestor_entry.h"

namespace nimbus {

class AncestorChain {
  public:
    typedef std::list<AncestorEntry> Pool;
    typedef std::list<Pool> Chain;

    AncestorChain();
    AncestorChain(const AncestorChain& other);
    virtual ~AncestorChain();

    Chain chain() const;
    const Chain* chain_p() const;

    void set_chain(Chain chain);

    static bool MergeAncestorChains(
        std::vector<boost::shared_ptr<const AncestorChain> > chains,
        boost::shared_ptr<AncestorChain> *result);

    bool LookUpVersion(logical_data_id_t l_id, data_version_t* version);

    AncestorChain& operator=(const AncestorChain& right);

  private:
    Chain chain_;

    static bool MergeTwoAncestorChains(
        boost::shared_ptr<const AncestorChain> c_1,
        boost::shared_ptr<const AncestorChain> c_2,
        boost::shared_ptr<AncestorChain> *result);
};


}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_ANCESTOR_CHAIN_H_
