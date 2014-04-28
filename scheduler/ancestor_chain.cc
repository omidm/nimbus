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

#include "scheduler/ancestor_chain.h"

using namespace nimbus; // NOLINT

AncestorChain::AncestorChain() {
}

AncestorChain::AncestorChain(const AncestorChain& other) {
  chain_ = other.chain_;
}

AncestorChain::~AncestorChain() {
}

AncestorChain::Chain AncestorChain::chain() const {
  return chain_;
}

const AncestorChain::Chain* AncestorChain::chain_p() const {
  return &chain_;
}

void AncestorChain::set_chain(Chain chain) {
  chain_ = chain;
}

AncestorChain& AncestorChain::operator=(const AncestorChain& right) {
  chain_ = right.chain_;
  return (*this);
}

bool AncestorChain::MergeAncestorChains(
    std::vector<boost::shared_ptr<const AncestorChain> > chains,
    boost::shared_ptr<AncestorChain> *result) {
  size_t count = chains.size();

  if (count == 0) {
    return false;
  }

  if (count == 1) {
    *result = boost::shared_ptr<AncestorChain>(new AncestorChain());
    (*result)->set_chain(chains[0]->chain());
    return true;
  }

  std::vector<boost::shared_ptr<const AncestorChain> > reduced;
  if ((count % 2) == 1) {
    reduced.push_back(chains[count - 1]);
    --count;
  }

  for (size_t i = 0; i < count; i = i + 2) {
    boost::shared_ptr<AncestorChain> merged;
    if (MergeTwoAncestorChains(chains[i], chains[i+1], &merged)) {
      reduced.push_back(merged);
    } else {
      return false;
    }
  }

  return MergeAncestorChains(reduced, result);
}


bool AncestorChain::MergeTwoAncestorChains(
    boost::shared_ptr<const AncestorChain> c_1,
    boost::shared_ptr<const AncestorChain> c_2,
    boost::shared_ptr<AncestorChain> *result) {
  boost::shared_ptr<AncestorChain> merged(new AncestorChain());

  const Chain* shallower = c_1->chain_p();
  const Chain* deeper = c_2->chain_p();
  if (c_1->chain_p()->size() > c_2->chain_p()->size()) {
    shallower = c_2->chain_p();
    deeper = c_1->chain_p();
  }

  Chain chain;
  // std::set<job_id_t> already_in_chain;
  boost::unordered_set<job_id_t> already_in_chain;

  size_t diff = deeper->size() - shallower->size();
  Chain::const_reverse_iterator d_iter = deeper->rbegin();
  Chain::const_reverse_iterator s_iter = shallower->rbegin();

  for (size_t i = 0; i < diff; i++) {
    chain.push_front(*d_iter);
    Pool::const_iterator it;
    for (it = d_iter->begin(); it != d_iter->end(); ++it) {
      already_in_chain.insert(it->id());
    }
    ++d_iter;
  }

  while ((d_iter != deeper->rend()) && (s_iter != shallower->rend())) {
    Pool pool;

    Pool::const_iterator it;
    for (it = d_iter->begin(); it != d_iter->end(); ++it) {
      if (!already_in_chain.count(it->id())) {
        pool.push_back(*it);
        already_in_chain.insert(it->id());
      }
    }
    for (it = s_iter->begin(); it != s_iter->end(); ++it) {
      if (!already_in_chain.count(it->id())) {
        pool.push_back(*it);
        already_in_chain.insert(it->id());
      }
    }

    chain.push_front(pool);

    ++d_iter;
    ++s_iter;
  }

  merged->set_chain(chain);
  *result  = merged;

  return true;
}

bool AncestorChain::AncestorChain::LookUpVersion(
    logical_data_id_t l_id,
    data_version_t *version) {
  bool found = false;
  Chain::iterator iter;
  for (iter = chain_.begin(); iter != chain_.end(); ++iter) {
    Pool pool = *iter;
    Pool::iterator it;
    for (it = pool.begin(); it != pool.end(); ++it) {
      if (it->version_map()->query_entry(l_id, version)) {
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }

  return found;
}


