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
  * Helper functions for escaping and unescaping.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "shared/escaper.h"

namespace nimbus {

bool IsEmptyString(std::string str) {
  if (str.length() == 0) {
      return true;
  }
  return false;
}

void EscapeString(std::string* input) {
  boost::algorithm::replace_all(*input, "%", "%0");
  boost::algorithm::replace_all(*input, ";", "%1");
  boost::algorithm::replace_all(*input, " ", "%2");
  boost::algorithm::replace_all(*input, "\n", "%3");
  boost::algorithm::replace_all(*input, "\t", "%4");
  boost::algorithm::replace_all(*input, "\r", "%5");
  boost::algorithm::replace_all(*input, ",", "%6");
}

void UnescapeString(std::string* input) {
  boost::algorithm::replace_all(*input, "%1", ";");
  boost::algorithm::replace_all(*input, "%2", " ");
  boost::algorithm::replace_all(*input, "%3", "\n");
  boost::algorithm::replace_all(*input, "%4", "\t");
  boost::algorithm::replace_all(*input, "%5", "\r");
  boost::algorithm::replace_all(*input, "%6", ",");
  boost::algorithm::replace_all(*input, "%0", "%");
}

}  // namespace nimbus


