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
 *   FILE: dbg.h
 * AUTHOR: Philip Levis, original work by Mike Castelle
 *  DESCR: Run-time configuration of debug output
 *
 * Debug output determined by DBG environment variable. dbg_modes.h has
 * definitions of the settings possible. One can specify multiple debugging
 * outputs by comma-delimiting (e.g. DBG=sched,timer). Compiling with
 * NDEBUG defined (e.g. -DNDEBUG) will stop all of the debugging
 * output, will remove the debugging commands from the object file.
 *
 * example usage: dbg(DBG_TIMER, "timer went off at %d\n", time);
 *
 */

#ifndef NIMBUS_LIB_DBG_H_
#define NIMBUS_LIB_DBG_H_

#if !defined(_NIMBUS_NO_DBG)

#include <stdio.h>
#include <stdarg.h>
#include "lib/nimbus.h"
#include "lib/dbg_modes.h"

extern "C" {
typedef struct dbg_mode {
  const char* d_name;
  uint64_t d_mode;
} nimbus_dbg_mode_names;

void dbg(nimbus_dbg_mode mode, const char *format, ...);
bool dbg_active(nimbus_dbg_mode mode);
void dbg_add_mode(const char *mode);
void dbg_add_modes(const char *modes);
void dbg_init(void);
void dbg_help(void);
void dbg_unset();
void dbg_set(nimbus_dbg_mode);


#else
/* No debugging */

#define dbg(...) { }
#define dbg_add_mode(...) { }
#define dbg_add_modes(...) { }
#define dbg_init() { }
#define dbg_help() { }
#define dbg_active(x) (FALSE)

#endif
} // NOLINT
#endif  // NIMBUS_LIB_DBG_H_
