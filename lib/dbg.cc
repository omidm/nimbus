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
 *   FILE: dbg.c
 * AUTHOR: Philip Levis
 *  DESCR: Variables and initialization of DBG routines.
 */

#if !defined(_NIMBUS_NO_DBG)

#include "lib/dbg.h"

extern "C" {
static nimbus_dbg_mode_names dbg_nametab[] = {
  DBG_NAMETAB
};

nimbus_dbg_mode dbg_modes = 0;

void dbg(nimbus_dbg_mode mode, const char *format, ...) {
  if (dbg_active(mode)) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    va_end(args);
  }
}

bool dbg_active(nimbus_dbg_mode mode) {
  return (dbg_modes & mode) != 0;
}

void dbg_unset() {
  dbg_modes = 0;
}

void dbg_set(nimbus_dbg_mode modes) {
  dbg_modes = modes;
}

void dbg_add_mode(const char *name) {
  int cancel;
  nimbus_dbg_mode_names *mode;

  if (*name == '-') {
    cancel = 1;
    name++;
  } else {
    cancel = 0;
  }
  for (mode = dbg_nametab; mode->d_name != NULL; mode++) {
    if (strcmp(name, mode->d_name) == 0) {
      break;
    }
  }
  if (mode->d_name == NULL) {
    fprintf(stderr, "Warning: Unknown debug option: "
            "\"%s\"\n", name);
    return;
  }

  if (cancel) {
    dbg_modes &= ~mode->d_mode;
  } else {
    dbg_modes |= mode->d_mode;
  }
}

void dbg_add_modes(const char *modes) {
  char env[256];
  char *name;

  strncpy(env, modes, sizeof(env));
  for (name = strtok(env, ","); name; name = strtok(NULL, ","))  { // NOLINT
    dbg_add_mode(name);
  }
}

void dbg_init(void) {
  const char *dbg_env;

  dbg_modes = DBG_NONE;

  dbg_env = getenv(DBG_ENV);
  if (!dbg_env) {
    dbg_modes = DBG_DEFAULT;
    return;
  }

  dbg_add_modes(dbg_env);
}

void dbg_help(void) {
  int i = 0;
  printf("Known dbg modes: ");

  while (dbg_nametab[i].d_name != NULL) {
    printf("%s", dbg_nametab[i].d_name);
    if (dbg_nametab[i + 1].d_name != NULL) {
      printf(", ");
    }
    i++;
  }
  printf("\n");
}
}

#endif  // _NIMBUS_NO_DBG
