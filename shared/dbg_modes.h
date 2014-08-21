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
 *   FILE: dbg_modes.h 
 * AUTHOR: Phil Levis (pal)
 *  DESCR: Definition of dbg modes and the bindings to DBG env settings. 
 */

#ifndef NIMBUS_SHARED_DBG_MODES_H_
#define NIMBUS_SHARED_DBG_MODES_H_

typedef uint64_t nimbus_dbg_mode;

#define DBG_MODE(x)     (1ULL << (x))

enum {
  DBG_ALL =             (~0ULL),        /* umm, "verbose"               */

/*====== Core modes =============*/
  DBG_BOOT =            DBG_MODE(0),    /* the boot sequence            */
  DBG_NET =             DBG_MODE(1),    /* network - connections, etc.  */
  DBG_PARSE =           DBG_MODE(2),    /* parsing packets              */
  DBG_SCHED =           DBG_MODE(3),    /* scheduler                    */
  DBG_DATA =            DBG_MODE(4),    /* data movement                */
  DBG_WORKER =          DBG_MODE(5),    /* worker execution             */
  DBG_HOSTS =           DBG_MODE(6),    /* host map/state updates       */
  DBG_DATA_OBJECTS =    DBG_MODE(7),    /* data manager/objects         */
  DBG_MEMORY =          DBG_MODE(8),    /* memory allocation            */
  DBG_TRANSLATE =       DBG_MODE(9),    /* data translation             */
  DBG_WORKER_FD =       DBG_MODE(10),   /* worker execution(frontend)   */
  DBG_WORKER_BD =       DBG_MODE(11),   /* worker execution(backend)    */
  DBG_APP       =       DBG_MODE(12),   /* worker execution(backend)    */
/*====== For application use =========*/
  DBG_PROJ =            DBG_MODE(57),   /* Projection module  -quh      */
  DBG_USR1 =            DBG_MODE(58),   /* User component 1             */
  DBG_USR2 =            DBG_MODE(59),   /* User component 2             */
  DBG_USR3 =            DBG_MODE(60),   /* User component 3             */
  DBG_TEMP =            DBG_MODE(61),   /* Temorpary testing use        */

  DBG_WARN =           DBG_MODE(62),   /* Warning condition              */
  DBG_ERROR =           DBG_MODE(63),   /* Error condition              */
  DBG_NONE =            0,              /* Nothing                      */

  DBG_DEFAULT =      DBG_ALL            /* default modes, 0 for none    */
};

#define DBG_NAMETAB \
        {"all",     DBG_ALL}, \
        {"boot",    DBG_BOOT   | DBG_ERROR | DBG_WARN}, \
        {"net",     DBG_NET    | DBG_ERROR | DBG_WARN}, \
        {"parse",   DBG_PARSE  | DBG_ERROR | DBG_WARN}, \
        {"sched",   DBG_SCHED  | DBG_ERROR | DBG_WARN}, \
        {"data",    DBG_DATA   | DBG_ERROR | DBG_WARN}, \
        {"worker",  DBG_WORKER | DBG_ERROR | DBG_WARN}, \
        {"hosts",   DBG_HOSTS  | DBG_ERROR | DBG_WARN}, \
        {"dobjects", DBG_DATA_OBJECTS | DBG_ERROR | DBG_WARN},  \
        {"memory",  DBG_MEMORY | DBG_ERROR | DBG_WARN},       \
        {"translate",  DBG_TRANSLATE | DBG_ERROR | DBG_WARN},       \
        {"app",     DBG_APP    | DBG_ERROR | DBG_WARN}, \
\
        {"proj",    DBG_PROJ | DBG_ERROR | DBG_WARN}, \
        {"usr1",    DBG_USR1 | DBG_ERROR | DBG_WARN}, \
        {"usr2",    DBG_USR2 | DBG_ERROR | DBG_WARN}, \
        {"usr3",    DBG_USR3 | DBG_ERROR | DBG_WARN}, \
        {"temp",    DBG_TEMP | DBG_ERROR | DBG_WARN}, \
        {"warn",    DBG_WARN | DBG_ERROR}, \
        {"error",   DBG_ERROR}, \
\
        {"none", DBG_NONE}, \
        { NULL, DBG_ERROR }

#define DBG_ENV         "DBG"

#endif  // NIMBUS_SHARED_DBG_MODES_H_
