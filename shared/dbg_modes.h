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
/*====== For application use =========*/
  DBG_USR1 =            DBG_MODE(59),   /* User component 1             */
  DBG_USR2 =            DBG_MODE(60),   /* User component 2             */
  DBG_USR3 =            DBG_MODE(61),   /* User component 3             */
  DBG_TEMP =            DBG_MODE(62),   /* Temorpary testing use        */

  DBG_ERROR =           DBG_MODE(63),   /* Error condition              */
  DBG_NONE =            0,              /* Nothing                      */

  DBG_DEFAULT =      DBG_ALL            /* default modes, 0 for none    */
};

#define DBG_NAMETAB \
        {"all",     DBG_ALL}, \
        {"boot",    DBG_BOOT   | DBG_ERROR}, \
        {"net",     DBG_NET    | DBG_ERROR}, \
        {"parse",   DBG_PARSE  | DBG_ERROR}, \
        {"sched",   DBG_SCHED  | DBG_ERROR}, \
        {"data",    DBG_DATA   | DBG_ERROR}, \
        {"worker",  DBG_WORKER | DBG_ERROR}, \
        {"hosts",   DBG_HOSTS  | DBG_ERROR}, \
\
        {"usr1",    DBG_USR1 | DBG_ERROR}, \
        {"usr2",    DBG_USR2 | DBG_ERROR}, \
        {"usr3",    DBG_USR3 | DBG_ERROR}, \
        {"temp",    DBG_TEMP | DBG_ERROR}, \
        {"error",   DBG_ERROR}, \
\
        {"none", DBG_NONE}, \
        { NULL, DBG_ERROR }

#define DBG_ENV         "DBG"

#endif  // NIMBUS_SHARED_DBG_MODES_H_
