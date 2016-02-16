#!/usr/bin/env bash

#
#  Copyright 2013 Stanford University.
#  All rights reserved.
# 
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions
#  are met:
# 
#  - Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 
#  - Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
# 
#  - Neither the name of the copyright holders nor the names of
#    its contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
# 
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#  FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
#  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
#  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
#  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
#  OF THE POSSIBILITY OF SUCH DAMAGE.
#

# Starts Nimbus workers on the machine this script is executed on.

# Author: Omid Mashayekhi <omidm@stanford.edu>

# **************************
# Text Reset
RCol='\x1B[0m'    
# Regular           
Bla='\x1B[0;30m';
Red='\x1B[0;31m';
Gre='\x1B[0;32m';
Yel='\x1B[0;33m';
Blu='\x1B[0;34m';
Pur='\x1B[0;35m';
Cya='\x1B[0;36m';
Whi='\x1B[0;37m';
# **************************

function print_usage {
  echo -e "${Blu}Usage:"
  echo -e "${Blu}./scripts/start-workers.sh"
  echo -e "                    [number-of-workers-to-launch]"
  echo -e "                    --flush [to redirect stdout/stderr to current console]"
  echo -e "                    ... <worker options> ... "
  echo -e "                    -l/--app_lib [path/to/application-bundle]"
  echo -e "                    ... <application options> ... "
  cd ${WORKER_DIR}; "./${WORKER_BIN}" -h 2>&1
  echo -e ">> worker listening ports are set to (--port/-p arg)+n, where n is [1 to worker_num]."
  echo -e "   default --port/-p is 5900, so default worker ports are 5901, 5902, ..."
  echo -e ">> controller listening port is set to 5900 by default."
  echo -e ">> controller ip is set to \"localhost\" by default."
  echo -e "\n>> to get the application options pass -h, after -l/--app_lib option."
  echo -e "     and for simplicity also use --flush option to dump the help in console."
  echo -e "${RCol}"
}


if [ -z "${NIMBUS_HOME}" ]; then
  export NIMBUS_HOME="$(cd "`dirname "$0"`"/..; pwd)"
fi

if [ -z "${DBG}" ]; then
  export DBG="error"
fi

if [ -z "${TTIMER}" ]; then
  export TTIMER="l1"
fi

WORKER_DIR="${NIMBUS_HOME}/nodes/nimbus_worker/"
WORKER_BIN="nimbus_worker"
LOG_DIR="${NIMBUS_HOME}/logs/workers"
LOG_DIR_OLD="${NIMBUS_HOME}/logs/workers-old"

ARGS="$@"
re='^[0-9]+$'

if ! [[ "$1" =~ $re ]]; then
  if [[ ${ARGS} = *--help* ]] || [[ ${ARGS} = *-h* ]]; then
    echo -e "${Blu}Launches the nimbus workers on the machine this script is executed on.${RCol}"
    print_usage
    exit 0
  else
    echo -e "${Red}ERROR: the first argument should be the number of workers to luanch!${RCol}"
    print_usage
    exit 1
  fi
fi

WORKER_NUM="$1"
shift
ARGS="$@"

if [[ ${ARGS} != *--app_lib* ]] && [[ ${ARGS} != *-l* ]]; then
  if [[ ${ARGS} = *--help* ]] || [[ ${ARGS} = *-h* ]]; then
    echo -e "${Blu}Launches the nimbus workers on the machine this script is executed on.${RCol}"
    print_usage
    exit 0
  else
    echo -e "${Red}ERROR: need to pass the application bundle through --app_lib(-l) option!${RCol}"
    print_usage
    exit 1
  fi
fi
# the args befor/after --applib(-l) are worker/application args.
pivot_idx=0
for arg in ${ARGS}; do
  if [[ "-l" == "${arg}" ]] || [[ "--app_lib" == "${arg}" ]]; then
    break
  fi
  pivot_idx=$((${pivot_idx}+1))
done

NEW_ARGS=""

idx=0
controller_ip_given=false
for arg in ${ARGS}; do
  if [ "--cip" == "${arg}" ] && [ "${idx}" -lt "${pivot_idx}" ]; then
    controller_ip_given=true
    break
  fi
  idx=$((${idx}+1))
done
if [ "${controller_ip_given}" == "false" ]; then
  NEW_ARGS=" --cip localhost "${NEW_ARGS}
fi

idx=0
controller_port_given=false
for arg in ${ARGS}; do
  if [ "--cport" == "${arg}" ] && [ "${idx}" -lt "${pivot_idx}" ]; then
    controller_port_given=true
    break
  fi
  idx=$((${idx}+1))
done
if [ "${controller_port_given}" == "false" ]; then
  NEW_ARGS=" --cport 5900 "${NEW_ARGS}
fi


idx=0
FLUSH=false
PORT_NUM=5900
while (( "$#" )); do
  if ([ "--port" == "$1" ] || [ "-p" == "$1" ]) && [ "${idx}" -lt "${pivot_idx}" ]; then
    shift
    idx=$((${idx}+1))
    if ! [[ "$1" =~ $re ]]; then
      echo -e "${Red}ERROR: the port argument should be a number, you gave \"$1\"!${RCol}"
      print_usage
      exit 1
    else
      PORT_NUM=$1
    fi
    shift
    idx=$((${idx}+1))
  elif [ "--app_lib" == "$1" ] || [ "-l" == "$1" ]; then
    NEW_ARGS="${NEW_ARGS} $1"
    shift
    idx=$((${idx}+1))
    NEW_ARGS="${NEW_ARGS} ${NIMBUS_HOME}/$1"
    shift
    idx=$((${idx}+1))
  elif [ "--flush" == "$1" ] && [ "${idx}" -lt "${pivot_idx}" ]; then
    FLUSH=true
    if [ "${WORKER_NUM}" != "1" ]; then
      echo -e "${Red}ERROR: if the --flush is active you can only launch one worker!${RCol}"
      exit 1
    fi
    shift
    idx=$((${idx}+1))
  else
    NEW_ARGS="${NEW_ARGS} $1"
    shift
    idx=$((${idx}+1))
  fi
done

echo -e "${Blu}NIMBUS_HOME  ... \"${NIMBUS_HOME}\"${RCol}"
echo -e "${Blu}DBG  ........... \"${DBG}\"${RCol}"
echo -e "${Blu}TTIMER  ........ \"${TTIMER}\"${RCol}"


if [ "${FLUSH}" == "false" ]; then
  if [ -e "${LOG_DIR}" ]; then
    echo -e "${Yel}WARNING: found old worker log folder (only one older log is kept!)${RCol}"
    rm -rf "${LOG_DIR_OLD}"
    mv "${LOG_DIR}" "${LOG_DIR_OLD}"
  fi
  mkdir -p "${LOG_DIR}"
fi


for i in $(seq ${WORKER_NUM}); do
  mkdir -p "${LOG_DIR}/$i"
  W_ARGS="--port $((${PORT_NUM}+${i})) ${NEW_ARGS}"
  if [ "${FLUSH}" == "false" ]; then
    cd ${WORKER_DIR}; "./${WORKER_BIN}" ${W_ARGS} 1>"${LOG_DIR}/$i/stdout" 2>"${LOG_DIR}/$i/stderr" &
    echo -e "${Gre}Launched worker with arguments \"${W_ARGS}\"; find stdout/stderr at: ${LOG_DIR}/$i/.${RCol}"
  else
    echo -e "${Gre}Launching worker with arguments \"${W_ARGS}\".${RCol}"
    cd ${WORKER_DIR}; "./${WORKER_BIN}" ${W_ARGS} 2>&1
  fi
done

