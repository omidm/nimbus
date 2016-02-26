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

# Common utililty functions used in the test scripts.

# Author: Omid Mashayekhi <omidm@stanford.edu>

# **************************
# Text Reset
RCol='\x1B[0m'    
# Regular           
Bla='\x1B[0;90m';
Red='\x1B[0;91m';
Gre='\x1B[0;92m';
Yel='\x1B[0;93m';
Blu='\x1B[0;94m';
Pur='\x1B[0;95m';
Cya='\x1B[0;96m';
Whi='\x1B[0;97m';
# **************************

if [ -z "${NIMBUS_HOME}" ]; then
  export NIMBUS_HOME="$(cd "`dirname "$0"`"/..; pwd)"
fi

# print_help desc usage flag
function print_help {
  if [[ "$3" == "-h" ]] || [[ "$3" == "--help" ]]; then
    echo -e "${Cya}$1"
    echo -e "\nUsage:"
    echo -e "$2"
    echo -e "${RCol}"
    exit 0
  fi
}

# start_experiment controller_args worker_args app_lib app_args
function start_experiment {
  make ${NIMBUS_HOME}/ clean-logs &> /dev/null
  ${NIMBUS_HOME}/scripts/stop-workers.sh &> /dev/null
  ${NIMBUS_HOME}/scripts/stop-controller.sh &> /dev/null
  ${NIMBUS_HOME}/scripts/start-controller.sh $1 &> /dev/null
  ${NIMBUS_HOME}/scripts/start-workers.sh $2 -l $3 $4 &> /dev/null
}

# wait_with_bar seconds
function wait_with_bar {
  local progress_bar="waiting ..."
  for i in $(seq $1)
  do
    echo -ne "${Yel}${progress_bar} \r${RCol}"
    progress_bar=${progress_bar}"."
    sleep 1
  done
}

# count_doun seconds
function count_down {
  for i in $(seq $1)
  do
    echo -ne "${Yel}$(($1-${i}+1)) second(s) remaining ... \r${RCol}"
    sleep 1
  done
}

function basic_completion_check {
  local count=$(cat ${NIMBUS_HOME}/logs/controller/stdout | grep -c "Simulation Terminated")
  if [ ${count} == "1" ]; then
    echo "true"
  else
    echo "false"
  fi
}

# wait_to_succeed check_function time_out
function wait_to_succeed {
  local start_time=$(date +%s)
  local end_time=$(date +%s)
  local progress_bar="waiting ..."
  local success=false
  while [ "$((${end_time}-${start_time}))" -lt "$2" ]; do
    success=$($1)
    if [ ${success} == "true" ]; then
      end_time=$(date +%s)
      break
    else
      end_time=$(date +%s)
      echo -ne "${Yel}${progress_bar} \r${RCol}"
      progress_bar=${progress_bar}"."
      sleep 1
    fi
  done
  
  ELAPSED=$((${end_time}-${start_time}))
  if [ ${success} != "true" ]; then
    echo -e "${Red}[ TIMEOUT ] did not succeed before time out!${RCol}"
    exit 1
  fi
}

# get_hash tag
function get_hash {
  local pattern="s/.*$1: //"
  local final_hash=$(grep "$1" logs/workers/*/stdout  | sed "${pattern}")
  if ! [ -z ${final_hash} ]; then
    echo -e "${Gre}[ SUCCESS ] got the hash value of ${final_hash}.${RCol}"
    HASH=${final_hash}
  else
    echo -e "${Red}[ FAILED  ] could not get the hash value!.${RCol}"
    exit 1
  fi
}

# check_hash tag correct_hash
function check_hash {
  local pattern="s/.*$1: //"
  local final_hash=$(grep "$1" logs/workers/*/stdout  | sed "${pattern}")
  if ! [ -z ${final_hash} ]; then
    HASH=${final_hash}
    if [ "${HASH}" == "$2" ]; then
      echo -e "${Gre}[ SUCCESS ] hash value matches!${RCol}"
    else
      echo -e "${Red}[ FAILED  ] hash value does not match! [${HASH} != $2]${RCol}"
      exit 1
    fi
  else
    echo -e "${Red}[ FAILED  ] could not get the hash value!.${RCol}"
    exit 1
  fi
}

function kill_one_worker {
  local WORKER_PIDS=$(ps -fu $USER| grep "nimbus_worker" | grep -v "grep" | awk '{print $2}')
  for pid in ${WORKER_PIDS}
  do
    echo -e "${Cya}Killing nimbus worker with pid: ${pid}${RCol}"
    kill ${pid}
    return
  done
  echo -e "${Red}[ FAILED  ] could not kill any worker.${RCol}"
  exit 1
  
}

function pause_controller {
  local CONTROLLER_PID=$(ps -fu $USER| grep "nimbus_controller" | grep -v "grep" | awk '{print $2}')
  if ! [ -z ${CONTROLLER_PID} ]; then
    echo -e "${Cya}Pausing nimbus controller with pid: ${CONTROLLER_PID}${RCol}"
    kill -STOP ${CONTROLLER_PID}
  else
    echo -e "${Red}[ FAILED  ] could not pause the controller.${RCol}"
    exit 1
  fi
}

function resume_controller {
  local CONTROLLER_PID=$(ps -fu $USER| grep "nimbus_controller" | grep -v "grep" | awk '{print $2}')
  if ! [ -z ${CONTROLLER_PID} ]; then
    echo -e "${Cya}Resuming nimbus controller with pid: ${CONTROLLER_PID}${RCol}"
    kill -CONT ${CONTROLLER_PID}
  else
    echo -e "${Red}[ FAILED  ] could not resume the controller.${RCol}"
    exit 1
  fi
}

# check_checkpoints min_count
function check_checkpoints {
  local count=$(cat ${NIMBUS_HOME}/logs/controller/stdout | grep -c "Checkpoint")
  if [ ${count} -lt "$1" ]; then
    echo -e "${Red}[ FAILED  ] controller did not make enough checkpoints, only ${count}!${RCol}"
    exit 1
  else
    echo -e "${Gre}[ SUCCESS ] controller made ${count} checkpoints.${RCol}"
  fi
}

# check_rewinding count
function check_rewinding {
  local count=$(cat ${NIMBUS_HOME}/logs/controller/stdout | grep -c "Rewind")
  if [ ${count} == "$1" ]; then
    echo -e "${Gre}[ SUCCESS ] controller performd $1 rewindings.${RCol}"
  else
    echo -e "${Red}[ FAILED  ] controller did not perform enough rewinding, only ${count}!${RCol}"
    exit 1
  fi
}

function clean_logs {
  make ${NIMBUS_HOME}/ clean-logs &> /dev/null
}

