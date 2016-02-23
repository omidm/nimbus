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

# Starts Nimbus controller on the machine this script is executed on.

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

function print_usage {
  echo -e "${Cya}Runs fault tolerance tests against stencil_1d application."
  echo -e "\nUsage:"
  echo -e "./scripts/test-stencil-ft.sh"
  echo -e "${RCol}"
}

ITER_NUM=500
CHUNK_NUM=16
BATCH_NUM=6
THREAD_NUM=16
TIME_OUT_T=30
BASE_LENGTH=20
FT_PERIOD=2
FAULT_DELAY=5
CHKP_NUM=2 # FAULT_DELAY / FT_PERIOD

CONTROLLER_ARGS="-t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1"
CONTROLLER_ARGS_FT="--aft --ft_period ${FT_PERIOD} -t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1"
WORKER_ARGS="4 --othread ${THREAD_NUM}"
APPLICATION_ARGS="-i ${ITER_NUM} --pn ${CHUNK_NUM} --cpp 1"



if [ -z "${NIMBUS_HOME}" ]; then
  export NIMBUS_HOME="$(cd "`dirname "$0"`"/..; pwd)"
fi

if [ -z "${DBG}" ]; then
  export DBG="error"
fi

if [ -z "${TTIMER}" ]; then
  export TTIMER="l1"
fi


# start_experiment controller_args worker_args app_args
function start_experiment {
  make ${NIMBUS_HOME}/ clean-logs &> /dev/null
  ${NIMBUS_HOME}/scripts/stop-workers.sh &> /dev/null
  ${NIMBUS_HOME}/scripts/stop-controller.sh &> /dev/null
  ${NIMBUS_HOME}/scripts/start-controller.sh $1 &> /dev/null
  ${NIMBUS_HOME}/scripts/start-workers.sh $2 -l applications/simple/stencil_1d/libstencil_1d.so $3 &> /dev/null
}

# wait_to_finish
function wait_to_finish {
  start_time=$(date +%s)
  end_time=$(date +%s)
  progress_bar="waiting ..."
  success=0
  while [ "$((${end_time}-${start_time}))" -lt "${TIME_OUT_T}" ]; do
    success=$(cat ${NIMBUS_HOME}/logs/controller/stdout | grep -c "Simulation Terminated")
    if [ ${success} == "1" ]; then
      end_time=$(date +%s)
      break
    else
      end_time=$(date +%s)
      echo -ne "${Yel}${progress_bar} \r${RCol}"
      progress_bar=${progress_bar}"."
      sleep 1
    fi
  done
  
  if [ ${success} != "1" ]; then
    echo -e "${Red}[ TIMEOUT ] experiment did not finish before time out!${RCol}"
    exit 1
  fi
}

function get_final_hash {
  sync
  local final_hash=$(grep "FINAL HASH" logs/workers/*/stdout  | sed 's/.*FINAL HASH: //')
  echo "${final_hash}"
}

# check_hash correct_hash
function check_hash {
  sync
  local final_hash=$(grep "FINAL HASH" logs/workers/*/stdout  | sed 's/.*FINAL HASH: //')
  if [ "${final_hash}" == "$1" ]; then
    echo -e "${Gre}[ SUCCESS ] hash value matches!${RCol}"
  else
    echo -e "${Red}[ FAILED  ] hash value does not match! [${final_hash} != $1]${RCol}"
    exit 1
  fi
}

function kill_one_worker {
  WORKER_PIDS=$(ps -fu $USER| grep "nimbus_worker" | grep -v "grep" | awk '{print $2}')
  for pid in ${WORKER_PIDS}
  do
    echo -e "${Yel}    killing nimbus worker with pid: ${pid}${RCol}"
    kill ${pid}
    break
  done
}

function pause_controller {
  CONTROLLER_PID=$(ps -fu $USER| grep "nimbus_controller" | grep -v "grep" | awk '{print $2}')
  echo -e "${Yel}    pausing nimbus controller with pid: ${CONTROLLER_PID}${RCol}"
  kill -STOP ${CONTROLLER_PID}
}

function resume_controller {
  CONTROLLER_PID=$(ps -fu $USER| grep "nimbus_controller" | grep -v "grep" | awk '{print $2}')
  echo -e "${Yel}    resuming nimbus controller with pid: ${CONTROLLER_PID}${RCol}"
  kill -CONT ${CONTROLLER_PID}
}


# check_checkpoints min_count
function check_checkpoints {
  local count=$(cat ${NIMBUS_HOME}/logs/controller/stdout | grep -c "Checkpoint")
  if [ ${count} -lt "$1" ]; then
    echo -e "${Red}[ FAILED  ] controller did not make enough checkpoints, only ${count}!${RCol}"
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



if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  print_usage
  exit 0
fi


echo -e "${Cya}Ruunnig the base experiment without faults:${RCol}"
while true; do
  start_experiment "${CONTROLLER_ARGS}" "${WORKER_ARGS}" "${APPLICATION_ARGS}"
  wait_to_finish
 
  ELAPSED=$((${end_time}-${start_time}))
  NORM_ELAPSED=${ELAPSED}
  if [ "${ELAPSED}" -gt "4" ]; then
    NORM_ELAPSED=$((${ELAPSED}-4))
  fi
  if [ ${ELAPSED} -lt "${BASE_LENGTH}" ]; then
    echo -e "${Yel}[ WARNING ] base experiment with ${ITER_NUM} iterations was only ${ELAPSED} (<${BASE_LENGTH}) seconds!${RCol}"
    ITER_NUM=$((((${BASE_LENGTH}-${ELAPSED})*(${ITER_NUM}/${NORM_ELAPSED}))+${ITER_NUM}))
    echo -e "${Yel}[ WARNING ] retrying the experiment with ${ITER_NUM} iterations to get the minimum length ... ${RCol}"
    APPLICATION_ARGS="-i ${ITER_NUM} --pn ${CHUNK_NUM} --cpp 1"
  else
    echo -e "${Gre}[ SUCCESS ] base experiment finished in ${ELAPSED} seconds.${RCol}"
    break
  fi
done
base_hash=$(get_final_hash)
if ! [ -z ${base_hash} ]; then
  echo -e "${Gre}[ SUCCESS ] got the hash value of ${base_hash}.${RCol}"
else
  echo -e "${Red}[ FAILED  ] could not get the hash value!.${RCol}"
  exit 1
fi


echo -e "${Cya}Ruunnig the experiment with faults:${RCol}"
start_experiment "${CONTROLLER_ARGS_FT}" "${WORKER_ARGS}" "${APPLICATION_ARGS}"

echo -e "${Yel}    waiting ${FAULT_DELAY} seconds before killing first worker...${RCol}"
sleep ${FAULT_DELAY}
pause_controller
check_checkpoints ${CHKP_NUM}
kill_one_worker
resume_controller

echo -e "${Yel}    waiting ${FAULT_DELAY} seconds before killing second worker...${RCol}"
sleep ${FAULT_DELAY}
pause_controller
check_checkpoints $((2*${CHKP_NUM}))
kill_one_worker
resume_controller

echo -e "${Yel}    waiting ${FAULT_DELAY} seconds before killing third worker...${RCol}"
sleep ${FAULT_DELAY}
pause_controller
check_checkpoints $((3*${CHKP_NUM}))
kill_one_worker
resume_controller

wait_to_finish
ELAPSED=$((${end_time}-${start_time}+3*${FAULT_DELAY}))
echo -e "${Gre}[ SUCCESS ] experiment with fault finished in ${ELAPSED} seconds.${RCol}"
check_rewinding 3
check_hash "${base_hash}"



make ${NIMBUS_HOME}/ clean-logs &> /dev/null
echo -e "${Gre}\n[ PASSED  ] stencil fault tolerance test passed successfuly!${RCol}"
exit 0


