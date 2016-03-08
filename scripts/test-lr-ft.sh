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

# Runs fault tolerance test against stencil 1D application.

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

source scripts/test-utils.sh


DESC="Runs fault tolerance tests against logistic regression application."
USAGE="./scripts/test-stencil-ft.sh"

print_help "${DESC}" "${USAGE}" "$1"

ITER_NUM=500
SAMPLE_NUM_M=".001"
BATCH_NUM=6
THREAD_NUM=16
TIME_OUT_T=60
BASE_LENGTH=20
FT_PERIOD=2
FAULT_DELAY=5
CHKP_NUM=2 # FAULT_DELAY / FT_PERIOD


CONTROLLER_ARGS="-t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1"
CONTROLLER_ARGS_FT="--aft --ft_period ${FT_PERIOD} -t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1"
WORKER_ARGS="4 --othread ${THREAD_NUM}"
APPLICATION_ARGS="-i ${ITER_NUM} --sn ${SAMPLE_NUM_M}"
APPLICATION_LIB="applications/ml/lr/liblr.so"


echo -e "${Cya}Ruunnig the base experiment without faults:${RCol}"
while true; do
  start_experiment "${CONTROLLER_ARGS}" "${WORKER_ARGS}" "${APPLICATION_LIB}" "${APPLICATION_ARGS}"
  wait_to_succeed basic_completion_check ${TIME_OUT_T}
 
  NORM_ELAPSED=${ELAPSED}
  if [ "${ELAPSED}" -gt "4" ]; then
    NORM_ELAPSED=$((${ELAPSED}-4))
  fi
  if [ ${ELAPSED} -lt "${BASE_LENGTH}" ]; then
    echo -e "${Yel}[ WARNING ] base experiment with ${ITER_NUM} iterations was only ${ELAPSED} (<${BASE_LENGTH}) seconds!${RCol}"
    ITER_NUM=$((((${BASE_LENGTH}-${ELAPSED})*(${ITER_NUM}/${NORM_ELAPSED}))+${ITER_NUM}))
    echo -e "${Yel}[ WARNING ] retrying the experiment with ${ITER_NUM} iterations to get the minimum length ... ${RCol}"
    APPLICATION_ARGS="-i ${ITER_NUM} --sn ${SAMPLE_NUM_M}"
  else
    echo -e "${Gre}[ SUCCESS ] base experiment finished in ${ELAPSED} seconds.${RCol}"
    break
  fi
done
get_hash "FINAL WEIGHT HASH"
correct_hash=${HASH}


echo -e "${Cya}Ruunnig the experiment with faults:${RCol}"
start_experiment "${CONTROLLER_ARGS_FT}" "${WORKER_ARGS}" "${APPLICATION_LIB}" "${APPLICATION_ARGS}"

echo -e "${Cya}Waiting ${FAULT_DELAY} seconds before killing the first worker...${RCol}"
count_down ${FAULT_DELAY}
pause_controller
check_checkpoints ${CHKP_NUM}
kill_one_worker
resume_controller

echo -e "${Cya}Waiting ${FAULT_DELAY} seconds before killing the second worker...${RCol}"
count_down ${FAULT_DELAY}
pause_controller
check_checkpoints $((2*${CHKP_NUM}))
kill_one_worker
resume_controller

echo -e "${Cya}Waiting ${FAULT_DELAY} seconds before killing the third worker...${RCol}"
count_down ${FAULT_DELAY}
pause_controller
check_checkpoints $((3*${CHKP_NUM}))
kill_one_worker
resume_controller

wait_to_succeed basic_completion_check ${TIME_OUT_T}
echo -e "${Gre}[ SUCCESS ] experiment with fault finished in $((${ELAPSED}+3*${FAULT_DELAY})) seconds.${RCol}"
check_rewinding 3
check_hash "FINAL WEIGHT HASH" "${correct_hash}"


clean_logs
echo -e "${Gre}\n[ PASSED  ] logistic regression fault tolerance test passed successfuly!${RCol}"
exit 0


