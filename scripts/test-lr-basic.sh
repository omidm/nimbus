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

# Checks controller and worker basic functionalities against stencil 1D
# application.

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


DESC="Runs comprehensive tests for reduction functionality"
DESC=${DESC}"\nagainst logistic regression application, tests include: "
DESC=${DESC}"\n   1. explicit and automatic reduction."
DESC=${DESC}"\n   2. reduction with and without combiner."
DESC=${DESC}"\n   3. multi-threaded controller and workers."
DESC=${DESC}"\n   4. different components of templates."
USAGE="./scripts/test-stencil-basic.sh"


print_help "${DESC}" "${USAGE}" "$1"

ITER_NUM=100
SAMPLE_NUM_M=".001"
BATCH_NUM=6
THREAD_NUM=16
TIME_OUT_T=20
APPLICATION_LIB="applications/ml/logistic_regression/liblr.so"

echo -e "${Cya}Ruunnig the base experiment:${RCol}"
echo -e "${Cya}CONTROLLER  : single-threaded w/o templates${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : explicit reduction using read/write set only${RCol}"
start_experiment "--dct" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --dar" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished base experiment in ${ELAPSED} seconds.${RCol}" 
get_hash "FINAL WEIGHT HASH"
correct_hash=${HASH}


echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : single-threaded w/o templates${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "--dct" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : single-threaded w/o templates${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "--dct" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/o templates${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dct" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/o templates${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dct" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dbm" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dbm" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and binding memoization${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dcb --dwt" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and binding memoization${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dcb --dwt" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dwt" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/o templates${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} --dwt" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/ worker template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM}" "1 --det" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/ worker template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM}" "1 --det" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM}" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, single-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM}" "1" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : two, single-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} -w 2 --split 2 1 1" "2" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : two, single-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} -w 2 --split 2 1 1" "2" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : four, single-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1" "4" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : four, single-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1" "4" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"



echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, multi-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM}" "1 --othread ${THREAD_NUM}" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc" 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : one, multi-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM}" "1 --othread ${THREAD_NUM}" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"


echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : four, multi-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/o combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1" "4 --othread ${THREAD_NUM}" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} --drc"
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"

echo -e "${Cya}Ruunnig experiment for:${RCol}"
echo -e "${Cya}CONTROLLER  : multi-threaded w/ controller template and cascaded binding${RCol}"
echo -e "${Cya}WORKERS     : four, multi-threaded w/ worker and execution template${RCol}"
echo -e "${Cya}APPLICATION : automatic reduction w/ combiner${RCol}"
start_experiment "-t ${THREAD_NUM} -a ${BATCH_NUM} -w 4 --split 4 1 1" "4 --othread ${THREAD_NUM}" "${APPLICATION_LIB}" "-i ${ITER_NUM} --sn ${SAMPLE_NUM_M} " 
wait_to_succeed basic_completion_check ${TIME_OUT_T} 
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 
check_hash "FINAL WEIGHT HASH" "${correct_hash}"


clean_logs

echo -e "${Gre}\n[ PASSED  ] all logistic regression basic tests passed successfuly!${RCol}"
exit 0


