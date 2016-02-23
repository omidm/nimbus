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

# Test written for water simulation against Nimbus.

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


DESC="Runs water simulation againts nimbus with two workers."
USAGE="./scripts/test-water-basic"

print_help "${DESC}" "${USAGE}" "$1"


TIME_OUT_T=100
CONTROLLER_ARGS="-w 2 --split 1 2 1"
WORKER_ARGS="2 --othread 2 "
APPLICATION_ARGS="--wl 0.35 -e 10"
APPLICATION_LIB="applications/physbam/water/libwater_app.so"


echo -e "${Yel}[ WARNING ] this test requires building physbam librray and applications first.${RCol}"
echo -e "${Yel}            if not built yet, you can issue \'make physbam\' in the nimbus root.${RCol}"
echo -e "${Yel}            it is not automated since physbam Makefile is flaky and it may build${RCol}"
echo -e "${Yel}            everything from scratch again, even after successful build.${RCol}"


echo -e "${Cya}Checking that old output files are removed:${RCol}"

make ${NIMBUS_HOME}/ clean-logs &> /dev/null

FOLDER=$(ls ${NIMBUS_HOME}/nodes/nimbus_worker/output/ 2>/dev/null)

if [ -z "${FOLDER}"]; then
  echo -e "${Gre}[ SUCCESS ] no old output file found${RCol}"
else
  echo -e "${Red}[ FAILED  ] seems that there are still old output files!${RCol}"
  exit 1
fi


echo -e "${Cya}Ruunnig water simulation with two workers against nimbus controller:${RCol}"

start_experiment "${CONTROLLER_ARGS}" "${WORKER_ARGS}" "${APPLICATION_LIB}" "${APPLICATION_ARGS}"
wait_to_succeed basic_completion_check ${TIME_OUT_T}
echo -e "${Gre}[ SUCCESS ] finished experiment in ${ELAPSED} seconds.${RCol}" 

echo -e "${Cya}Checking the simulation output:${RCol}"

FOLDER_COUNT=$(ls ${NIMBUS_HOME}/nodes/nimbus_worker/output/ 2>/dev/null | wc -w)

if [ "${FOLDER_COUNT}" == "12" ]; then
  echo -e "${Gre}[ SUCCESS ] output file seems complete!${RCol}"
else
  echo -e "${Red}[ FAILED  ] output file does not seem complete!${RCol}"
  exit 1
fi


echo -e "${Cya}Displaying the results:${RCol}"
if ! [ -z ${SSH_TTY} ]; then
  echo -e "${Yel}[ WARNING ] this is a ssh session, some visual functionality may not work.${RCol}"
fi

${NIMBUS_HOME}/applications/physbam/physbam-app/opengl_3d/opengl_3d ${NIMBUS_HOME}/nodes/nimbus_worker/output/ &> /dev/null &
function window_check {
  local WINDOWID=$(xdotool search --name opengl_3d*)
  if ! [ -z ${WINDOWID} ]; then
    echo "true"
  else
    echo "false"
  fi
}
wait_to_succeed window_check 10

WINDOWID=$(xdotool search --name opengl_3d*)
xdotool key --window ${WINDOWID} "p"

echo -e "${Gre}[ VISUAL CHECK ] hit p (play/stop play), s (step forward), ctrl+s (step backward) ${RCol}"
echo -e "${Gre}[ VISUAL CHECK ] by default it starts playing simulation and there are 10 frames!${RCol}"

echo -e "${Gre}[ VISUAL CHECK ] window closes and the logs are cleaned in 10 seconds!${RCol}"
wait_with_bar 10
if [ -z ${SSH_TTY} ]; then
  xdotool windowactivate --sync ${WINDOWID} key --clearmodifiers --delay 100 alt+F4
else
  wmctrl -c "opengl"
fi

clean_logs
exit 0

