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

function print_usage {
  echo -e "${Cya}Runs water simulation againts nimbus with two workers."
  echo -e "\nUsage:"
  echo -e "./scripts/test-water.sh"
  echo -e "${RCol}"
}


TIME_OUT_T=100
CONTROLLER_ARGS="-w 2 --split 1 2 1"
WORKER_ARGS="2 --othread 2 "
WATER_ARGS="--wl 0.35 -e 10"

if [ -z "${NIMBUS_HOME}" ]; then
  export NIMBUS_HOME="$(cd "`dirname "$0"`"/..; pwd)"
fi

if [ -z "${DBG}" ]; then
  export DBG="error"
fi

if [ -z "${TTIMER}" ]; then
  export TTIMER="l1"
fi

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

${NIMBUS_HOME}/scripts/stop-workers.sh &> /dev/null
${NIMBUS_HOME}/scripts/stop-controller.sh &> /dev/null
${NIMBUS_HOME}/scripts/start-controller.sh ${CONTROLLER_ARGS} &> /dev/null
${NIMBUS_HOME}/scripts/start-workers.sh ${WORKER_ARGS} -l applications/physbam/water/libwater_app.so ${WATEER_ARGS} &> /dev/null

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

if [ ${success} == "1" ]; then
  echo -e "${Gre}[ SUCCESS ] finished in $((${end_time}-${start_time})) seconds.${RCol}"
else
  echo -e "${Red}[ TIMEOUT ] experiment did not finish before time out!${RCol}"
  exit 1
fi


echo -e "${Cya}Checking the simulation output:${RCol}"

FOLDER_COUNT=$(ls ${NIMBUS_HOME}/nodes/nimbus_worker/output/ 2>/dev/null | wc -w)

if [ "${FOLDER_COUNT}" == "12" ]; then
  echo -e "${Gre}[ SUCCESS ] output file seems complete!${RCol}"
else
  echo -e "${Red}[ FAILED  ] output file does not seem complete!${RCol}"
  exit 1
fi


echo -e "${Cya}Displaying the results:${RCol}"

${NIMBUS_HOME}/applications/physbam/physbam-app/opengl_3d/opengl_3d ${NIMBUS_HOME}/nodes/nimbus_worker/output/ &> /dev/null &

WINDOWID=$(xdotool search --name opengl_3d*)

TIME_OUT_T=10
start_time=$(date +%s)
end_time=$(date +%s)
progress_bar="waiting to open the opengl_3d window ..."
while [ "$((${end_time}-${start_time}))" -lt "${TIME_OUT_T}" ]; do
  WINDOWID=$(xdotool search --name opengl_3d*)
  if ! [ -z ${WINDOWID} ]; then
    end_time=$(date +%s)
    break
  else
    end_time=$(date +%s)
    echo -ne "${Yel}${progress_bar} \r${RCol}"
    progress_bar=${progress_bar}"."
    sleep 1
  fi
done
if ! [ -z ${WINDOWID} ]; then
  echo -e "${Gre}[ SUCCESS ] a window dispalying simulation should have been opened!${RCol}"
else
  echo -e "${Red}[ TIMEOUT ] window did not open before time out!${RCol}"
  exit 1
fi

xdotool key --window ${WINDOWID} "p"

echo -e "${Gre}[ VISUAL CHECK ] hit p (play/stop play), s (step forward), ctrl+s (step backward) ${RCol}"
echo -e "${Gre}[ VISUAL CHECK ] by default it starts playing simulation and there are 10 frames!${RCol}"


echo -e "${Gre}[ VISUAL CHECK ] window closes and the logs are cleaned in 10 seconds!${RCol}"
sleep 10
xdotool windowactivate --sync ${WINDOWID} key --clearmodifiers --delay 100 alt+F4
make ${NIMBUS_HOME}/ clean-logs &> /dev/null

# if running the test over ssh -X use this command to close!
# wmctrl -c "opengl"

exit 0

