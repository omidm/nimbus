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

# This script makes the nimbus_worker process running on the current machine
# straggle by a factor of 10. It launches cpu bound processes on all available
# cores and enforces the cpu share with cgroup.
#
# If the first argument is 't', it kills all background processes, and
# essentially removes the straggler. The default straggling ration is 10x;
# meaning 9 shares of cpu goes to the background process and 1 share goes to
# nimbus_worker. You can change the ratio by passing an integer (>1) as the
# first argument.  
#
# you might need to install cgroup:
#    $ sudo apt-get install cgroup-bin
# For more information refer to README-FAKE-STRAGGLER file in this directory.

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


if [ "$1" = t ]; then
  echo -e "${Gre}Killing all bg_process's ...${RCol}"
  COUNT=$(sudo killall -v bg_process 2>&1 | grep -v "no process found" | wc -l)
  echo -e "${Gre}Killed ${COUNT} background processes.${RCol}"
  exit 0
fi

CORE_NUM=$(grep -c ^processor /proc/cpuinfo)
if [ -z ${CORE_NUM} ]; then
  echo -e "${Red}Could not detect number of cores.${RCol}"
  exit 1
fi
echo -e "${Gre}Detected ${CORE_NUM} cores.${RCol}"

PID=`ps aux | grep nimbus_worker | grep -v grep | awk '{print $2}'`
if [ -z ${PID} ]; then
  echo -e "${Red}Could not detect PID of nimbus_worker process!${RCol}"
  exit 1
fi
echo -e "${Gre}PID:  ${PID}${RCol}"

SPID=`ps -T -p ${PID} | grep -v SPID | awk '{print $2}'`
if [ -z ${SPID} ]; then
  echo -e "${Red}Could not detect SPID's of nimbus_worker threads!${RCol}"
  exit 1
fi
echo -e "${Gre}SPID's:${RCol}"
echo ${SPID}

RATIO="10"
re='^[0-9]+$'
if [[ $1 =~ ${re} ]] ; then
  RATIO=$1
fi
SHARE=$((1024 / ($RATIO - 1)))
if [ -z ${SHARE} ]; then
  echo -e "${Red}Could not compute the share.${RCol}"
  exit 1
fi
echo -e "${Gre}RATIO: ${RATIO} - SHARE: ${SHARE}(/1024)${RCol}"

sudo cgcreate -g cpu:/cpulimited
sudo cgcreate -g cpu:/lesscpulimited
sudo cgset -r cpu.shares=${SHARE} cpulimited

for t in ${SPID}; do
  sudo cgclassify -g cpu:/cpulimited ${t}
done


REL_PATH=`dirname $0`
cd ${REL_PATH}/bg_process/
for i in `seq 1 ${CORE_NUM}`; do
  ls;sudo cgexec -g cpu:lesscpulimited ./bg_process -t 1 -s 1 &
done
cd - > /dev/null

