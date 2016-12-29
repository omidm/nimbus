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
# straggle by a factor of 10. It assumes that the machine has 8 hyper-threaded
# processors (16 concurrent threads) and the worker is using only 4 processors
# (8 threads), launched for example with "taskset" command.

# you might need to install cgroup:
#    $ sudo apt-get install cgroup-bin
# For more information refer to README-FAKE-STRAGGLER file in this directory.

# Author: Omid Mashayekhi <omidm@stanford.edu>

if [ "$1" = t ]; then
  sudo killall -v bg_process
  exit 0
fi

PID=`ps aux | grep nimbus_worker | grep -v grep | awk '{print $2}'`
echo "PID:  " $PID

SPID=`ps -T -p $PID | grep -v SPID | awk '{print $2}'`
echo "SPID: " $SPID

sudo cgcreate -g cpu:/cpulimited
sudo cgcreate -g cpu:/lesscpulimited
sudo cgset -r cpu.shares=114 cpulimited

for t in $SPID; do
  sudo cgclassify -g cpu:/cpulimited $t
done

for i in {1..8}; do
  cd bg_process/; sudo cgexec -g cpu:lesscpulimited taskset -c 0-3,8-11 ./bg_process -t 1 -s 1 &
done

