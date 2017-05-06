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

# This script generates an initial scaffold for a Nimbus application.
# Usage:
#           ./generate-scaffold.sh [options]
# Options:
#           -h [ --help ] produce help message
#           -p [ --path ] [REQUIRED] the path to a folder that will contain the
#                         application files. The path could be absolute or
#                         relative compared to where this script is called
#                         from. The folder could not exist.
#           -n [ --name ] [REQUIRED] the name of the application.
#
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
  echo -e "${Cya}This script generates an initial scaffold for a Nimbus application."
  echo -e "Usage:"
  echo -e "          ./generate-scaffold.sh [options]"
  echo -e "Options:"
  echo -e "          -h [ --help ] produce help message"
  echo -e "          -p [ --path ] [REQUIRED] the path to a folder that will contain"
  echo -e "                        the application files. The path could be absolute"
  echo -e "                        or relative compared to where this script is called"
  echo -e "                        from. The folder could not already exist."
  echo -e "          -n [ --name ] [REQUIRED] the name of the application."
  echo -e "${RCol}"
}

if [ -z "${NIMBUS_HOME}" ]; then
  export NIMBUS_HOME="$(cd "`dirname "$0"`"/..; pwd)"
fi

ARGS="$@"

if [[ ${ARGS} = *--help* ]] || [[ ${ARGS} = *-h* ]]; then
  print_usage
  exit 0
fi


PATH=""
NAME=""
while (( "$#" )); do
  if [ "--path" == "$1" ] || [ "-p" == "$1" ]; then
    shift
    PATH="$1"
  elif [ "--name" == "$1" ] || [ "-n" == "$1" ]; then
    shift
    NAME="$1"
  else
    shift
  fi
done

if [ "${PATH}" == "" ]; then
  echo -e "${Red}ERROR: provide a path for the application!${RCol}"
  echo -e "${Red}Use -h option for the usage.${RCol}"
  exit 1
fi

if [ "${NAME}" == "" ]; then
  echo -e "${Red}ERROR: provide a name for the application!${RCol}"
  echo -e "${Red}Use \"-h\" option for the usage.${RCol}"
  exit 1
fi



