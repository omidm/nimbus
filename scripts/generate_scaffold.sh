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

if [ -z "${NIMBUS_ROOT}" ]; then
  export NIMBUS_ROOT="$(cd "`dirname "$0"`"/..; pwd)"
fi

SRC=${NIMBUS_ROOT}/applications/scaffold/

ARGS="$@"

if [[ ${ARGS} = *--help* ]] || [[ ${ARGS} = *-h* ]]; then
  print_usage
  exit 0
fi

# get the arguments
ADIR=""
NAME=""
while (( "$#" )); do
  if [ "--path" == "$1" ] || [ "-p" == "$1" ]; then
    shift
    ADIR="$1"
  elif [ "--name" == "$1" ] || [ "-n" == "$1" ]; then
    shift
    NAME="$1"
  else
    shift
  fi
done

if [ "${ADIR}" == "" ]; then
  echo -e "${Red}ERROR: provide a path for the application!${RCol}"
  echo -e "${Red}Use -h option for the usage.${RCol}"
  exit 1
fi

if [ "${NAME}" == "" ]; then
  echo -e "${Red}ERROR: provide a name for the application!${RCol}"
  echo -e "${Red}Use \"-h\" option for the usage.${RCol}"
  exit 1
fi

# create the directory
if [ -d ${ADIR} ] || [ -f ${ADIR} ]; then
  echo -e "${Red}ERROR: the folder/file already exists at: ${ADIR}${RCol}"
  exit 1
fi

mkdir ${ADIR} &> /dev/null
if [ ! $? == 0 ]; then
  echo -e "${Red}ERROR: could not create a directory at: \"${ADIR}\"${RCol}"
  exit 1
fi

# create sub directories
dirs="protobuf_source protobuf_compiled" 
for d in ${dirs}; do
  mkdir ${ADIR}/${d} &> /dev/null
  if [ ! $? == 0 ]; then
    echo -e "${Red}ERROR: could not create a directory at: \"${ADIR}/${d}\"${RCol}"
    exit 1
  fi
done

# copy the files
files="app.cc app.h job.cc job.h data.cc data.h utils.cc utils.h Makefile"
for f in ${files}; do
  cp ${SRC}/${f} ${ADIR}
  if [ ! $? == 0 ]; then
    echo -e "${Red}ERROR: could not copy: \"${SRC}/${f}\" to \"${ADIR}\"${RCol}"
    exit 1
  fi
done

files="parameter_msg.proto vector_msg.proto Makefile"
for f in ${files}; do
  cp ${SRC}/protobuf_source/${f} ${ADIR}/protobuf_source/
  if [ ! $? == 0 ]; then
    echo -e "${Red}ERROR: could not copy: \"${SRC}/${f}\" to \"${ADIR}\"${RCol}"
    exit 1
  fi
done

files=".gitignore"
for f in ${files}; do
  cp ${SRC}/protobuf_compiled/${f} ${ADIR}/protobuf_compiled/
  if [ ! $? == 0 ]; then
    echo -e "${Red}ERROR: could not copy: \"${SRC}/${f}\" to \"${ADIR}\"${RCol}"
    exit 1
  fi
done

# set the nimbus root
sed -i.bak "s+NIMBUS_ROOT = ../../+NIMBUS_ROOT = ${NIMBUS_ROOT}+g" ${ADIR}/Makefile
rm ${ADIR}/Makefile.bak

sed -i.bak "s+NIMBUS_ROOT = ../../../+NIMBUS_ROOT = ${NIMBUS_ROOT}+g" ${ADIR}/protobuf_source/Makefile
rm ${ADIR}/protobuf_source/Makefile.bak

# set the app  name
sed -i.bak "s+ScaffoldApp+${NAME}+g" ${ADIR}/*.h
sed -i.bak "s+ScaffoldApp+${NAME}+g" ${ADIR}/*.cc
rm ${ADIR}/*.bak

# set the app library name
LIBNAME=`echo "${NAME}" | awk '{print tolower($0)}'`
sed -i.bak "s+TARGET = libscaffold.so+TARGET = lib${LIBNAME}.so+g" ${ADIR}/Makefile
rm ${ADIR}/Makefile.bak




