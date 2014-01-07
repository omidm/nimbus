#!/bin/bash

# **************************
# Text Reset
RCol='\x1B[0m'    
# Regular           
Bla='\x1B[0;30m';
Red='\x1B[0;31m';
Gre='\x1B[0;32m';
Yel='\x1B[0;33m';
Blu='\x1B[0;34m';
Pur='\x1B[0;35m';
Cya='\x1B[0;36m';
Whi='\x1B[0;37m';
# **************************

echo -e "${Pur}Testing the scheduler version 1 against stencil application version multi and two workers. ${RCol}"

NIMBUS_ROOT="../../"
NIMBUS_LIB="libnimbus.so"
SCHEDULER_PATH="../../test/scheduler_v1/"
SCHEDULER="scheduler"
APPLICATION_PATH="../../application/stencil_1d_multi/"
APPLICATION="libstencil_1d_multi.so"
WORKER_PATH="../../test/stencil_worker/"
WORKER="worker_multi"
TEMP_FILE="$(pwd)/_nimbus_temp_file_.txt"

# check whether to make clean or not, by default it will make clean.
CLEAN="make clean"
if [ "$1" = n ]
then
  CLEAN=""
fi

echo -e "${Pur}Building Nimbus library ...${RCol}"
# build the nimbus library
cd ${NIMBUS_ROOT}; rm -f ${NIMBUS_LIB}; $CLEAN &> ${TEMP_FILE}; make -j 12 &> ${TEMP_FILE};
if [ -f ${NIMBUS_LIB} ]
then
  echo -e "${Gre}SUCCESS: built the Nimbus library successfully. ${RCol}"
  cd - &> ${TEMP_FILE};
else
  echo -e "${Red}ERROR: could not build Nimbus library. ${RCol}"
  exit;
fi

echo -e "${Pur}Building the application ...${RCol}"
# build the application
cd ${APPLICATION_PATH}; rm -f ${APPLICATION}; $CLEAN &> ${TEMP_FILE}; make &> ${TEMP_FILE};
if [ -f ${APPLICATION} ]
then
  echo -e "${Gre}SUCCESS: built the application successfully. ${RCol}"
  cd - &> ${TEMP_FILE};
else
  echo -e "${Red}ERROR: could not build application. ${RCol}"
  exit;
fi

echo -e "${Pur}Building the scheduler ...${RCol}"
# build the scheduler
cd ${SCHEDULER_PATH}; rm -f ${SCHEDULER}; $CLEAN &> ${TEMP_FILE}; make &> ${TEMP_FILE};
if [ -f ${SCHEDULER} ]
then
  echo -e "${Gre}SUCCESS: built the scheduler successfully. ${RCol}"
  cd - &> ${TEMP_FILE};
else
  echo -e "${Red}ERROR: could not build the scheduler. ${RCol}"
  exit;
fi

echo -e "${Pur}Building the worker ...${RCol}"
# build the worker
cd ${WORKER_PATH}; rm -f ${WORKER}; $CLEAN &> ${TEMP_FILE}; make &> ${TEMP_FILE};
if [ -f ${WORKER} ]
then
  echo -e "${Gre}SUCCESS: built the worker successfully. ${RCol}"
  cd - &> ${TEMP_FILE};
else
  echo -e "${Red}ERROR: could not build the worker. ${RCol}"
  exit;
fi


echo -e "${Pur}Startig the test ...${RCol}"

echo -e "${Pur}launching scheduler version I ...${RCol}"
../scheduler_v1/scheduler 2  > scheduler.txt &
sleep 1
 
echo -e "${Pur}launching stencil worker multi number 1 ...${RCol}"
../stencil_worker/worker_multi 1 > worker1.txt &
sleep 1
 
echo -e "${Pur}launching stencil worker multi number 2 ...${RCol}"
../stencil_worker/worker_multi 2 > worker2.txt
 
echo "last worker terminated, checking the results ... ";

output=$(./error_check.py scheduler.txt)
if [ ${output} -ne "0" ]
then
  echo -e "${Red}ERROR: found ${output} ERROR and WARNING in scheduler. ${RCol}"
else
  echo -e "${Gre}SUCCESS: No Error or WARNING in scheduler. ${RCol}"
fi

output=$(./error_check.py worker1.txt)
if [ ${output} -ne "0" ]
then
  echo -e "${Red}ERROR: found ${output} ERROR and WARNING in worker 1.${RCol}"
else
  echo -e "${Gre}SUCCESS: No Error or WARNING in worker 1. ${RCol}"
fi

output=$(./error_check.py worker2.txt)
if [ ${output} -ne "0" ]
then
  echo -e "${Red}ERROR: found ${output} ERROR and WARNING in worker 2. ${RCol}"
else
  echo -e "${Gre}SUCCESS: No Error or WARNING in worker 2. ${RCol}"
fi



output=$(./output_check.py worker1.txt "-111985389, 217213596, -305294628, 362485704, ")
if [ ${output} != "True" ]
then
  echo -e "${Red}ERROR: output does not match for worker 1. ${RCol}"
  exit;
else
  echo -e "${Gre}SUCCESS: output matched for worker 1. ${RCol}"
fi

output=$(./output_check.py worker2.txt "-376143576, 340031724, -256983732, 138097656, ")
if [ ${output} != "True" ]
then
  echo -e "${Red}ERROR: output does not match for worker 2. ${RCol}"
  exit;
else
  echo -e "${Gre}SUCCESS: output matched for worker 2. ${RCol}"
fi

echo -e "${Gre}TEST PASSED! ${RCol}"


rm -f *.txt ${TEMP_FILE}

