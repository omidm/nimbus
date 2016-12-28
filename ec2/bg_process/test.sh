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

if [ "$1" = t ]; then
  killall -v bg_process
  exit 0
fi

echo -e "${Pur}Launching $1 [arg 1] back ground processes each with $2 [arg 2] threads and array size of $3 [arg 3] (K)  tied to cores $4 [arg 4] ...${RCol}"
for i in `seq 1 $1`;
do
  taskset -c $4 ./bg_process -t $2 -s $3 &
done

