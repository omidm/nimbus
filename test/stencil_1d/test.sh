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

if [ "$1" = t ]
then
killall -v worker
exit
fi

# echo -e "${Pur}launching scheduler version I ...${RCol}"
# ../scheduler_v1/scheduler 2  > scheduler.txt &
# sleep 1
 
echo -e "${Gre}export DBG=errors...${RCol}"
export DBG=error

echo -e "${Pur}Launching the workers ...${RCol}"
for i in `seq 1 $1`;
do
  ./worker -sip localhost -sport 5900 -port 590$i &
done

