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

echo -e "${Gre}export TTIMER=none...${RCol}"
export TTIMER=l1

echo -e "${Pur}Launching $1 workers  each with $2 threads...${RCol}"
echo -e "${Pur}Extra arguments: $3 ${RCol}"
for i in `seq 1 $1`;
do
  # ./worker $3 --cip localhost --cport 5900 --othread $2 -p 590$i --pnx 8 --pny 8 --pnz 8 --ppnx 8 --ppny 8 --ppnz 8 -s 160 --psl 0 -e 4 --wl 0.35 &
  # ./worker $3 --cip localhost --cport 5900 --othread $2 -p 590$i --pnx 4 --pny 4 --pnz 4 --ppnx 4 --ppny 4 --ppnz 4 -s 80 --psl 0 -e 4 --wl 0.35 &
  ./worker $3 --cip localhost --cport 5900 --othread $2 -p 590$i --psl 0 --wl 0.35 -e 10 & 
done

