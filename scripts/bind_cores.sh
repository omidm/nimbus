#!/bin/bash

# all the threads with the name as first argument will be bound to cores in the
# second argument using taskset
#
# example: ./bind_cores.sh scheduler 0-3 
#          binds all the threads of process scheduler to core 0-3
#

/bin/ps -eLf -fu $USER | grep $1 | grep -v grep |awk '{print $4}' | while read pid;
do
  echo $pid
  taskset -cp $2 $pid
done

