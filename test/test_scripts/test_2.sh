#!/bin/bash

echo "testing the scheduler version 1 against stencil application version single and one worker ..."

echo "launching scheduker version I ...";
../scheduler_v1/scheduler 1  > scheduler.txt &
sleep 2

echo "launching stencil worker single number 1 ...";
../stencil_worker/worker_single 1 > worker1.txt

echo "wait for termination ... ";

#read the files for correct output.

rm -f *.txt


