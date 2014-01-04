#!/bin/bash

echo "testing the scheduler version 1 against stencil application version multi and two workers ..."

echo "launching scheduker version I ...";
../scheduler_v1/scheduler 2  > scheduler.txt &
sleep 2

echo "launching stencil worker multi number 1 ...";
../stencil_worker/worker_multi 1 > worker1.txt &
sleep 2

echo "launching stencil worker multi number 2 ...";
../stencil_worker/worker_multi 2 > worker2.txt

echo "wait for termination ... ";

#read the files for correct output.

rm -f *.txt

