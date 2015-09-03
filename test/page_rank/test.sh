#!/usr/bin/env bash

rm scheduler_log
../scheduler_v2/scheduler --port 5000 -w 1 &>> scheduler_log & 
./worker --port 5001 --cport 5000 --cip 127.0.0.1 --input graph-1p/ --iterations 0
