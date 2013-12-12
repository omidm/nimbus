#!/bin/bash

#export LD_LIBRARY_PATH=../../lib/:../../application/projection-two-worker/:$LD_LIBRARY_PATH
mpirun -np 1 ./projection_multiple 1 : -np 1 ./projection_multiple 2
