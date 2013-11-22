#!/bin/bash

export LD_LIBRARY_PATH=../../lib/:../../application/test_quhang/projection/:$LD_LIBRARY_PATH
mpirun -np 1 ./projection_multiple 1 : -np 1 ./projection_multiple 2
