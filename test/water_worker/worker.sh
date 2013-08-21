#!/bin/bash

export LD_LIBRARY_PATH=../../lib/:../../application/water_test_single/:$LD_LIBRARY_PATH
./worker
