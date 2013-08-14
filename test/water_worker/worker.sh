#!/bin/bash

export LD_LIBRARY_PATH=../../lib/:../../application/water_test/:$LD_LIBRARY_PATH
./worker
