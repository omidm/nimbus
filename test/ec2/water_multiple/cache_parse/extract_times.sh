#!/usr/bin/env bash

numworkers=$1
regex="lock"

for ((i=1; i<=$numworkers; ++i ))
do
    echo "Extracting pattern $regex from ../output/${i}_cache_time.txt"
    grep $regex ../output/${i}_cache_time.txt > ${i}_block.txt
done
