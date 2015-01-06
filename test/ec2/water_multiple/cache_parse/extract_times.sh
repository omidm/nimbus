#!/usr/bin/env bash

numworkers=$1
regex="lock"
logpath="../output"

for (( i=1; i<=$numworkers; ++i ))
do
    echo "Extracting pattern $regex from ${logpath}/${i}_cache_time.txt"
    grep $regex ${logpath}/${i}_cache_time.txt > ${i}_block.txt
done
