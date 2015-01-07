#!/usr/bin/env bash

log_path="/home/nimbus/hangdata_v2"
num_workers=$1

for (( i=1; i<=$num_workers; i++))
do
    infile=${log_path}/${i}_cache_time.txt
    echo ./parse_cache_times.py -i $infile -o ${i}_parse.txt
    ./parse_cache_times.py -i ${infile} -o ${i}_parse.txt
done

echo ./plot_cache_times.py 1 $num_workers
./plot_cache_times.py 1 $num_workers
