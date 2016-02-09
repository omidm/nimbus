#!/usr/bin/env bash

numworkers=$1
regex="lock"

for ((i=1; i<=$numworkers; ++i ))
do
    echo "Histograms for worker ${i}"
    ./parse_to_hist.py -i ${i}_block.txt -o ${i}_cdf.png -p cdf
    ./parse_to_hist.py -i ${i}_block.txt -o ${i}_pdf.png -p pdf
done
