#!/bin/bash

if [ -z $3 ]; then echo "Usage: $0 format start end"; exit; fi

start_frame=$2
end_frame=$3

for ((i=$start_frame; i<=$end_frame; i++)); do
	if [ ! -e `printf $1 $i` ]; then
		echo $i;
	fi
done
