#!/bin/bash

if [ $# -ne 2 ]; then
  echo "*** PhysBAM Log searcher"
  echo "    Usage: physbam_log_searcher <directory-path> <search-string>"
  exit
fi

dirpath=$1
phrase=$2
counter=0

for i in {1..512}; do
   if grep -q $phrase "$dirpath/mpi$i.log"; then
     echo Node $i has $phrase
     let "counter+=1"
   fi
done

echo Total: $counter file had $phrase 
