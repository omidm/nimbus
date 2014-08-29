#!/bin/bash

echo -n "Checking files follow naming guidelines .{c,h,cc}: "
if ls ${PWD}/*.{C,H,cpp,cxx} 2> /dev/null | grep -e ".C" -e ".H" -e ".cpp" -e ".cxx" &> /dev/null; then
  echo "\n  WARNING: Incorrectly named files exist:"
  ls ${PWD}/*.{C,H,cpp,cxx} 2> /dev/null | grep -e ".C" -e ".H" -e ".cpp" -e ".cxx"
else
  echo "OK"
fi

