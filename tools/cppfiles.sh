#!/bin/bash

if ls ${PWD}/*.{C,H,cpp,cxx} 2> /dev/null | grep -e ".C" -e ".H" -e ".cpp" -e ".cxx" &> /dev/null; then
  echo "WARNING: Incorrectly named files exist:"
  ls ${PWD}/*.{C,H,cpp,cxx} 2> /dev/null | grep -e ".C" -e ".H" -e ".cpp" -e ".cxx"
else
  echo "All file names good."
fi

