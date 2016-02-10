#!/bin/bash
function test {
    "$@"
    status=$?
    if [ $status -ne 0 ]; then
        echo "error with $1"
    fi
    return $status
}
cd particle_simple_rw
cd Build/Debug
rm particle_simple_rw-debug
make
cd ../..
test ./particle_simple_rw-debug 
