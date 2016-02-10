#!/bin/bash

filename="Makefile"

while getopts "f:" option
do
    case $option in
	f ) filename=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))

if [ -z "$1" ]; then
    echo "Usage: create_lib_makefile.sh [-f <makefile>] <library name> <cpp files> ..."
    exit;
fi

# Preserve contents above CUT OFF POINT from previous Makefile
if [ -e $filename ]; then
    grep "^# CUT OFF POINT" $filename > /dev/null
    if [ $? == 0 ]; then
	mv $filename /tmp/old_makefile
	cat /tmp/old_makefile | awk 'BEGIN { echo=1 } /^\# CUT OFF POINT/ {echo=0} {if (echo==1) print $0}' > $filename
    else
	rm -f $filename
	touch $filename
    fi
else
    touch $filename
fi

cat >> $filename << EOF
# CUT OFF POINT (lines above will be preserved, lines below will be overwritten)
# Set any necessary USE_* variables (e.g. USE_OPENGL=yes) above CUT OFF POINT.

TARGETS = $1
TARGET_TYPE = STATIC_LIBRARY

LOCAL_SRC = \\
EOF

shift
while [ -n "$1" ]; do
    echo "	$1 \\" >> $filename
    shift
done

cat >> $filename << EOF

include \$(PHYSBAM)/Public_Library/Makefile.common
EOF
