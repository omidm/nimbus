#!/bin/bash

DIR=`dirname $0`
FILE=/tmp/headers-`perl -e 'print rand();'`

( cd $PUBLIC ; find -name '*.h' ) | sed 's@\./@@' > $FILE;
echo $* | xargs perl -i $DIR/fix_headers.pl $FILE 2>&1

rm $FILE
