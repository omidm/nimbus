#!/bin/bash

DIR=`dirname $0`
FILE=/tmp/headers-`perl -e 'print rand();'`

( cd $PUBLIC ; find -name '*.h' ) | sed 's@\./@@' > $FILE;
find $PUBLIC -name '*.h' -or -name '*.cpp' | xargs perl -i $DIR/fix_headers.pl $FILE 2>&1
find $PROJ -name '*.h' -or -name '*.cpp' | xargs perl -i $DIR/fix_headers.pl $FILE 2>&1
find $PHYSBAM/Tools -name '*.h' -or -name '*.cpp' | xargs perl -i $DIR/fix_headers.pl $FILE 2>&1
find $PHYSBAM/Tests -name '*.h' -or -name '*.cpp' | xargs perl -i $DIR/fix_headers.pl $FILE 2>&1

rm $FILE
