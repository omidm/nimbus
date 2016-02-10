import sys
import os
import string
import time
from stat import *

if len(sys.argv) < 5:
    print "<generated_scene_folder> <output_folder> <start frame> <end frame> <ray-tracing binary>"
    sys.exit()

for i in range(string.atoi(sys.argv[3]), string.atoi(sys.argv[4])):
    checkfile =  sys.argv[2] + "render.%(#)05d.png" % {"#" : i}

    resubmitflag = 0

    try:
        st = os.stat(checkfile)
    except OSError:
        print checkfile + " is missing, resubmitting..."
        resubmitflag = 1
    else:
        if st[ST_SIZE] == 0:
            print checkfile + " is missing, resubmitting..."
            resubmitflag = 1

    if resubmitflag == 1:
        os.system("qsub " + sys.argv[1] + "tmp_submit_script" + str(i) + ".sh")
