#!/usr/bin/python
import os
import sys
import time

os.umask(2) # make everyone able to read stuff

############################################
# Add new test files right here
###########################################
tests=["fluid.py","solids.py","arb.py","rigid_bodies.py","deformable_rigid.py"]

# set display to local host
os.environ["DISPLAY"]="localhost:0"

output_directory,data_directory=None,None
try:
    data_directory,output_directory_base=sys.argv[1],sys.argv[2]
except:
    print "Usage: %s <physbam pub data> <output directory base>"%(sys.argv[0])
    sys.exit(1)

# setup data directory
if not os.path.exists(data_directory) or not os.path.isdir(data_directory):
    print "%s is not a directory"%(data_directory)
    sys.exit(1)

# symlink public_data if it isn't there
try:
   os.unlink("Public_Data")
except:
   pass
os.symlink(data_directory,"Public_Data")

# setup output directory
times=time.gmtime()
time_string="%04d%02d%02d-%02d%02d%02d"%(times[0],times[1],times[2],times[3],times[4],times[5])
time_output_directory=os.path.join(output_directory_base,time_string) # Time string YYYYMMDD-HHMM
if os.path.exists(time_output_directory):
    print "Already have test at this time!"
    sys.exit(1)
os.mkdir(time_output_directory)

# process each test
for test in tests:
    test_dir=os.path.join(time_output_directory,test.replace(".py",""))
    os.mkdir(test_dir)
    exitcode=os.system("python Scripts/test/%s %s"%(test,test_dir))


