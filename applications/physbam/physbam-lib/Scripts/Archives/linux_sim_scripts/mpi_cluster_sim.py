#!/usr/bin/python

import os
import sys

#submit_path=os.path.dirname(sys.argv[0])
#if submit_path==".": submit_path=os.getcwd()
#print submit_path

if len(sys.argv)<4:
    print "Usage: %s <# nodes> <name of sim> <binary of sim> [arguments to sim]"%sys.argv[0]
    print " "
    print "Note: the sim will be run with the current directory being in the directory"
    print "where the binary lives"
    print " "
    print "You can use # nodes 1 to run single threaded"
    sys.exit(1)
    
nodes_str=sys.argv[1]
name=sys.argv[2]
binary=sys.argv[3]
command=sys.argv[4:]

#############################
# find parameters and set current directory
#############################

run_directory=os.path.dirname(binary)
print run_directory
if run_directory=="." or run_directory=="": run_directory=os.getcwd()
node_count=int(nodes_str)
sim_command=binary+" "+" ".join(command)

os.chdir(run_directory)

#############################
# make a mon session
#############################
from pd.common import CONNECT
mon=None
mon_sessiond=None
try:
    mon=CONNECT.mon_client()
    mon_session=mon.Create_Session(os.environ["USER"])['id']
    mon.Label_Session(mon_session,name)
except:
    mon.close()
else:
    mon.close()
if not mon_session:
    print "Failed to register mon session, server possibly not running"
    sys.exit(1)

#############################
# write the script file
#############################
script="""#!/bin/tcsh
cd %s
export PMON_SESSION_ID=%d
"""%(run_directory,mon_session)

# if MPI
mpi_string="" # default to no parallel environment
if node_count>1:
    mpi_string="-pe mpi %d"%node_count
    script+="""
cat $PE_HOSTFILE  | awk '{print $1 " cpu=" $2}' > hostfile.${JOB_ID}
cat hostfile.${JOB_ID}

lamboot hostfile.${JOB_ID}
echo Running on nodes
lamnodes 
/usr/local/adm/pd/pd/mon/MON_WRAPPER.py mpirun -np %d %s

lamhalt
"""%(node_count,sim_command)
else:
    script+="""
echo running on `hostname`
/usr/local/adm/pd/pd/mon/MON_WRAPPER.py %s
"""%(sim_command)

print "-------------Script------------------------------------------------"
print script
print "--------------------------------------------------------------------"

#############################
# do the submit
#############################

queue_name="cluster"
if node_count>1: queue_name="mpi"
stdin,stdout=os.popen4("qsub -shell y -S /bin/bash -cwd %s -q %s -N %s"%(mpi_string,queue_name,name))
stdin.write(script)
stdin.close()
print stdout.read()


