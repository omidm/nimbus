#!/usr/bin/python
######################################################################
# Copyright 2005, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#######################################################################
# Standard test wrapper script
# run with run_standard_tests.py <octree,3d> <example#> <resolution#>
#######################################################################
import os
import sys
import re
import time
import popen2

# Define tests
test_names={1:"Flat surface",2:"Falling drop",3:"4 sources",4:"Sphere splashing",5:"Glass filling"}
binaries={"octree":"./fluids_octree","3d":"./fluids_3d"}

# returns memory usage given a pid
def Memory_Usage(pid):
    scale = {'kB':1/1024.0,'mB': 1,'KB':1/1024.0,'MB':1}
    try:
        fp=open("/proc/%d/status"%pid)
        while 1:
            i=fp.readline()
            if i=="":break
            elif i.startswith("VmSize:"):
                label,size,unit=i.split(None,3)
                return float(size)*scale[unit]
    except: return -1

# test output directory
def Make_Example_Path(example,resolution):
    year,month,day,hour,minute,second,weekday,yearday,dst=time.localtime()
    data_directory=os.path.join("data","%s_Ex_%02d_Res_%02d_%s_%04d%02d%02d_%02d%02d"%(os.environ['HOSTNAME'],example,resolution,os.environ['HOSTNAME'],year,month,day,hour,minute))
    try: os.makedirs(data_directory)
    except: pass

    info=["Host: %s"%os.environ["HOSTNAME"],"Example: %d (%s)"%(example,test_names[example]),"Resolution: %d"%resolution,"Run Date: %s"%time.ctime()]
    open(os.path.join(data_directory,"info.py"),"w").write(repr(info))

    return data_directory

def Run_Example(type,example,resolution):
    path=Make_Example_Path(example,resolution)
    platform_modifier=""
    if os.environ.has_key('PLATFORM'): platform_modifier="_"+os.environ['PLATFORM']
    cmd="%s%s -resolution %d -o %s %d"%(binaries[type],platform_modifier,resolution,path,example)
    print "------------------------------------------------------------------------------------"
    print "Test %d Resolution %d (%s)"%(example,resolution,test_names[example])
    print "------------------------------------------------------------------------------------"
    print "Running '%s'"%cmd

    # Copy scene
    try: shutil.copy(os.path.join("Standard_Tests","%d.scene"%example),os.path.join(path,"render.scene"))
    except: pass

    # start command with pipe and setup output directory
    pipe=popen2.Popen4(cmd)
    input,pid=pipe.fromchild,pipe.pid
    out=open(os.path.join(path,"output.txt"),"w")
    memory_samples={}
    # detect frames
    frame_regexp=re.compile("END Frame (\d+).+(\d+\.\d+) s")

    while 1:
        line=input.readline()
        if line=="": break
        out.write(line)
        line=line.strip()
        frame_match=frame_regexp.match(line)
        if frame_match:
            try:
                frame=int(frame_match.group(1))
                frametime=float(frame_match.group(2))
                print "%5d\t\t%f s"%(frame,frametime)
                memory_samples[frame]=Memory_Usage(pid)
                open(os.path.join(path,"memory.py"),"w").write(repr(memory_samples))
                sys.stdout.flush()
            except: pass
        open(os.path.join(path,"memory.py"),"w").write(repr(memory_samples))

# Get arguments
try:
    type,example,resolution=sys.argv[1],int(sys.argv[2]),int(sys.argv[3])
    if not binaries.has_key(type): raise Exception
    if not test_names.has_key(example): raise Exception
except:
    print "Usage: run_standard_tests.py <type '3d' or 'octree'> <example> <resolution>"
    sys.exit(0)


Run_Example(type,example,resolution)
