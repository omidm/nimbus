#!/usr/bin/python
######################################################################
# Copyright 2005, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#######################################################################
# Sample sim wrapper script
#######################################################################
import time
import sys
import os
import re
import SIM_WRAPPER_ENGINE

#######################################################################
# Allocate wrapper script
#######################################################################
cmd=sys.argv[1:]
wrapper=SIM_WRAPPER_ENGINE.SIM_WRAPPER(cmd)
wrapper.poll_interval=5 # poll/report every 5 seconds

#######################################################################
# Clean up the variables that SIM_WRAPPER sets
#######################################################################
def Report_Start(pid,rpc,session):
    rpc.Remove_Status_If_Exists(session,"memory")
    rpc.Remove_Status_If_Exists(session,"last_frame")
    rpc.Remove_Status_If_Exists(session,"last_frame_time")
    rpc.Remove_Status_If_Exists(session,"last_frame_duration")
# Register with the wrapper class
wrapper.start_callback=Report_Start

#######################################################################
# Register Memory Reporting
#######################################################################
# this will be called every time the reporting thread wakes up
max_memory=0
def Report_Memory(pid,rpc,session):
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
    usage=Memory_Usage(pid)
    if usage>0:
        global max_memory
        max_memory=max(max_memory,usage)
        rpc.Update_Status(session,"memory",usage)
        rpc.Update_Status(session,"max_memory",max_memory)
# Register with the wrapper class
wrapper.polled_events.append(Report_Memory)

#######################################################################
# Register Frame Reporting
#######################################################################
# this is applied to each line of input read from the simulation stdout
frame_regexp=re.compile("END Frame (\d+)[ ]+(\d+\.\d+) s")
# this is called with the match structure of applying the regular expression 
def Match_Frame(frame_match):
    frame=int(frame_match.group(1))
    frametime=float(frame_match.group(2))
    return (frame,frametime) # pack and return data
def Report_Frame(pid,rpc,session,data):
    frame,frametime=data # unpack data
    # send the info to the server
    rpc.Update_Status(session,"last_frame",frame)
    rpc.Update_Status(session,"last_frame_duration",frametime)
    rpc.Update_Status(session,"last_frame_time",time.time())
# Register with the wrapper class
wrapper.simulation_events.append((frame_regexp,Match_Frame,Report_Frame))

#######################################################################
# Run the script
#######################################################################
wrapper.run()
