#!/usr/bin/python
######################################################################
# Copyright 2005, Frank Losasso, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#######################################################################
# Classes that help create a wrapper script
#######################################################################
import os
import sys
import re
import time
import threading
import popen2
import traceback
import signal
from pd.common import CONNECT

#######################################################################
# Class DUMMY_RPC for testing
#######################################################################
class DUMMY_RPC:
    def Update_Status(self,session,key,value):
        print "Updating on session=%d  %s=%s"%(session,key,repr(value))
    def Remove_Status(self,session,key):
        print "Updating on session=%d  clear %s"%(session,key)

#######################################################################
# Class PRODUCER_CONSUMER_QUEUE
#######################################################################
class PRODUCER_CONSUMER_QUEUE:
    def __init__(self):
        self.events=[]
        self.mutex=threading.Lock()

    def Produce(self,data):
        self.mutex.acquire()
        self.events.append(data)
        self.mutex.release()

    def Consume(self):
        self.mutex.acquire()
        element=self.events.pop()
        self.mutex.release()
        return element

    def Consume_All(self):
        self.mutex.acquire()
        elements=self.events
        self.events=[]
        self.mutex.release()
        return elements

#######################################################################
# Exception EXIT_REPORTING_EXCEPTION
#######################################################################
class EXIT_REPORTING_EXCEPTION(Exception):
    def __init__(self):
        pass
#######################################################################
# Class MON_WRAPPER
#######################################################################
class MON_WRAPPER:
    def __init__(self,cmd):
        self.cmd=cmd
        self.start_callback=lambda pid,rpc,session: None
        self.simulation_events=[]
        self.polled_events=[]
        self.simulation_event_queue=PRODUCER_CONSUMER_QUEUE()
        self.poll_interval=20 # 5 seconds
        self.start_time=time.time()
        self.seconds_since_last_poll=0
        self.session=int(os.environ["PMON_SESSION_ID"])
        self.rpc_connection=None
        self.rpc_connection_active=0
        self.exited=0

    def run(self):
        # set up signal handling
        signal.signal(signal.SIGTERM, self.kill)
        signal.signal(signal.SIGINT, self.kill)
        
        # make the reporting thread
        self.reporting_thread=threading.Thread()
        self.reporting_thread.run=self.reporting
        self.reporting_thread.start()

        # launch the process and capture the pid
        pipe=popen2.Popen4(self.cmd)
        self.input,self.pid=pipe.fromchild,pipe.pid

        # report that we have started the program
        rpc_connection=None
        try:
            rpc_connection=CONNECT.mon_client()
            self.report_start(self.pid,rpc_connection,self.session);
        except:
            print "MON_WRAPPER_ENGINE: Failed to report start"
            rpc_connection.close()
        else:
            rpc_connection.close()

        # parse and pipe input
        while 1:
            line=""
            try:
                line=self.input.readline()
            except IOError:
                continue
            if line=="": break # cmd died, break to find out why
            print line,
            sys.stdout.flush()
            line=line.strip()
            for regexp,match_func,act_func in self.simulation_events:
                match=regexp.match(line)
                try:
                    event_data=match_func(match)
                except:
                    pass
                else:
                   self.simulation_event_queue.Produce((act_func,event_data))
        # end of process detect exit code
        exit_code=pipe.wait()
        print "\ngot exit code %d - should exit in at most a couple of seconds"%exit_code
        self.simulation_event_queue.Produce((self.report_exit,(exit_code,time.time())))
        self.reporting_thread.join()

    def kill(self,signal_number,stack_frame):
        if not self.exited:
            self.exited=1
            os.kill(self.pid,signal_number)

    def report_start(self,pid,rpc,session):
        rpc.Update_Status(session,"start_time",self.start_time)
        rpc.Update_Status(session,"pid",self.pid)
        rpc.Remove_Status_If_Exists(session,"exit_code")
        rpc.Remove_Status_If_Exists(session,"exit_time")
        rpc.Remove_Status_If_Exists(session,"total_duration")
        self.start_callback(pid,rpc,session)
        
    def report_exit(self,pid,rpc,session,data):
        exit_code,exit_time=data
        try:
            rpc.Update_Status(session,"exit_code",exit_code)
            rpc.Update_Status(session,"exit_time",exit_time)
            rpc.Update_Status(session,"total_duration",exit_time-self.start_time)
            rpc.Remove_Status_If_Exists(session,"pid")
        except:
            pass
        raise EXIT_REPORTING_EXCEPTION
    
    def reporting(self):
        try:
            while 1:
                time.sleep(.5)
                self.seconds_since_last_poll+=.5
                if self.seconds_since_last_poll<self.poll_interval and not self.exited: continue
                self.seconds_since_last_poll=0
                print "MON ------ Polling ----------------"
                try:
                    if not self.rpc_connection_active:
                        print "MON ------- Connection inactive reconnecting ---------"
                        self.rpc_connection=CONNECT.mon_client(10)
                        self.rpc_connection_active=1
                    #rpc_connection=DUMMY_RPC()
                except:
                    print "MON ------ Failed to get rpc_connection for reporting --------"
                    traceback.print_exc()
                    self.rpc_connection.close()
                    self.rpc_connection=None
                    self.rpc_connection_active=0
                    time.sleep(10)
                    pass
                else:
                    events=self.simulation_event_queue.Consume_All()
                    for event_function,event_data in events:
                        #print "events func=%s data=%s"%(repr(event_function),repr(event_data))
                        try:
                            ret=event_function(self.pid,self.rpc_connection,self.session,event_data)
                            #print "ret=%s"%ret
                        except EXIT_REPORTING_EXCEPTION:
                            raise EXIT_REPORTING_EXCEPTION
                        except:
                            print "MON ------ Error doing reporting event -------"
                            traceback.print_exc()
                            self.rpc_connection=None
                            self.rpc_connection_active=0
                            continue
                    for poll_function in self.polled_events:
                        try:
                            poll_function(self.pid,self.rpc_connection,self.session)
                        except:
                            print "MON ------ Error doing poll function ------"
                            traceback.print_exc()
        except EXIT_REPORTING_EXCEPTION:
            self.rpc_connection.close()
            self.rpc_connection=None
            pass
