#!/usr/bin/python
import sys
import os
import dialog
import time
import CLIENT_LIBRARY


guys=0
curr=time.time()
while 1:
    guys+=1
    CLIENT_LIBRARY.client.Update_Status(2,"foo","bar")
    if guys % 200==0:
        old_time=curr
        curr=time.time()
        print "%s %s"%(repr(guys),repr(curr))
        print "Query/Sec = %f"%(float(guys)/(curr-old_time))
        guys=0
        
