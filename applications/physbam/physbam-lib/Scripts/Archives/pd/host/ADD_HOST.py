#!/usr/bin/python
import pd.common.SOCKET
import pd.common.CONNECT
import os
import sys

if len(sys.argv)<2 or sys.argv[1]=="-h":
    print "Usage:  ADD_HOST.py host_name\n"
    sys.exit(0)

client=pd.common.CONNECT.host_client().Add_Host(sys.argv[1])
