#!/usr/bin/python
import pd.common.SOCKET
import pd.common.CONNECT
import os
import sys
import dialog

def Print_Usage():
    print "Usage:  UNLOCK_HOST.py host_name\n"
    sys.exit(0)

hostname=""
force=0
if len(sys.argv)==2:
    if sys.argv[1][0]=="-": Print_Usage();
    hostname=sys.argv[1]
else: Print_Usage()

client=pd.common.CONNECT.host_client()

print "Unlocking host %s"%hostname
client.Unlock_Host(hostname)
