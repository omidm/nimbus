#!/usr/bin/python
import pd.common.SOCKET
import pd.common.CONNECT
import os
import sys
import dialog

def Print_Usage():
    print "Usage:  LOCK_HOST.py [-f] host_name\n  -f forces the lock even if there are claims on the machine"
    sys.exit(0)

hostname=""
force=0
if len(sys.argv)==2:
    if sys.argv[1][0]=="-": Print_Usage();
    hostname=sys.argv[1]
elif len(sys.argv)==3:
    if sys.argv[1]!="-f": Print_Usage();
    force=1
    hostname=sys.argv[2]
else: Print_Usage()

client=pd.common.CONNECT.host_client()
hosts=client.Host_List()
if not hosts.has_key(hostname):
    print "Host \"%s\" does not exist"%hostname
    sys.exit(0)

if len(hosts[hostname]["claims"])!=0 and not force:
    print "The host has the following claims:"
    for i in hosts[hostname]["claims"]:
        print hosts[hostname]["claims"][i]
    print "  Use -f to force the lock on this host"
    sys.exit(0)

print "Locking host %s"%hostname
client.Lock_Host(hostname)
