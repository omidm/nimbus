#!/usr/bin/python
import pd.common.SOCKET
import pd.common.CONNECT
import os
import sys
try:
   hosts=sys.argv[1:]
   if len(hosts)<1: raise Exception
except:
   print "Usage: %s <host1> [host2 host3 ...]"%sys.argv[0]
   sys.exit(1)

client=pd.common.CONNECT.host_client()

for i in hosts:
    print "Deleting %s"%i
    client.Delete_Host(i)
