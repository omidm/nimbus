#!/usr/bin/python
import os
import sys
import re

cmd=sys.argv[1]

def parse_path(raw_path):
    path=raw_path
    if raw_path[0]!="/":
        pwd=os.environ['PWD'] # use pwd so we get /solver/vol1/v/foo instead of /solver/vol1/v--foo
        path=os.path.normpath(pwd+"/"+path)
        print "renormalized path to %s"%path
    r=re.compile("^/solver/(vol[0-9])/.+$")
    match=r.match(path)
    if not match:
        print "path does not start properly with /solver/vol{0,1,2,3}"
        sys.exit(1)
    volume=match.group(1)
    return volume,path

from pd.common import CONNECT

if cmd=="create":
    volume,path=parse_path(sys.argv[2])
    disk=CONNECT.disk_client(volume)
    print disk.Create(path)
    print disk.Set_Quota(path,"100G")
elif cmd=="destroy":
    volume,path=parse_path(sys.argv[2])
    disk=CONNECT.disk_client(volume)
    print disk.Destroy(path)
elif cmd=="setquota":
    volume,path=parse_path(sys.argv[2])
    quota=sys.argv[3]
    disk=CONNECT.disk_client(volume)
    print disk.Set_Quota(path,quota)
elif cmd=="list":
    volume=sys.argv[2]
    disk=CONNECT.disk_client(volume)
    print "%-40s %10s %10s %10s"%("Volume","Used","Avail","Repr")
    print " "
    for vol,used,avail,rep,mount in disk.List():
        print "%-40s %10s %10s %10s"%(vol,used,avail,rep)
else:
    print """
Usage:
    vfs create <path>
    vfs destroy <path>
    vfs setquota <path> <quota>
    vfs list <volume>
"""
    

