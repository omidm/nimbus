#!/usr/bin/python
from pd.common import SOCKET
from pd.common import CONNECT
import sys
import time
import os
import socket

# get arguments
try:
    executable,usernames=sys.argv[0],sys.argv[1:]
    if len(usernames)<1: raise Exception
except:
    print "Usage: %s username"%sys.argv[0]
    sys.exit(0)

# read message
if sys.stdin.isatty():
    print "Type message to send. (^d to send, ^c to cancel)"
    message=sys.stdin.read()
else:
    message=sys.stdin.read()

# try to send it
client=0
try:
    client=CONNECT.send_client()
    client.Send_Text(usernames,message)
except SOCKET.COMMAND_EXCEPTION,e:
    print "ERROR: %s"%e
    client.close()
    sys.exit(1)
else:
    client.close()
