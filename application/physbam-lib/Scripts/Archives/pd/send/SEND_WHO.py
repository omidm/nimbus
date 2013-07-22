#!/usr/bin/python
from pd.common import SOCKET
from pd.common import CONNECT
import sys
import time
import os
import socket

# try to send it
client=0
try:
    client=CONNECT.send_client()
    users=client.Users()
    print "Logged in users:"
    for user,clients in users.items():
        print "%s -> %s"%(user,clients)
                
except SOCKET.COMMAND_EXCEPTION,e:
    print "ERROR: %s"%e
    client.close()
    sys.exit(1)
else:
    client.close()
