#!/usr/bin/python
from pd.common import SOCKET
from pd.common import CONNECT
import sys
import time
import os
import socket

filename=None

# get arguments
try:
    if sys.argv[1]=='-f':
        if len(sys.argv)<4: raise Exception
        filename=sys.argv[2]
        sys.argv=sys.argv[0:1]+sys.argv[3:]

    usernames=sys.argv[1:]
    if len(usernames)<1: raise Exception
except:
    print "Usage: %s [-f file] username [username ...]"%sys.argv[0]
    sys.exit(1)

if not filename:
    print "Select image rectangle to send"
    
    filename=os.tmpnam()+".jpg"
    os.system("import %s -quality 100"%filename)

try:
    data = open(filename, "rb").read()
except IOError:
    print "Image file %s not found" % imageFile
    sys.exit(0)

# try to send it
client=0
try:
    client=CONNECT.send_client()
    client.Send_Picture(usernames,data)
except SOCKET.COMMAND_EXCEPTION,e:
    print "ERROR: %s"%e
    client.close()
    sys.exit(1)
else:
    client.close()

    

