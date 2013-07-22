#!/usr/bin/python
import time
import traceback
import os
import sys

from pd.common import CONNECT
from helpers import *



def Read_Back(raw_data):
    #print "GOT RAW DATA "+repr(data)
    seqid,exception,data=raw_data
    if data[0]=="SEND":
         type,sender,to,data=data
         if data[0]=="MESSAGE":
             sys.stdout.write("\r")
             log_formatted(sys.stdout,sender,to)
             sys.stdout.write(data[1])
             sys.stdout.write("\n")
             print "Cmd> ",
         elif data[0]=="PICTURE":
             sys.stdout.write("\r")
             log_formatted(sys.stdout,sender,to)
             tmpname=os.tmpnam()+".jpg"
             sys.stdout.write("Image:    %s\n"%tmpname)
             print "Cmd> ",

             pid=os.fork()
             if pid==0:
                 tmpname=os.tmpnam()+".jpg"
                 open(tmpname,"wb").write(data[1])
                 os.execvp("display",("display","-title","Picture from %s"%(format_instance(sender)),tmpname))
    

client=0
#client=CONNECT.send_client()
try:
    #foo=a.A()
    client=CONNECT.send_client()
    client.Set_Callback(Read_Back)
    #client.Register_Client(os.environ["USER"],os.environ["HOSTNAME"])
    
    while 1:
        print "Cmd> ",
        cmd=sys.stdin.readline()
        cmd=cmd[:-1].split(" ")
        if cmd[0]=="help":
            print "Commands:"
            print "---------"
            print "who"
            print "send <username> <msg>"
            print "help"
        if cmd[0]=="who":
            users=client.Users()
            for user,clients in users.items():
                print "%s -> %s"%(user,clients)
        if cmd[0]=="send":
            recipient=cmd[1]
            client.Send_Text([recipient]," ".join(cmd[1:]))
except:
    traceback.print_exc()
    client.close()
else:
    client.close()

#    #print client.Users()
#    #client.Send_Message(["aselle"],"fuck you")
#    #print "Foo"
#    #client.Send_Message(["aselle","shinar"],"fuck you")
#    #client.Test_Except()
#    #while 1:
#    #    client.Function_Call("Host_List",tuple())
#        
