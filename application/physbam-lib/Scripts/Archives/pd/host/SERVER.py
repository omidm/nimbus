#!/usr/bin/python
from pd.common import CONFIG
from pd.common import SOCKET
import os
import mutex
import time
import threading
import socket

# Host representation: dictionary having hostname, user

class HOSTS_SERVER:
    def __init__(self):
        self.hosts_filename=CONFIG.hosts_filename
        self.commands=["Host_List","Add_Host","Delete_Host","Set_User"]
        self.mutex=threading.Lock()
        self.next_claim_id=1
        self.Read_Host_List()

    # Public
    def Read_Host_List(self):
        self.hosts=eval(open(self.hosts_filename).read())

    def Write_Host_List(self):
        open(self.hosts_filename,"w").write(repr(self.hosts))

    # Public
    def Host_List(self):
        return self.hosts

    def Add_Host(self,host):
        if self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Host already added")
        self.hosts[host]={"user":None}
        self.Write_Host_List()
    
    def Delete_Host(self,host):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        del self.hosts[host]
        self.Write_Host_List()
    
    def Lock_Host(self,host):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        if self.hosts[host]["locked"]: raise SOCKED.COMMAND_EXCEPTION("Host already locked")
        self.hosts[host]["locked"]=1
        self.Write_Host_List()
    
    def Unlock_Host(self,host):
        if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
        if not self.hosts[host]["locked"]: raise SOCKED.COMMAND_EXCEPTION("Host already unlocked")
        self.hosts[host]["locked"]=0
        self.Write_Host_List()
    
    def Set_User(self,host,user):
       if not self.hosts.has_key(host): raise SOCKET.COMMAND_EXCEPTION("Invalid host")
       self.hosts[host]["user"]=user
    
if __name__ == "__main__":
    server=HOSTS_SERVER()
    SOCKET.SERVER(socket.gethostbyname(CONFIG.pdhosts_server_host),CONFIG.pdhosts_server_port,server,
                  (CONFIG.server_private_key_file,CONFIG.server_certificate_file,CONFIG.ca_certificate_file))
