#!/usr/bin/python
from pd.common import CONFIG
from pd.common import SOCKET
import os
import mutex
import time
import threading
import socket

# Host representation: dictionary having hostname, user

class SEND_SERVER:
    def __init__(self):
        #self.hosts_filename=CONFIG.hosts_filename
        self.commands=["Register_Client","Users","Send_Text","Send_Picture"]
        self.lookup_hosts=True
        self.mutex=threading.Lock()
        self.next_claim_id=1
        self.users={}
        self.clients={}
        self.clientid_to_username={}

        #self.Read_Host_List()

    def Client_Connect(self,x):
        self.clients[x.host]=x
        
    def Client_Disconnect(self,x):
        del self.clients[x.host]
        if self.clientid_to_username.has_key(x.host):
            user=self.clientid_to_username[x.host]
            del self.clientid_to_username[x.host]
            self.users[user].remove(x.host)
            print "Unregistered user=%s client=%s"%(user,x.host)

    def Registered(self,client_id):
        return self.clientid_to_username.has_key(client_id)

    def Register_Client(self,client_id,user,host):
        if self.Registered(client_id): raise SOCKET.COMMAND_EXCEPTION("Connection already registered for user %s"%self.clientid_to_username[client_id])
        if not self.users.has_key(user): self.users[user]=[]
        self.clientid_to_username[client_id]=user
        self.users[user].append(client_id)
        print "Registered user=%s client=%s"%(user,client_id)

    def Send(self,client_id,users,data):
        if not self.Registered(client_id): raise SOCKET.COMMAND_EXCEPTION("Your client is not registered")
        not_found_users=[]
        users_and_clientids=[]
        for user in users:
            if not self.users.has_key(user) or len(self.users.keys())==0: not_found_users.append(user)
            else: users_and_clientids.extend(map(lambda x: (user,x),self.users[user]))
        if len(not_found_users)>0: raise SOCKET.COMMAND_EXCEPTION("No registration for users: %s"%",".join(not_found_users))
        for user,clientid in users_and_clientids:
            client=self.clients[clientid]
            client.queueWrite((-100,None,("SEND",(self.clientid_to_username[client_id],client_id),users_and_clientids,data)))

    def Send_Text(self,client_id,users,message):
        return self.Send(client_id,users,("MESSAGE",message))

    def Send_Picture(self,client_id,users,picture):
        return self.Send(client_id,users,("PICTURE",picture))
    
    def Users(self,client):
        return self.users

    
if __name__ == "__main__":
    server=SEND_SERVER()
    SOCKET.SERVER(socket.gethostbyname(CONFIG.pdsend_server_host),CONFIG.pdsend_server_port,server) #,(CONFIG.server_private_key_file,CONFIG.server_certificate_file,CONFIG.ca_certificate_file))
