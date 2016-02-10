#!/usr/bin/python

from pd.common import CONFIG
from pd.common import SOCKET
import os
import mutex
import time
import threading
import pickle
        
class SERVER:
    def __init__(self):
        self.sessions={}
        self.next_id=1
        # Define RPC interface
        self.mutex=threading.Lock()
        self.commands=["Register_Client","Session_Info","Session_List","Create_Session","Delete_Session","Label_Session","Update_Status","Remove_Status_If_Exists"]
        self.backup_interval=30
        # Save loop
        self.session_file=CONFIG.pdmon_session_file
        self.backup_thread=threading.Thread()
        self.backup_thread.run=self.Backup
        self.backup_thread.start()
        # read in data
        try:
            self.sessions=pickle.load(open(self.session_file,"r"))
            if(len(self.sessions.keys())>0):
                self.next_id=max(self.sessions.keys())+1
        except:
            pass
        print "MON_SERVER: Next id starts at %d"%self.next_id

    def Client_Connect(self,x):
        pass
        
    def Client_Disconnect(self,x):
        pass

    def Register_Client(self,client_id,user,host):
        pass

    # private
    def Validate_Session_Id(self,session_id):
        if type(session_id)!=int: raise SOCKET.COMMAND_EXCEPTION("Invalid session id")
        elif not self.sessions.has_key(session_id): raise SOCKET.COMMAND_EXCEPTION("Invalid session id %d"%session_id)

    # PUBLIC ROUTINES
    def Session_Info(self,client,session_id):
        self.Validate_Session_Id(session_id)
        return self.sessions[session_id]
    
    def Session_List(self,client):
        return self.sessions
        
    def Create_Session(self,client,username):
        session_id,directory=None,None
        session_id=self.next_id
        self.next_id+=1
        info={"id":session_id, "label": "<unnamed>","username": username,"created_date":time.time(),"last_update": 0,"user_status":{}}
        self.sessions[session_id]=info
        return info

    def Delete_Session(self,client,session_id):
        self.Validate_Session_Id(session_id)
        del self.sessions[session_id]

    def Label_Session(self,client,session_id,label):
        self.Validate_Session_Id(session_id)
        self.sessions[session_id]["label"]=label

    def Update_Status(self,client,session_id,key,value):
        self.Validate_Session_Id(session_id)
        self.sessions[session_id]["user_status"][key]=value
        self.sessions[session_id]["last_update"]=time.time()

    def Remove_Status_If_Exists(self,client,session_id,key):
        self.Validate_Session_Id(session_id)
        try:
            del self.sessions[session_id]["user_status"][key]
        except:
            pass

    def Backup(self):
        while 1:
            time.sleep(self.backup_interval)
            print "Backing up..."
            try:
                self.mutex.acquire()
                pickle.dump(self.sessions,open(self.session_file,"w"))
            finally:
                self.mutex.release()


import socket
if __name__ == "__main__":
    server=SERVER()
    SOCKET.SERVER(socket.gethostbyname(CONFIG.pdmon_server_host),CONFIG.pdmon_server_port,server)
#    SOCKET.SERVER(socket.gethostbyname(CONFIG.pdmon_server_host),CONFIG.pdmon_server_port,server,
#                  (CONFIG.server_private_key_file,CONFIG.server_certificate_file,CONFIG.ca_certificate_file))
