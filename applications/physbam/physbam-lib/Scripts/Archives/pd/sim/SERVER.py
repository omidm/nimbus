#!/usr/bin/python

# Host representation: hostname, max_cpus, max_memory
# Session representation: id, label, username, created_date, state={active,inactive,done} machine, claim_id, memory, cpus, user_status (where state and status are dictionaries)

from pd.common import CONFIG
from pd.common import SOCKET
import os
import mutex
import time
import threading
        
class SERVER:
    def __init__(self):
        self.session_directory=CONFIG.session_directory
        self.sessions={}
        self.next_id=1
        self.server_id="pd/%s"%os.environ["HOSTNAME"]
        self.hosts_client=client=SOCKET.CLIENT(CONFIG.pdhosts_server_host,CONFIG.pdhosts_server_port,
                  (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file))

        print "our test %s"%self.hosts_client.Host_List()
        self.Read_All_Sessions()

        # Define RPC interface
        self.mutex=threading.Lock()
        self.commands=["Session_Info","Session_List","Create_Session","Activate_Session","Deactivate_Session","Label_Session","Update_State","Host_List","Update_Status","Remove_Status_If_Exists","Session_Directory"]

    # PRIVATE ROUTINES
    def Read_Session(self,session_id):
        # read main status
        try:
            info=eval(open(os.path.join(self.session_directory,str(session_id),"etc","info.py")).read())
            self.sessions[info["id"]]=info
        except:
            pass

    def Read_All_Sessions(self):
        for directory in os.listdir(self.session_directory):
            self.Read_Session(directory)
        self.next_id=reduce(max,self.sessions.keys()+[self.next_id])+1

    def Write_Session(self,session_id):
        open(os.path.join(self.session_directory,str(session_id),"etc","info.py"),"w").write(repr(self.sessions[session_id]))


    def Validate_Session_Id(self,session_id):
        if type(session_id)!=int: raise SOCKET.COMMAND_EXCEPTION("Invalid session id")
        elif not self.sessions.has_key(session_id): raise SOCKET.COMMAND_EXCEPTION("Invalid session id %d"%session_id)

    # PUBLIC ROUTINES
    def Session_Info(self,session_id):
        self.Validate_Session_Id(session_id)
        return self.sessions[session_id]
    
    def Session_List(self):
        return self.sessions
        
    def Create_Session(self,username,memory,cpus):
        session_id,directory=None,None
        while 1:
            session_id=self.next_id
            self.next_id+=1
            directory=os.path.join(self.session_directory,str(session_id))
            if not os.path.exists(directory): break
        etc_directory=os.path.join(directory,"etc")
        info={"id":session_id, "label": "<unnamed>","username": username,"created_date":time.time(),"state":"inactive","machine":None,"claim_id":None,"memory":memory,"cpus":cpus,"user_status":{}}
        os.umask(0)
        os.mkdir(directory,01775) # create directory
        os.mkdir(etc_directory,01775) # create etc directory
        self.sessions[session_id]=info
        self.Write_Session(session_id)
        return info

    def Activate_Session(self,session_id,desired_hostname):
        self.Validate_Session_Id(session_id)
        if self.sessions[session_id]["machine"]!=None:
            raise SOCKET.COMMAND_EXCEPTION("session is already bound to machine %s"%self.sessions[session_id]["machine"])
        if self.sessions[session_id]["state"]=="active":
            raise SOCKET.COMMAND_EXCEPTION("session already activated but no machine: PANIC")
        claim_id=self.hosts_client.Claim_Host(desired_hostname,self.server_id,self.sessions[session_id]["username"],self.sessions[session_id]["cpus"],self.sessions[session_id]["memory"])
        self.sessions[session_id]["machine"]=desired_hostname
        self.sessions[session_id]["claim_id"]=claim_id
        self.sessions[session_id]["state"]="active"
        self.Write_Session(session_id)

    def Deactivate_Session(self,session_id,state):
        self.Validate_Session_Id(session_id)
        if state == "active": raise SOCKET.COMMAND_EXCEPTION("state must not be active")
        if self.sessions[session_id]["state"]!="active": raise SOCKET.COMMAND_EXCEPTION("session already inactive")
        if self.sessions[session_id]["machine"]==None: raise SOCKET.COMMAND_EXCEPTION("session is  not bound to machine but session is inactive: PANIC")
        self.hosts_client.Release_Host(self.sessions[session_id]["machine"],self.sessions[session_id]["claim_id"])
        self.sessions[session_id]["claim_id"]=None
        self.sessions[session_id]["machine"]=None
        self.sessions[session_id]["state"]=state
        self.Write_Session(session_id)

    def Label_Session(self,session_id,label):
        self.Validate_Session_Id(session_id)
        self.sessions[session_id]["label"]=label
        self.Write_Session(session_id)

    def Update_State(self,session_id):
        self.Validate_Session_Id(session_id)
        self.Read_Session(session_id)

    def Host_List(self):
        return self.hosts_client.Host_List()

    def Update_Status(self,session_id,key,value):
        self.Validate_Session_Id(session_id)
        self.sessions[session_id]["user_status"][key]=value
        self.Write_Session(session_id)

    def Remove_Status_If_Exists(self,session_id,key):
        self.Validate_Session_Id(session_id)
        try:
            del self.sessions[session_id]["user_status"][key]
            self.Write_Session(session_id)
        except:
            pass

    def Session_Directory(self,session_id):
        self.Validate_Session_Id(session_id)
        directory=os.path.join(self.session_directory,str(session_id))
        return directory

import socket
if __name__ == "__main__":
    server=SERVER()
    SOCKET.SERVER(socket.gethostbyname(CONFIG.pdsim_server_host),CONFIG.pdsim_server_port,server,
                  (CONFIG.server_private_key_file,CONFIG.server_certificate_file,CONFIG.ca_certificate_file))
