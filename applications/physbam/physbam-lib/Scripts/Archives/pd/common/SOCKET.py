#!/usr/bin/python
import signal
import threading
import sys
import pickle
import socket
import select
import os
import time
import traceback
from OpenSSL import SSL
from OpenSSL import crypto

def daemonize():
    sys.stdin=None
    log=open("log.txt","w")
    sys.stdout=log
    sys.stderr=log
    PID = os.fork()
    if PID != 0:
        print "shutting down parent process:",PID
        sys.exit()
    os.setsid()
    PID = os.fork()
    if PID != 0:
        print "shutting down parent process:",PID
        sys.exit()

class CONNECTION_LOST_ERROR(Exception):
    def __init__(self,value=""):
        self.value=str(value)
        
    def __str__(self):
        return self.value

class THREAD_EXIT_ERROR(Exception):
    def __init__(self,value=""):
        self.value=str(value)
        
    def __str__(self):
        return self.value

class COMMAND_EXCEPTION(Exception):
    def __init__(self,value,detail=""):
        self.value=value
        self.detail=detail

    def __str__(self):
        return "%s %s"%(self.value,self.detail)

class READER:
    def __init__(self,delimiter="\n"):
        #self.socket=socket
        self.buffer=""
        self.delimiter=delimiter
        
    def get(self,sock):
        try:
           chunk=sock.recv(819200)
        except (SSL.WantReadError, SSL.WantWriteError, SSL.WantX509LookupError):
           pass
        except SSL.ZeroReturnError:
           raise CONNECTION_LOST_ERROR()
        except SSL.Error, errors:
           raise CONNECTION_LOST_ERROR(errors)
        else:
           if chunk=="": raise CONNECTION_LOST_ERROR
           self.buffer+=chunk
        items=[]
        while 1:
            index=self.buffer.find(self.delimiter)
            if index != -1:
                line=self.buffer[:index+len(self.delimiter)]
                self.buffer=self.buffer[index+len(self.delimiter):]
                items.append(line)
            else:
                break
        return items

class WRITER:
    def __init__(self):
        self.buffer=""
        self.mutex=threading.Lock()
        self.thread_exit=False # used by server

    def put(self,sock):
        if len(self.buffer)==0: return
        try:
            try:
                self.mutex.acquire()
                sent=0
                #print "trying to write '%s'"%self.buffer
                idx=self.buffer.find("\n")
                if idx!=-1:
                    sent=sock.send(self.buffer[:idx+1])
                else: sent=sock.send(self.buffer)
                #print "Wrote %d bytes " %sent
                if sent==0: raise CONNECTION_LOST_ERROR()
                self.buffer=self.buffer[sent:]
            finally:
                self.mutex.release()
        except (SSL.WantReadError, SSL.WantWriteError, SSL.WantX509LookupError):
           pass
        except SSL.ZeroReturnError:
           raise CONNECTION_LOST_ERROR()
        except SSL.Error, errors:
           raise CONNECTION_LOST_ERROR(errors)
           


class SERVER_RUN_CMD_THREAD(threading.Thread):
    def __init__(self,client_socket,client_host,server):
        threading.Thread.__init__(self)

        self.socket=client_socket
        self.host=client_host
        self.server=server

        self.reader=READER(".\n")
        self.writer=WRITER()
        self.start()

    def queueWrite(self,response): # response is of form (seqid,exception,text)
        #print "trying to queue"
        pickled_result=pickle.dumps(response)
        #print "sending... result=%s"%repr(response)
        self.writer.mutex.acquire()
        #print "got mutex"
        self.writer.buffer+=(pickled_result+"\n")
        #print "added to bufefr"
        self.writer.mutex.release()
        #print "release mutex"

    def run(self):
        # add myself to the connection list
        self.server.mutex.acquire()
        self.server.Client_Connect(self)
        self.server.mutex.release()
        
        
        try:
            while 1:
                socklist=[self.socket]
                write_socklist=socklist
                # if nothing needs to be written don't select on it
                self.writer.mutex.acquire()
                if len(self.writer.buffer)==0: write_socklist=[]
                if self.writer.thread_exit: raise CONNECTION_LOST_ERROR
                self.writer.mutex.release()
                # selec to see what we need to read
                #print "select put %s %s %s"%(repr(socklist),repr(write_socklist),repr(socklist))
                [read_ready,write_ready,except_ready]=select.select(socklist,write_socklist,socklist,.1)
                #print "select got %s %s %s"%tuple(map(repr,[read_ready,write_ready,except_ready]))
                if len(except_ready)>0:
                    print "got exception on socket"%except_ready[0]
                    break
                for r in read_ready:
                    # Get the python command object
                    pickle_chunks=self.reader.get(r)
                    for pickle_chunk in pickle_chunks:    
                        #print "pickle_chunk=%s"%pickle_chunk
                        thing=None
                        result=(None,None,None) # (seqid,exception, result)
                        seqid=-1 # unknown sequence id
                        try:
                            try: thing=pickle.loads(pickle_chunk[:-1])
                            except: raise COMMAND_EXCEPTION("Non-pickled encountered from %s"%repr(self.host))
                            # make sure command is of form (command, args)
                            command,args=None,None
                            try:
                                seqid,command,args=thing
                                #print "Client %s:%d Command %s%s"%(self.host[0],self.host[1],command,
                                #repr(args))
                            except: raise COMMAND_EXCEPTION("Invalid command packing")
                            # Make sure command exists
                            if not self.server.commands.has_key(command): raise COMMAND_EXCEPTION("Command \"%s\" not found"%repr(command))
                            try:
                                self.server.mutex.acquire()
                                print "%s: Ran command %s"%(self.host,command)
                                result=(seqid,None,self.server.commands[command](self.host,*args))
                            finally:
                                self.server.mutex.release()
                        except COMMAND_EXCEPTION,e:
                            result=(seqid,(e.value,e.detail),None)
                        except:
                            result=(seqid,("Other Exception",str(sys.exc_info()[1])+"\n"+"\n".join(traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_traceback))),None)
                        
                        self.queueWrite(result)
                for w in write_ready:
                    self.writer.put(w)
            

        except CONNECTION_LOST_ERROR:
            print "%s: Connection lost"%repr(self.host)
            pass
        print "%s: DISCONNECTED"%repr(self.host)
        # remove myself from the connection list
        self.server.mutex.acquire()
        self.server.Client_Disconnect(self)
        self.server.mutex.release()

        self.socket.shutdown(2) # shutdown both end
        self.socket.close()
        self.socket=None

        
# check to make sure certificate is good
# NOTE THIS IS OUTSIDE THE CLASSES otherwise circular reference count hell happens
def verify_client_certificate(connection,certificate,errnum,depth,ok):
    #print "Got certificate: %s ok=%d"%(certificate.get_subject(),ok)
    return ok

class HEARTBEAT(threading.Thread):
    def __init__(self,server):
        threading.Thread.__init__(self)
        self.server=server
        self.start()

    def run(self):
        while 1:
            print "Heartbeat at time %s"%(time.time())
            self.server.mutex.acquire()
            for host,conn in self.server.connections.items():
                conn.queueWrite((-100,None,("HEARTBEAT")))
                print "%s"%repr(host)
            self.server.mutex.release()
            time.sleep(5)
        
        

class SERVER:
    def __init__(self,bindhost,port,rpc_server,secure=None):
        self.rpc_server=rpc_server
        self.commands=dict(map(lambda x: (x,getattr(rpc_server,x)),rpc_server.commands))
        self.mutex=rpc_server.mutex
        self.exiting=False
        try:
            self.lookup_hosts=rpc.lookup_hosts
        except:
            self.lookup_hosts=False
        self.connections={}

        self.raw_listen_socket=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        # Make listen socket
        if secure:
            server_pkey_file,server_cert_file,ca_cert_file=secure
            ctx=SSL.Context(SSL.SSLv23_METHOD)
            ctx.set_verify(SSL.VERIFY_PEER,verify_client_certificate)
            ctx.use_privatekey_file(server_pkey_file)
            ctx.use_certificate_file(server_cert_file)
            ctx.load_verify_locations(ca_cert_file)
            self.listen_socket=SSL.Connection(ctx,self.raw_listen_socket)
        else:
            self.listen_socket=self.raw_listen_socket
        print "have bindhost=%s and port=%d"%(bindhost,port)
        self.listen_socket.bind((bindhost,port))
        self.listen_socket.listen(5)

        # socket
        def abort(signal,frame):
            self.listen_socket.close()
            sys.exit(0)
        signal.signal(signal.SIGHUP,abort)

        # Server loop
        #self.heartbeat=HEARTBEAT(self)
        try:
            while 1:
                (client_socket,addr)=self.listen_socket.accept()
                print "%s: CONNECTED"%repr(addr)
                try:
                    addr=socket.gethostbyaddr(addr[0])[0],addr[1]
                except:
                    pass
                SERVER_RUN_CMD_THREAD(client_socket,addr,self)
        finally:
            for i,thread in self.connections.items():
                thread.writer.mutex.acquire()
                thread.writer.thread_exit=True
                thread.writer.mutex.release()
                thread.join()
            self.listen_socket.close()
            #self.listen_socket.shutdown(2)

    def Client_Connect(self,x):
        self.connections[x.host]=x
        self.rpc_server.Client_Connect(x)
        
    def Client_Disconnect(self,x):
        del self.connections[x.host]
        self.rpc_server.Client_Disconnect(x)
        
class CLIENT_RUN_THREAD(threading.Thread):
    def __init__(self,hostname,port,secure,timeout):
        threading.Thread.__init__(self)

        self.hostname=hostname
        self.port=port
        self.secure=secure
        self.timeout=timeout
        
        self.reader=READER(".\n")
        self.writer=WRITER()
        self.response_condition=threading.Condition()
        self.waiting_seqid=0
        self.exit=False
        self.waiting_response=None
        self.callback=lambda x: None

        self.connect()

        self.start()

    def connect(self):
        # Initialize context
        self.socket=None
        if self.secure:
            client_pfile,client_cert_file,ca_cert_file=self.secure
            #server_pkey=crypto.load_privatekey(crypto.FILETYPE_PEM,server_ptext)
            #server_cert=crypto.load_certificate(crypto.FILETYPE_PEM,server_cert_text)
            #ca_cert=crypto.load_certificate(crypto.FILETYPE_PEM,ca_cert_text)
            ctx = SSL.Context(SSL.SSLv23_METHOD)
            ctx.set_verify(SSL.VERIFY_PEER,verify_client_certificate) # Demand a certificate
            #ctx.use_privatekey(server_pkey)
            #ctx.use_certificate(server_cert)
            #ctx.use_verify_locations(ca_cert)
            ctx.use_privatekey_file(client_pfile)
            ctx.use_certificate_file(client_cert_file)
            ctx.load_verify_locations(ca_cert_file)
            self.socket=SSL.Connection(ctx,socket.socket(socket.AF_INET,socket.SOCK_STREAM))
        else:
            self.socket=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
        if self.timeout and socket!=None:
            try:
                socket.settimeout(self.timeout)
            except: pass
        self.socket.connect((self.hostname,self.port))

        register_cmd=(0,"Register_Client",(os.environ["USER"],os.environ["HOSTNAME"]))
        self.queueCmd(register_cmd)

        
    def run(self):
        try: # wait for thread exit
            while 1:
                try: # wait for socket close
                    while 1:
                        socklist=[self.socket]
                        write_socklist=socklist
                        # if nothing needs to be written don't select on it
                        self.writer.mutex.acquire()
                        if len(self.writer.buffer)==0: write_socklist=[]
                        need_exit=self.exit
                        self.writer.mutex.release()
                        if need_exit:
                            raise THREAD_EXIT_ERROR
                        # selec to see what we need to read
                        [read_ready,write_ready,except_ready]=select.select(socklist,write_socklist,socklist,.1)
                        #print "select got %s %s %s"%tuple(map(repr,[read_ready,write_ready,except_ready]))
                        if len(except_ready)>0:
                            #print "got exception on socket"%except_ready[0]
                            break
                        for r in read_ready:
                            # Get the python command object
                            pickle_chunks=self.reader.get(r)
                        
                            for pickle_chunk in pickle_chunks:
                                 seqid,exception,value=None,None,None
                                 #try:
                                 data=pickle.loads(pickle_chunk)                    #print repr(data)
                                 seqid,exception,value=data
                                 #print "Got seqid=%d except=%s value=%s"%(seqid,repr(exception),repr(value))
                                 
                                 self.response_condition.acquire()
                                 #print "waiting seq id is %d and we got %d"%(self.waiting_seqid,seqid)
                                 if self.waiting_seqid==seqid:
                                     self.waiting_response=data
                                     data=None
                                     #print "in child have %s"%repr(self.waiting_response)
                                     self.response_condition.notify()
                                 self.response_condition.release()
                                 if seqid!=0 and data!=None: # else condition for non expected things
                                     self.callback(data)
                                     
                                 #except:
                                 #    raise COMMAND_EXCEPTION("Failed to read return value")                    
        
                        for w in write_ready:
                            self.writer.put(w)
        
                except CONNECTION_LOST_ERROR:
                    self.socket.close()
                    pass
                else:
                    self.socket.close()
                    
                self.reconnect()
                
        except THREAD_EXIT_ERROR:
            pass

    def reconnect(self):
        while 1:
            print "Connection lost, trying to reconnect in 1 second"
            time.sleep(1)
            try:
                self.connect()
                break # if successful get out of reconnect loop
            except:
                pass
        

    def queueCmd(self,cmd): # response is of form (seqid,exception,text)
        pickled_cmd=pickle.dumps(cmd)
        #print "sending... cmd=%s"%repr(cmd)
        self.writer.mutex.acquire()
        self.writer.buffer+=(pickled_cmd+"\n")
        self.writer.mutex.release()

        

class CLIENT:
    # usage is (hostname, port, (privatekey_text,certificate_text,ca_certificate_text))
    def __init__(self,hostname,port,secure=None,timeout=0):
        self.local_methods={"close":self.close,"__nonzero__":None,"__del__":self.__del__,"Set_Callback":self.Set_Callback}
        self.seqid=100
        self.runthread=CLIENT_RUN_THREAD(hostname,port,secure,timeout)
        #self.reader=SOCKET_READER(".\n")

    def close(self):
        #print "closing guy"
        self.runthread.writer.mutex.acquire()
        self.runthread.exit=True
        self.runthread.writer.mutex.release()

    def __del__(self):
        self.close()

    def Test_Except(self):
        raise COMMAND_EXCEPTION("a","b")

    def Set_Callback(self,x):
        self.runthread.callback=x

    def Function_Call(self,command,args):
        #return
        
        pickled_command=pickle.dumps((self.seqid,command,args))
        self.runthread.queueCmd((self.seqid,command,args))

        # wait for response from communication thread
        self.runthread.response_condition.acquire()
        self.runthread.waiting_seqid=self.seqid
        self.runthread.response_condition.wait()
        seqid,exception,value=self.runthread.waiting_response
        self.runthread.waiting_response=None
        self.runthread.response_condition.release()
        
        # next sequence id
        self.seqid+=1

        if exception!=None:
            raise COMMAND_EXCEPTION(exception[0],exception[1])

        return value

    def __getattr__(self,method_name):
        #print "__getattr__(%s)"%method_name
        if self.local_methods.has_key(method_name): return self.local_methods[method_name]
        return lambda *x: self.Function_Call(method_name,x)

