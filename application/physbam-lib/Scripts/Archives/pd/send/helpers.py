#!/usr/bin/python 
import time
def format_instance(x):
    username,client=x
    addr,port=client
    return "%s@%s:%d"%(username,addr,port)

def formatted_tos(x):
    users={}
    for i in x:
        username,client=i
        users[username]=1
    return users.keys()

def log_formatted(fp,sender,to):
    fp.write("--------------------------------------------------------------------\n")
    
    fp.write("Date:     %s\n"%time.ctime())
    fp.write("From:     %s\n"%format_instance(sender))
    tos=formatted_tos(to)
    fp.write("To:       %s\n"%tos[0])
    for other in tos[1:]:
        fp.write("          %s\n"%other)

