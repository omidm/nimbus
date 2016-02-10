#!/usr/bin/python
import time
import traceback
import os
import sys
import pygtk
import gtk
import gobject

from pd.common import CONNECT
from SEND_READ import SEND_READ_WINDOW
from SEND_WRITE import SEND_WRITE_WINDOW
from SEND_VIEW_IMAGE import SEND_VIEW_IMAGE_WINDOW
from helpers import *

# Figure out good loggind directory
year,month,day,hour,minute,second,wday,yday,dst=time.localtime()
hostname=os.environ["HOSTNAME"].split(".")[0]
home=os.environ["HOME"]
fname="%04d-%02d-%02d_%s_%d"%(year,month,day,hostname,os.getpid())
logdir=os.path.join(home,".send")
logpath=os.path.join(logdir,fname)
# Make log file
if not os.path.isdir(logdir): os.mkdir(logdir)
log_fp=open(logpath,"w+")
log_fp.write("Starting log...\n")

def Read_Back(raw_data):
    seqid,exception,data=raw_data
    if data[0]=="SEND":
         type,sender,to,data=data
         if data[0]=="MESSAGE":
             log_formatted(log_fp,sender,to)
             log_fp.write(data[1])
             log_fp.write("\n")
             log_fp.flush()

             try:
                 gtk.gdk.threads_enter()
                 SEND_READ_WINDOW(client,format_instance(sender),map(format_instance,to),data[1])
             finally:
                 gtk.gdk.threads_leave()
             
         elif data[0]=="PICTURE":
             try:
                 gtk.gdk.threads_enter()
                 tmpname=os.tmpnam()+".jpg"
                 fp=open(tmpname,"wb")
                 fp.write(data[1])
                 fp=None
                 SEND_VIEW_IMAGE_WINDOW(client,format_instance(sender),map(format_instance,to),tmpname)
             finally:
                 gtk.gdk.threads_leave()

             log_formatted(log_fp,sender,to)
             log_fp.write("Image:    %s\n"%tmpname)
             log_fp.flush()
             
    

label=None

client=None
try:
    client=CONNECT.send_client()
    gtk.gdk.threads_init()
    client.Set_Callback(Read_Back)
    #SEND_WRITE_WINDOW(client,["aselle"],"hi")
    gtk.main()
    

except:
    traceback.print_exc()
    client.close()
else:
    client.close()
