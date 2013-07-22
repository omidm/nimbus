#!/usr/bin/python
import pygtk
import gtk
import gtk.glade
import smtplib
import time
import os
import sys
#import SEND_SEND_MESSAGE
from SEND_WRITE import SEND_WRITE_WINDOW

import pango
import random
import math
#from pd.send import UTIL

def hsv_to_rgb(h,s,v):
    h,s,v=float(h),float(s),float(v)
    hi=int(h/60)%6
    f=h/60-hi
    p=v*(1-s)
    q=v*(1-f*s)
    t=v*(1-(1-f)*s)
    if hi==0: return v,t,p
    elif hi==1: return q,v,p
    elif hi==2: return p,v,t
    elif hi==3: return p,q,v
    elif hi==4: return t,p,v
    elif hi==5: return v,p,q



def get_text(parent,label):
    dlg=gtk.Dialog(label,parent,gtk.DIALOG_MODAL,("Cancel",0,"Ok",1))
    label=gtk.Label(label)
    entry=gtk.Entry()
    dlg.vbox.pack_start(label,True,True,0)
    dlg.vbox.pack_start(entry,True,True,0)
    label.show()
    entry.show()
    
    if dlg.show()==0:
        return None
    else:
        return entry.get_text()

def username_only(addr):
    return addr.split("@")[0]

class SEND_READ_WINDOW:
    def __del__(self):
        #print "deleting window"
        pass
        
    def forward(self,widget):
        pass
        #text=get_text(self.window,"Forward to:")
        #if text:
        #    try:
        #        usernames=text.split(",")
        #        message="%s wrote:\n> "%self.user+self.message.replace("\n","\n> ")
        #        hosts=Users_To_Hosts(usernames)
        #        if len(hosts)>0:
        #            for user,host in hosts:
        #                UTIL.Send_Message(host,message)
        #    except:
        #        pass

    def compute_window_dimensions(self,message):
        lines=message.split("\n")
        max_line_length=reduce(max,map(len,lines))
        return (max(max_line_length*8+20,512),max(len(lines)*14+30,150))

    def reply(self,w):
        quoted_message="\n".join(map(lambda x: "> "+x,self.message.split("\n")))
            
        SEND_WRITE_WINDOW(self.client,[username_only(self.sender)],("Message wrote %s\n"%username_only(self.sender))+quoted_message)

    def close(self,w):
        self.window.destroy()
        
    def email(self,w):
        username=os.environ['USER']
        from_addr="%s@graphics.stanford.edu"%username_only(self.sender)
        to_addr=["%s@graphics.stanford.edu"%os.environ["USER"]]
        #print "Emailing from: %s to %s"%(from_addr,to_addr)
        msgtext="Date: %s\nFrom: %s\nTo:    %s\n\n%s"%(time.ctime(),self.sender,(", ".join(self.to)),self.message)
        smtplib.SMTP("gpo.stanford.edu").sendmail(from_addr,to_addr,"Subject: [Send] Message\n\n%s"%(msgtext))
        
    def replyClose(self,w):
        self.reply(w)
        self.window.destroy()

    def wordWrap(self,w):
        if w.get_active():
            self.text.set_wrap_mode(gtk.WRAP_WORD)
        else:
            self.text.set_wrap_mode(gtk.WRAP_NONE)
        
        
    def __init__(self,client,sender,to,message):
        sys.stdout.write('\a')
        sys.stdout.flush()

        self.client=client
        self.sender=sender
        self.to=to
        #self.hostname=hostname
        self.message=message


        self.tree=gtk.glade.XML("/usr/local/adm/pd/pd/send/send.glade","receiveMessageWindow")
        signals={"on_replyButton_clicked":self.reply,
                 "on_replyCloseButton_clicked":self.replyClose,
                 "on_emailButton_clicked":self.email,
                 "on_closeButton_clicked":self.close,
                 "on_wordWrapButton_toggled":self.wordWrap}
        
        self.tree.signal_autoconnect(signals)
        # get widgets

        self.window=self.tree.get_widget("receiveMessageWindow")
        self.dateField=self.tree.get_widget("datefield")
        self.fromField=self.tree.get_widget("fromfield")
        self.toField=self.tree.get_widget("tofield")
        self.text=self.tree.get_widget("text")


        # size window
        dimensions=self.compute_window_dimensions(message)
        self.window.resize(dimensions[0],dimensions[1])  

        # populate fields
        self.dateField.set_text(time.ctime())
        self.fromField.set_text(self.sender)

        self.to_usernames=[]
        for i in self.to:
            username=username_only(i)
            if username not in self.to_usernames:
                self.to_usernames.append(username)
        #print self.to_usernames
        self.toField.set_text("\n".join(self.to_usernames))

        textbuffer = self.text.get_buffer()
        self.text.set_editable(False);
        self.text.modify_font(pango.FontDescription("Luxi Mono 10"))
        hue=random.randint(0,360)
        r,g,b=map(int,hsv_to_rgb(hue,.9,50))
        self.text.modify_base(gtk.STATE_NORMAL,gtk.gdk.color_parse("#%02X%02X%02X"%(r,g,b)))
        r,g,b=map(int,hsv_to_rgb(hue,.2,255))
        self.text.modify_text(gtk.STATE_NORMAL,gtk.gdk.color_parse("#%02X%02X%02X"%(r,g,b)))
        textbuffer.set_text("%s"%(self.message))

