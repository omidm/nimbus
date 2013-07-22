#!/usr/bin/python
import pygtk
import gtk
import gtk.glade
import smtplib
import time
import os
import sys
import pango
import random
import math
import shutil

# TODO: factor
def username_only(addr):
    return addr.split("@")[0]


class SEND_VIEW_IMAGE_WINDOW:
    def close(self,w):
        self.window.destroy()

    def file_ok(self,w):
        fname=self.filesel.get_filename()
        shutil.copy(self.fname,fname)
        self.filesel.destroy()
        

    def save(self,w):
        self.filesel=gtk.FileSelection(title="Save image to...")
        filename=self.filesel.get_filename()
        self.filesel.connect("destroy",lambda w: self.filesel.destroy())
        self.filesel.ok_button.connect("clicked",self.file_ok)
        self.filesel.cancel_button.connect("clicked",lambda w:self.filesel.destroy())
        self.filesel.set_filename("image.jpg")
        self.filesel.show()
        
    def __init__(self,client,sender,to,fname):
        self.client=client
        self.sender=sender
        self.to=to
        self.fname=fname

        self.tree=gtk.glade.XML("/usr/local/adm/pd/pd/send/send.glade","receiveImageWindow")
        signals={"on_saveButton_clicked":self.save,
                 "on_closeButton_clicked":self.close}
        self.tree.signal_autoconnect(signals)
        # get widgets
        self.dateField=self.tree.get_widget("datefield")
        self.fromField=self.tree.get_widget("fromfield")
        self.toField=self.tree.get_widget("tofield")
        self.image=self.tree.get_widget("image")
        self.window=self.tree.get_widget("receiveImageWindow")

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

        # load image
        self.image.set_from_file(fname)
        pixbuf=self.image.get_pixbuf()
        size=(pixbuf.get_width(),pixbuf.get_height())
        self.window.resize(size[0]+60,size[1]+140)
        self.window.show()
        
