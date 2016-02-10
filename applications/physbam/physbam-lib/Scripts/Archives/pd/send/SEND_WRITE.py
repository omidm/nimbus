#!/usr/bin/python
from pd.common import SOCKET
import pygtk
import gtk
import smtplib
import os
import socket
import sys
import pango
import traceback

class SEND_WRITE_WINDOW:
    def compute_window_dimensions(self,message):
        lines=message.split("\n")
        max_line_length=reduce(max,map(len,lines))
        return (max(512,max_line_length*8+90),max((len(lines)+3)*14+90,150))

    def cancel(self,w):
        self.window.destroy()

    def send(self,w):

        
        text=self.textbuffer.get_text(self.textbuffer.get_start_iter(),self.textbuffer.get_end_iter())
        try:
            self.client.Send_Text(self.to,text)
        except:
            dlg=gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_OK)
            dlg.set_markup("Could not send error: \n%s"%traceback.format_exc)
            dlg.show()
        self.window.destroy()
    
    def __init__(self,client,to,message):
        self.to=to
        self.message=message
        self.client=client

        # set gui
        self.tree=gtk.glade.XML("/usr/local/adm/pd/pd/send/send.glade","sendMessageWindow")
        signals={"on_cancelButton_clicked":self.cancel,
                 "on_sendButton_clicked":self.send}
        self.tree.signal_autoconnect(signals)

        # get widgets
        self.window=self.tree.get_widget("sendMessageWindow")
        self.toField=self.tree.get_widget("tofield")
        self.text=self.tree.get_widget("text")


        self.toField.set_text("\n".join(self.to))
        # resize
        dimensions=self.compute_window_dimensions(message)
        self.window.resize(dimensions[0],dimensions[1])  
        
        self.textbuffer = self.text.get_buffer()
        self.text.modify_font(pango.FontDescription("Courier Bold 10"))
        self.text.modify_base(gtk.STATE_NORMAL,gtk.gdk.color_parse("#443355"))
        self.text.modify_text(gtk.STATE_NORMAL,gtk.gdk.color_parse("white"))

        self.textbuffer.set_text("\n\n"+message);
        self.textbuffer.place_cursor(self.textbuffer.get_start_iter());
                
        self.window.show()

