#!/usr/bin/python
import os
import sys


import curses
import curses.wrapper
from curses.textpad import Textbox

class Widget:
    def __init__(self,parent,x,y,w,h,title):
        self.xprime=x
        self.yprime=y
        self.wprime=w
        self.hprime=h
        self.x=x+1
        self.y=y+1
        self.h=h-2
        self.w=w-2
        self.title=title
        self.pwin=parent.subwin(self.hprime,self.wprime,self.yprime,self.xprime)
        self.win=parent.subwin(self.h,self.w,self.y,self.x)
        
        open("/tmp/fdsaf","w").write(repr((self.h,self.w,self.y,self.x)))
        
        
    def draw(self,is_active):
        open("/tmp/ha","w").write("ha")
        style=curses.A_BOLD|curses.color_pair(2)
        if is_active:
            style=curses.A_REVERSE|curses.A_BOLD|curses.color_pair(2)
        self.pwin.addstr(0,0,"%-*s"%(self.wprime,self.title),style)
        #self.pwin.addstr(self.hprime-1,0,"+"+"-"*(self.wprime-2))# *self.w)
        #for j in range(self.h):
        #    self.pwin.addch(j+1,0,"|")
        #    self.pwin.addch(j+1,self.wprime-1,"|")
        self.pwin.refresh()
        
    def key(s):
        pass

class ButtonWidget(Widget):
    def __init__(self,parent,x,y,w,h,title,callback):
        Widget.__init__(self,parent,x,y,w,h,title)
        self.callback=callback
    def key(self,k):
        if k==ord('\n'):
            self.callback()

class EditWidget(Widget):
    def __init__(self,parent,x,y,w,h,title,default):
        Widget.__init__(self,parent,x,y,w,h,title)
        self.default=default
        self.editor=Textbox(self.win)
        self.win.addstr(0,0,default)
    def key(self,k):
        if k==ord('\n'):
            self.editor.edit()

class ListWidget(Widget):
    def __init__(self,parent,x,y,w,h,title,lst,callback,default):
        Widget.__init__(self,parent,x,y,w,h,title)
        self.lst=lst
        self.pad_width=reduce(max,map(len,self.lst))
        self.active=0
        for i in xrange(len(self.lst)):
            if self.lst[i]==default: self.active=i
        self.callback=callback

    def draw(self,is_active):
        Widget.draw(self,is_active)
        self.win.addstr(0,0,"Fujc")
        for i in range(self.h):
            style=curses.A_NORMAL
            base=0
            if len(self.lst)>self.h:
                base=self.active
            if base+i==self.active: # self.active:
                style=curses.A_REVERSE
            txt=""
            if base+i<len(self.lst):
                txt=self.lst[base+i]
            self.win.addstr(i,0,"%-*s"%(self.w-1,txt),style) # "%-*s"%(self.pad_width,self.lst[i]),style)
        self.win.refresh()

    def key(self,k):
        if k==curses.KEY_DOWN:
            self.active=(self.active+1)%len(self.lst)
        if k==curses.KEY_UP:
            self.active-=1
            if self.active<0: self.active=len(self.lst)-1
        self.callback(self.lst[self.active])


class GUI:
    def __init__(self):
        self.widgets=[]
        self.active=0
        self.quit=False


    def start(self):
        stdscr=curses.wrapper(self.run)
        
    def run(self,stdscr):

        self.initWidgets(stdscr)
        curses.init_pair(1,curses.COLOR_WHITE,curses.COLOR_RED)
        curses.init_pair(2,curses.COLOR_WHITE,curses.COLOR_BLUE)
        #stdscr.addstr(0,0,"%-80s"%self.title,curses.color_pair(1)|curses.A_BOLD)

        while 1:
            stdscr.addstr(0,0,"%-80s"%self.title,curses.color_pair(1)|curses.A_BOLD)
            for i in xrange(len(self.widgets)):
                self.widgets[i].draw(i==self.active)

            c=stdscr.getch()
            #stdscr.addstr(0,60,"%6s"%repr(c))
            #stdscr.addstr(16,0,'\t')
            stdscr.refresh()
            if c==ord('\t'):
                self.active=(self.active+1)%len(self.widgets)
            else:
                self.widgets[self.active].key(c)
            if self.quit:
                return
