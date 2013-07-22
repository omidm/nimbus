#!/usr/bin/python
######################################################################
# Copyright 2006, Andrew Selle
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# log_parse.py - parse the log file
######################################################################
import xml.parsers.expat
import pickle
import re

print_prints=False
print_stats=False
print_scopes=False

frame_re=re.compile("Frame ([0-9]+)")
# colors
RED=chr(27) + '[1;31m'
#BRIGHTRED=chr(27) + '[1;31m'
GREEN=chr(27) + '[1;32m'
BLUE=chr(27) + '[1;34m'
CLEAR=chr(27) + '[00m'

# scope information
root=None
class SCOPE:
    def __init__(self,id,name=None):
        self.id=id
        self.name=name
        self.count=0
        self.time=0
        self.children={}
        self.children_order=[]
        self.stats={}

    def show(self,depth=0):
        def robust_divide(x,y):
            try:
                return x/y
            except:
                return 0
        print "%s%*s%-*s %s%10.04f s%s (%5d time%s %10.04f s/time"%(GREEN,2*depth,"",80-2*depth,self.id,RED,self.time,CLEAR,self.count,(') ','s)')[self.count!=1],robust_divide(self.time,self.count))
        for id in self.children_order:
            self.children[id].show(depth+1)
        for id,stats in self.stats.items():
            if type(stats)==float:
                print "%*s%s = %f"%(2*depth+1,"",id,stats)
            elif type(stats)==int:
                print "%*s%s = %d"%(2*depth+1,"",id,stats)

    def find(self,label):
        if self.id==label: return self
        for id,child in self.children.items():
            value=child.find(label)
            if value: return value
        return None
    
    def sum_stat(self,label):
        sum=0
        if self.stats.has_key(label): sum+=self.stats[label]
        for id,child in self.children.items():
            sum+=child.sum_stat(label)
        return sum
    
    def find_time(self,label):
        if self.id==label:
            return self.time
        else:
            total_time=0.
            for id,child in self.children.items():
                total_time+=child.find_time(label)
            return total_time

def find_scope(scope_id,parent):
    global root
    if not parent:
        root=SCOPE(scope_id)
        return root
    else:
        if not parent.children.has_key(scope_id):
            parent.children[scope_id]=SCOPE(scope_id)
            parent.children_order.append(scope_id)
        return parent.children[scope_id]

def float_or_int(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

scope_names,scopes,one_line,timing,printing=[],[],True,0,False
def start_element(name,attrs):
    global timing,one_line,printing,enabled
    def handle_newline():
        global one_line
        if one_line:
            one_line=False
            print " "
    if name=="scope":
        scope_name=attrs['name']
        scope_id=attrs.get('id',scope_name)
        scope_names.append(scope_name)
        parent=None
        if len(scopes)>0: parent=scopes[-1]
        if scope_id=='FRAME':
            match=re.match('^Frame (\d+)$',scope_name)
            if int(match.group(1))>=start_frame:
                enabled=True
        scopes.append(find_scope(scope_id,parent))
        if print_scopes:
            handle_newline()
            one_line=True
            indent=len(scope_names)
            print ("%s%*s%-*s %s"%(GREEN,2*indent,"",80-2*indent,attrs["name"],CLEAR)),
    elif name=="time":
        indent=len(scope_names)
        timing=attrs["value"]
    elif name=="print" and print_prints:
        handle_newline()
        printing=True
    elif name=="stat":
        if not scopes[-1].stats.has_key(attrs["name"]):
            scopes[-1].stats[attrs["name"]]=0
        scopes[-1].stats[attrs["name"]]+=float_or_int(attrs["value"])
        if print_stats:
            handle_newline()
            indent=len(scope_names)
            print "%*s%s%s = %s%s"%(2*indent+2,"",BLUE,attrs["name"],attrs["value"],CLEAR)
        
def char_data(data):
    if printing:
        print "%*s%s"%(2*len(scope_names)+2,"",data)
    
class ABORT(Exception):
    pass

def end_element(name):
    global timing,one_line,printing
    if name=="scope":
        if print_scopes:
            if one_line:
                print "%s%s s%s"%(RED,timing,CLEAR)
            else:
                print "%s%*s%-*s  %s%s s%s"%(GREEN,2*len(scope_names),"",80-2*len(scope_names),"END "+scope_names[-1],RED,timing,CLEAR)
            one_line=False
        if enabled:
            scopes[-1].time+=float(timing)
            scopes[-1].count+=1
        match=scope_names[-1].startswith('Frame') # see if we are a frame
        poped_scope_name=scope_names.pop()
        scopes.pop()
        # if we have gotten to end frame abort
        if match:
            m=frame_re.match(poped_scope_name)
            if int(m.group(1))>=end_frame:
                raise ABORT
    elif name=="print":
        printing=False
    
def Parse(filename,start_frame_input=None,end_frame_input=100000000):
    global start_frame
    global end_frame
    global enabled
    if end_frame_input is None: end_frame=100000000
    else: end_frame=end_frame_input
    start_frame=start_frame_input
    enabled=start_frame is None
    p=xml.parsers.expat.ParserCreate()
    p.StartElementHandler=start_element
    p.EndElementHandler=end_element
    p.CharacterDataHandler = char_data
    fp=open(filename)
    try:
        p.ParseFile(fp)
    except xml.parsers.expat.ExpatError,e:
        print "XML Parse Error on line %d column %d -- code %d"%(e.lineno,e.offset,e.code)
    except ABORT:
        print "Finished with %d frames"%end_frame
        
    return root
