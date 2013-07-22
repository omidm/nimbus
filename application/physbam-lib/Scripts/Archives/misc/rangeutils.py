#!/usr/bin/python
import os
import re

RED=chr(27) + '[1;31m'
#BRIGHTRED=chr(27) + '[1;31m'
GREEN=chr(27) + '[1;32m'
BLUE=chr(27) + '[1;34m'
CLEAR=chr(27) + '[00m'


def format_sub_range(items):
    last,start=None,None
    subranges=[]
    for i in items:
        if last != i-1:
            if start!=None: subranges.append((start,last))
            start=i
        last=i
    if start!=None: subranges.append((start,last))
    #print subranges
    #print map(format_range_tuple,subranges)
    return ",".join(map(lambda x:format_range_tuple(False,x),subranges))

def format_range_tuple(use_wildcards,range_tuple):
    #print "range tuple %s"%repr(range_tuple)
    if range_tuple:
        start,end=range_tuple
        if start==end: return str(start)
        elif use_wildcards: return str("*")
        else: return "%d-%d"%(start,end)
    else:
        return ""

class FileRange:
    def __init__(self,key,num):
        self.items={}
        self.missing=[]
        self.pattern=key
        self.padding=0
        self.number_format="%d"
        #print "creating %s num %s"%(repr(key),repr(num))
        if num != None: self.extrema=(num,num)
        else: self.extrema=None
        self.use_wildcards=False
        
    def append(self,num,padding=0):
        self.items[num]=True
        if num!=None:
            currmin,currmax=self.extrema
            self.extrema=min(currmin,num),max(currmax,num)
        if padding:
            self.padding=max(padding,self.padding)

    def compute_padding(self):
        if self.padding>0:
            self.number_format="%%0%dd"%self.padding

    def compute_missing(self):
        self.missing=[]
        if self.extrema:
            start,end=self.extrema
            for i in range(start,end+1):
                if not self.items.has_key(i): self.missing.append(i)
    
    def format(self,use_wildcards,show_missing):
        missing_str=""
        if show_missing:
            self.compute_missing();
            if len(self.missing)>0:
                missing_str="\t\t\t\t\t%s(MISSING %s)%s"%(RED,format_sub_range(self.missing),CLEAR)
        return self.pattern[0]+format_range_tuple(use_wildcards,self.extrema)+self.pattern[1]+missing_str

    def filename(self,num):
        if num==None:
            return self.pattern[0]+self.pattern[1]
        else:
            return self.pattern[0]+(self.number_format%num)+self.pattern[1]

    def filenames(self):
        return map(self.filename,self.items.keys())

            
class FileRanges:
    def __init__(self,basedir):
        # all ranges
        self.ranges={}
        # get sorted list of files 
        names=os.listdir(basedir)
        names.sort()
    
        # cut based on range
        number_regexp=re.compile("^-?[0-9]+$")
        def newcut(x):
            blocks=x.split(".")
            map(lambda x: number_regexp.match(x) != None,blocks)
            padded=0
            for i in xrange(len(blocks)-1,-1,-1):
                if number_regexp.match(blocks[i])!=None:
                    if i==len(blocks)-1:
                        if blocks[i][0]=="0": padded=len(blocks[i])
                        return ".".join(blocks[:i])+".",int(blocks[i]),"",padded
                    else:
                        if blocks[i][0]=="0": padded=len(blocks[i])
                        return ".".join(blocks[:i])+".",int(blocks[i]),"."+".".join(blocks[i+1:]),padded
            return x,None,"",padded
        cut_names=map(newcut,names)
    
        # extract items
        current_prefix=None
        for i in cut_names:
            prefix,num,suffix,padded=i
            key=(prefix,suffix)
            #print i
            # make new entry for prefix
            if not self.ranges.has_key(key):
                self.ranges[key]=FileRange(key,num)
            self.ranges[key].append(num,padded)
        for i in self.ranges.values(): i.compute_padding()

    def format(self,use_wildcards=False,show_missing=False):
        return "\n".join(map(lambda x: x.format(use_wildcards,show_missing),self.ranges.values()))
            
    
