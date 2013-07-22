#!/usr/bin/python
import os
import sys
from elementtree.ElementTree import parse

if len(sys.argv)!=2:
    print>>sys.stderr, "Usage: formatlog <dir|logfile>"
    sys.exit(1)
log=sys.argv[1]
if os.path.isdir(log): log=os.path.join(log,'log.txt')
tree=parse(open(log))
root=tree.getroot()

RED=chr(27) + '[1;31m'
#BRIGHTRED=chr(27) + '[1;31m'
GREEN=chr(27) + '[1;32m'
BLUE=chr(27) + '[1;34m'
CLEAR=chr(27) + '[00m'

def display(root,indent):
    print "%s%*s%-*s %s s%s"%(GREEN,2*indent,"",80-2*indent,root.attrib["name"],root[-1].attrib["value"],CLEAR)
    #print "%*s%s"%(5,"","hiu")
    #if len(root)==1: print " %s"%(root[-1].attrib["value"])
    #print " "
    for child in root:
        if child.tag=="time":
            pass
        elif child.tag=="stat":
            print "%*s%s%s = %s%s"%(2*indent+2,"",BLUE,child.attrib["name"],child.attrib["value"],CLEAR)
            pass
        elif child.tag=="print":
            print "%*s%s%s%s"%(2*indent+2,"",RED,child.text,CLEAR)
            pass
        elif child.tag=="error":
            print child.text
            pass
        else:
            display(child,indent+1)

display(root,0)
