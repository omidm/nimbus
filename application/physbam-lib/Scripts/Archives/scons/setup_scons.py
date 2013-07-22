#!/usr/bin/python
import os
import shutil
import sys

def relative_path(base,directory):
    b,d=os.path.abspath(base).split(os.sep),os.path.abspath(directory).split(os.sep)
    i=0
    while i<len(b) and i<len(d) and b[i]==d[i]: i+=1
    path=map(lambda x:'..',b[i:])+d[i:]
    if len(path)==0: return "."
    return os.path.join(*path)

def symlink_wrapper(x,y):
    if os.path.islink(y): os.unlink(y)
    try:
        os.symlink(x,y)
    except OSError:
        pass
if not os.__dict__.has_key("symlink"): operation=shutil.copy
else: operation=symlink_wrapper

scripts_directory=os.path.abspath(os.path.dirname(sys.argv[0]))
physbam_directory=os.path.abspath(os.path.join(scripts_directory,'..','..','..'))
assert(scripts_directory==os.path.join(physbam_directory,'Scripts','Archives','scons'))

for root,dirs,files in os.walk(scripts_directory):
    destdir=os.path.join(physbam_directory,relative_path(scripts_directory,root))
    if os.path.exists(destdir):
        for file in files:
            if not file.startswith("SC") or file=="README": continue
            src,dest=os.path.join(root,file),os.path.join(destdir,file)
            src,dest=relative_path(destdir,os.path.normpath(src)),os.path.normpath(dest)
            print "%s -> %s"%(src,dest)
            os.chdir(destdir)
            operation(src,dest)

