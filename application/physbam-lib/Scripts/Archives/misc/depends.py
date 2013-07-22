#!/usr/bin/python
import os
import re
import sys
from optparse import OptionParser

include_re=re.compile("^#\s*include\s+<(.+)>")

usage="usage: %prog [options] <log fileanme>"
parser=OptionParser(usage)
parser.add_option("-g","--generate",action="store_true",dest="generate",default=False)
parser.add_option("-t","--tree",action="store_true",dest="tree",default=False)
parser.add_option("-f","--file-dependency",action="store",type="string",dest="file_depends",default=None)
parser.add_option("-d","--dir-dependency",action="store",type="string",dest="dir_depends",default=None)
(options,args)=parser.parse_args()


#def test(x):
#    m=include_re.match(x)
#    match=False
#    if m: match=True
#    print "%s -> %s"%(x,match)
#test("#include <stdio.h>")
#test("#    include <stdio.h> // lame")
#test(" #    include <stdio.h>")

class DependencyWalker:
    def dir_func(self,arg,dir,files):
        def cpp_or_h(x): return x.endswith(".cpp") or x.endswith(".h")
        self.files.extend(filter(cpp_or_h,map(lambda x: os.path.join(dir,x).replace("./",""),files)))
        
    def __init__(self):
        self.files=[]
        os.path.walk(".",self.dir_func,None)
        self.parents={} # maps from filename to dependencies

    def get_all_deps(self):
        todo_list=[];todo_list.extend(self.files)
        while len(todo_list)>0:
            guy=todo_list.pop()
            todo_list.extend(filter(lambda x: not self.parents.has_key(x) and os.path.exists(x),self.get_deps(guy)))

        self.find_dir_depends()

    def get_deps(self,file):
            included_files=[]
    
            #print "File: %s"%file
            try:
                fp=open(file)
            except:
                self.parents[file]=[]
                return []
            while 1:
                line=fp.readline()
                if line=="": break
    
                m=include_re.match(line)
                if m:
                    included_files.append(m.group(1))
            self.parents[file]=included_files
            return included_files

    def find_roots(self):
        for k,deps in self.parents.items():
            if len(deps)==0: print k

    def file_to_dir(self,x):
        items=x.split("/")
        if len(items)!=2: return None
        return items[0]

    def find_dir_depends(self):
        self.dirs={}
        for k,deps in self.parents.items():
            dir=self.file_to_dir(k)
            if not self.dirs.has_key(dir):
                self.dirs[dir]={}
            for i in deps:
                dep_dir=self.file_to_dir(i)
                if dep_dir and dep_dir != dir:
                    self.dirs[dir][dep_dir]=True

    def make_dir_graph(self):
        fp=open("foo.dot","w")
        fp.write("digraph {\n")
        for dir,dep_dirs in self.dirs.items():
            for dep_dir in dep_dirs:
                fp.write("  \"%s\" -> \"%s\"\n"%(dep_dir,dir))
        fp.write("}\n")

    def find_root(self):
        for dir,dep_dirs in self.dirs.items():
            print "Diretory %s has %d dependencies"%(dir,len(dep_dirs))
            print "\n".join(map(lambda x:"    "+x,dep_dirs))

    def Save(self,filename):
        open(filename,"w").write(repr((self.dirs,self.files,self.parents)))
        
    def Load(self,filename):
        self.dirs,self.files,self.parents=eval(open(filename).read())
        

os.chdir(os.path.join(os.environ["PHYSBAM"],"Public_Library"))
dep=DependencyWalker()
if options.generate:
    dep.get_all_deps()
    dep.Save("depends.cache")
else:
    dep.Load("depends.cache")
    #print dep.dirs
#dep.find_roots()
#dep.find_dir_depends()
#dep.make_dir_graph()

if options.file_depends:
    if not options.tree:
        print "\n".join(dep.parents[options.file_depends])
    else:
        def recurse_depend(x,seen,level):
            if not dep.parents.has_key(x) or seen.has_key(x): return
            seen[x]=True
            for i in dep.parents[x]:
                print  "%*s %s"%(level," ",i)
                recurse_depend(i,seen,level+1)
        recurse_depend(options.file_depends,{},0)

if options.dir_depends:
    if not options.tree:
        print "\n".join(dep.dirs[options.dir_depends].keys())
    else:
        def recurse_depend(x,seen,level):
            if not dep.dirs.has_key(x) or seen.has_key(x): return
            seen[x]=True
            for i in dep.dirs[x].keys():
                print  "%*s %s"%(level," ",i)
                recurse_depend(i,seen,level+1)
        recurse_depend(options.dir_depends,{},0)

        
        
    
#dep.find_root()


                                    

            
