#!/usr/bin/python

import sys
import os 
import re

class DAG:
    def __init__(self,name):
        self.name=name
        self.children=[]
        global display_hash
        display_hash=dict()
        global cyclic_check
        cyclic_check=[]
        #self.display_exclusion_regex=re.compile(".*(LOG|READ_WRITE|ARRAY|VECTOR|UTILITIES|MATRIX)")

    def display(self,offset=0):
        global display_hash
        global cyclic_check
        if(self.name in display_hash):
            if(display_hash[self.name]):
		return
        if(self.name in cyclic_check):
            print cyclic_check.keys()
            return

        #match=self.display_exclusion_regex.match(self.name)
        #if(match): return
        for i in range(0,offset): print " ",
        print os.path.basename(self.name)
        display_hash[self.name]=True
        cyclic_check.append(self.name)
        for child in self.children:
            child.display(offset+1)
        cyclic_check.pop() # check that this gives proper behavior
        #display_hash[self.name]=False

    def search_for_cyclic_dependencies(self):
        global cyclic_hash
        if(self.name in cyclic_hash):
            if(cyclic_hash[self.name]):
                print self.name,", "
                return True

        cyclic_hash[self.name]=True
        for child in self.children:
            if(self.search_for_cyclic_dependencies()):
                print self.name,", "
                return True
        cyclic_hash[self.name]=False

class DEPENDENCIES_DAG:
    def __init__(self):
        self.library_path=os.environ["PHYSBAM"]+"/Public_Library/"
        print self.library_path
        self.library_include_regex=re.compile("^#include *<([^>]+)>$")
        self.local_include_regex=re.compile("^#include *\"([^\"]+)\"$")
        self.header_regex=re.compile("^(.*)\.h$")
        #self.exclusion_regex=re.compile(".*((OC|QUAD)TREE|RLE|DYADIC|VOF_ADVECTION|MPI)")
        self.exclusion_regex=re.compile("asdfasdfasdfafdasdfasdf")
        self.if_regex=re.compile("^#if (.*)$")
        self.ifdef_regex=re.compile("^#ifdef (.*)$")
        self.ifndef_regex=re.compile("^#ifndef (.*)$")
        self.endif_regex=re.compile("^#endif")
        self.files_to_parse=[]
        self.dag_entry_for_file=dict()
        self.defined_macros=dict()
        self.is_parsable_scope=[]
        self.is_parsable=True

    def Add_File_To_Queue(self,file,node,is_cpp):
#        if((not is_cpp) & (os.path.exists(file))):
        if(os.path.exists(file)):
            if(file not in self.dag_entry_for_file):
                child=DAG(file)
            else: child=self.dag_entry_for_file[file]
            child.parent=node
            if(child not in node.children): node.children.append(child)

        if((file not in self.dag_entry_for_file) & (os.path.exists(file))):
            self.files_to_parse.append(file)
            self.dag_entry_for_file[file]=child

    def Add_Include_Files_To_Queue(self,file,node):
        #print file
        local_path=os.path.dirname(file)+"/"
        try:
            file_handle=open(file)
        except:
            print "unable to open file %s"%file
            return

        for line in file_handle:
            self.Parse_Macros(line,file)
            if(not self.is_parsable): continue

            include_file=""
            match=self.library_include_regex.match(line)
            if(match): include_file=self.library_path+match.group(1)
            else:
                match=self.local_include_regex.match(line)
                if(match): include_file=local_path+match.group(1)

            if(match):
                #if(self.exclusion_regex.match(include_file)):
                #    print "EXCLUDING %s"%include_file
                #    continue

                self.Add_File_To_Queue(include_file,node,False)
                if(os.path.exists(include_file)):
                    #if this is a header file, add the corresponding cpp file
                    match=self.header_regex.match(include_file)
                    if(match):
                        cpp_file=match.group(1)+".cpp"
                        if(cpp_file!=file):
                            self.Add_File_To_Queue(cpp_file,node,True)

        file_handle.close()

    def Parse_Macros(self,line,file):
        macro_match=self.ifdef_regex.match(line)
        if(macro_match):
            if(macro_match.group(1) in self.defined_macros): self.is_parsable_scope.append(self.is_parsable)
            else:
                #print "Found undefined ifdef macro %s"%macro_match.group(1)
                self.is_parsable_scope.append(self.is_parsable)
                self.is_parsable=False
        else:
            macro_match=self.ifndef_regex.match(line)
            if(macro_match):
                if(macro_match.group(1) in self.defined_macros):
                    #print "Found undefined ifndef macro %s"%macro_match.group(1)
                    self.is_parsable_scope.append(self.is_parsable)
                    self.is_parsable=False
                else: self.is_parsable_scope.append(self.is_parsable)
            else:
                macro_match=self.if_regex.match(line) # assume all #if statements are true...
                if(macro_match):
                    self.is_parsable_scope.append(self.is_parsable)
                else:
                    macro_match=self.endif_regex.match(line)
                    if(macro_match):
                        #print file
                        self.is_parsable=self.is_parsable_scope.pop()
                        #print "Set is_parsable to ",self.is_parsable

       

    def Construct_Dag(self,root_file):
        root_node=DAG(root_file)
        self.files_to_parse.append(root_file)
        self.dag_entry_for_file[root_file]=root_node
        while(len(self.files_to_parse)):
            current_file=self.files_to_parse.pop()
            current_node=self.dag_entry_for_file[current_file]
            self.Add_Include_Files_To_Queue(current_file,current_node)
        print "number of depending files=%d"%len(self.dag_entry_for_file.keys())
        #print self.dag_entry_for_file.keys()
        root_node.display()
        #root_node.search_for_cyclic_dependencies()

def main():
    dependencies_dag=DEPENDENCIES_DAG()
    dependencies_dag.Construct_Dag(sys.argv[1])

if __name__=="__main__":
    main()
