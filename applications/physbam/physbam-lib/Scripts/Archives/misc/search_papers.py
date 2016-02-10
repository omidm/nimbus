#!/usr/bin/python
import re
import sys
import pygtk
import gtk
import urllib
import os

RED=chr(27) + '[1;31m'
#BRIGHTRED=chr(27) + '[1;31m'
GREEN=chr(27) + '[1;32m'
BLUE=chr(27) + '[1;34m'
CLEAR=chr(27) + '[00m'

name_re=re.compile("^(.+)-\s*([A-Za-z]+\s?[A-Za-z]*)\.?\s*$")
journal_re=re.compile("^([A-Za-z]+).*Volume\s([0-9]+).*Issue\s([0-9]+).*$")

if len(sys.argv)==2:
   name=sys.argv[1]
else:
   name=None
print name

counts={}
all_papers=[]
journal=""
volume=""
issue=""

for line in open(os.path.join(os.environ["PHYSBAM"],"Public_Library","Documentation","PAPERS_TO_PRINT.txt")).readlines():
    if line=="": break
    line=line[:-1]

    m=name_re.match(line)
    j=journal_re.match(line)
    
    if m:
        if name != None:
           if m.group(2)==name:
              print "%20s"%(m.group(1))
        else:
            print "[%s%15s%s] %s"%(BLUE,m.group(2),CLEAR,m.group(1))
        if  not counts.has_key(m.group(2)):
            counts[m.group(2)]=0
        counts[m.group(2)]+=1
        all_papers.append((m.group(1).strip(),m.group(2),journal,volume,issue))
            
    elif j:
        journal=j.group(1)
        volume=j.group(2)
        issue=j.group(3)
        print "%s%s %s%s %s%s%s"%(RED,journal,GREEN,volume,BLUE,issue,CLEAR)
    else:
        print "%s%s%s"%(GREEN,line,CLEAR)

#print all_papers

class GUI:
   def delete_event(self,widget,event,data=None):
      gtk.main_quit()
      return False

   def click_event(self,widget,event,data=None):
      selection=widget.get_selection()
      pathlist=[]
      def foreach_cb(model,path,iter,pathlist):
         pathlist.append(model.get(iter,0,1,2,3,4))
      selection.selected_foreach(foreach_cb,pathlist)
      if len(pathlist)>0:
         title="\""+pathlist[0][0]+"\""
         if pathlist[0][2]=='JCP':
            cmd="firefox http://www.sciencedirect.com/science?_ob=PublicationURL\&_cdi=6863\&_auth=y\&_acct=C000012078\&_version=1\&_urlVersion=0\&_userid=145269\&_pubType=J\&md5=cdcd6e727bb4b0c150afd6ca922b5e89 &"
            print "%s%s %s%s %s%s%s: %s"%(RED,pathlist[0][2],GREEN,pathlist[0][3],BLUE,pathlist[0][4],CLEAR,title)
         else:
            print title
            cmd="firefox http://www.google.com/search?%s &"%urllib.urlencode({"q":title})
         os.system(cmd)


   
   def __init__(self,papers):
      self.papers=papers
      self.window=gtk.Window(gtk.WINDOW_TOPLEVEL)
      self.window.set_title("Search for Papers")

      self.window.connect("delete_event", self.delete_event)
      
      self.liststore=gtk.ListStore(str,str,str,str,str)

      for paper in self.papers:
         self.liststore.append(paper)

      self.tvcolumn_title=gtk.TreeViewColumn("Title")
      self.tvcolumn_name=gtk.TreeViewColumn("name")
      self.listview=gtk.TreeView(self.liststore)
      self.listview.append_column(self.tvcolumn_title)
      self.listview.append_column(self.tvcolumn_name)
      self.cell1=gtk.CellRendererText()
      self.cell2=gtk.CellRendererText()
      self.tvcolumn_title.pack_start(self.cell1,True)
      self.tvcolumn_title.add_attribute(self.cell1,'text',0)
      self.tvcolumn_name.pack_start(self.cell2,True)
      self.tvcolumn_name.add_attribute(self.cell2,'text',1)
      self.listview.connect("row-activated", self.click_event)

      self.window.add(self.listview)

      self.window.show_all()
      

guys=counts.items()
guys.sort(lambda x,y: y[1]-x[1])
print "%-15s %-5s"%("Name","Count")
print "%-15s %-5s"%("-"*15,"-"*5)
print "\n".join(map(lambda x: "%-15s %5d"%x,guys))

#g=GUI(["BC","C"])
g=GUI(all_papers)
gtk.main()
