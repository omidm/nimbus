#!/usr/bin/python
import os
import sys
import time
import curses
import curses.wrapper
import shutil

from GUI import ButtonWidget
from GUI import ListWidget
from GUI import EditWidget
from GUI import Widget
from GUI import GUI

def choose_list(name,items):
    i=1
    print "%s"%name
    print "-"*len(name)
    for item in items:
        print "%3d %s"%(i,repr(item))
        i+=1
    print "Choose: ",
    sys.stdin.readline()


models="/solver/vol3/hair1/data"

#choose_list("Model",os.listdir(models))



class SimGUI(GUI):
    def __init__(self):
        GUI.__init__(self)
        self.modeldir="/solver/vol3/hair1/data"
        self.widgets=[]
        self.model=None
        self.sim=None
        self.ok=False

        self.outputName=None
        self.outputBinary=None
        self.outputModel=None
        self.outputSim=None

    def initWidgets(self,stdscr):
        self.modelWidget=ListWidget(stdscr,0,1,24,8,"Model",os.listdir("/solver/vol3/hair1/data"),self.updateSimWidget,"body2")
        self.simWidget=ListWidget(stdscr,25,1,50,8,"Sim",["none"],self.updateSim,"")
        self.nameWidget=EditWidget(stdscr,0,9,80,3,"Name","untitled")
        self.binaryWidget=EditWidget(stdscr,0,13,80,3,"Binary","$PHYSBAM/Projects/solids_3d/solids_3d_nocona")
        self.processorWidget=EditWidget(stdscr,0,17,10,3,"Procs","1")
        self.exampleWidget=EditWidget(stdscr,12,17,10,3,"Example","1")
        self.updateSimWidget(self.modelWidget.lst[self.modelWidget.active])
        self.title="Create a Hair Sim"

        self.widgets=[self.modelWidget,self.simWidget,self.nameWidget,self.binaryWidget,self.processorWidget,self.exampleWidget,
                      ButtonWidget(stdscr,50,20,10,2,"Cancel",self.cancelCallback),
                      ButtonWidget(stdscr,60,20,10,2,"Ok",self.okCallback)]

    def editCallback(self):
        os.system("vi")
        #self.stdscr.refresh()
        

    def cancelCallback(self):
        self.quit=True
    def okCallback(self):
        self.outputName=self.nameWidget.editor.gather()
        self.outputBinary=self.binaryWidget.editor.gather()
        self.outputModel=self.model
        self.outputSim=self.sim
        self.outputProcessors=int(self.processorWidget.editor.gather())
        self.outputExample=int(self.exampleWidget.editor.gather())
        
        self.quit=True
        self.ok=True

    def updateSim(self,param):
        if self.sim!=param:
            self.sim=param

    def updateSimWidget(self,param):
        if self.model!=param:
            self.model=param
            self.simWidget.lst=os.listdir(os.path.join(self.modeldir,self.model))

            self.simWidget.active=0



gui=SimGUI()
gui.start()

if gui.ok:
    datestr="%04d%02d%02d-%02d%02d"%(time.localtime()[:5])
    rundir=os.path.join("/solver/vol3/hair1/sims/%s-%s"%(datestr,gui.outputName))

    data_directory=os.path.join(gui.modeldir,gui.model)
    model=gui.model
    sim=gui.sim
    paramfile=os.path.join(data_directory,sim,"parameters")
    os.system("vi %s"%paramfile)

    if os.path.exists(rundir):
        print "Path exists already try new sim name"
    else:
        os.mkdir(rundir)
        os.system("%s --copy %s"%(gui.outputBinary,rundir))
        shutil.copy("/data/aselle/PhysBAM/Scripts/linux_sim_scripts/mpi_cluster_sim.py",rundir)
        os.chdir(rundir)
        cmd=("python mpi_cluster_sim.py %d %s ./solids_3d -example HAIR_SIM_TESTS -d %s -modelname %s -hairsim %s %d "%(gui.outputProcessors,gui.outputName,data_directory,model,sim,gui.outputExample))
        print cmd
        os.system(cmd)


