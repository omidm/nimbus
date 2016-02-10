#!/usr/bin/python
######################################################################
# Copyright 2005, Andrew Selle.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#######################################################################
# Standard test wrapper script
# run with run_standard_tests.py <octree,3d> <example#> <resolution#>
#######################################################################
import os
import sys
import re
import time
import shutil
from pylab import *

# Define tests
test_names={1:"Flat surface",2:"Falling drop",3:"4 sources",4:"Sphere splashing",5:"Glass filling"}
binaries={"octree":"./fluids_octree","3d":"./fluids_3d"}

######################################################################
# class RUN
#  Data structure representing a single run's statistics
######################################################################
class RUN:
    def __init__(self,label,directory):
        self.label=label
        self.directory=directory
        # Frame and substep detection regexps
        self.substep_regexp=re.compile("END substep (\d+).+(\d+\.\d+) s")
        self.frame_regexp=re.compile("END Frame (\d+).+(\d+\.\d+) s")
        # series is all the data series hash of hashes... i.e. self.series = {"times": {1:3.4, 2:4.5}, "substeps": {1:10, 2:20} }
        # stats is a  hash of stat tuples (min, avg, max)
        self.series={}
        self.series_stats={}
        self.load()
        self.has_rendering=False

    ######################################################################
    # load series stats
    ######################################################################
    def series_stats_table(self):
        #return "foo"
        return [["Label","Min","Avg","Max"]],map(lambda x: (x[0],"%.2f"%x[1][0],"%.2f"%x[1][1],"%.2f"%x[1][2]),self.series_stats.items())


    ######################################################################
    # Generate stats for a esries
    ######################################################################
    def stats(self,x):
        values=x.values()
        return min(values),sum(values)/len(values),max(values)

    ######################################################################
    # parse output file
    ######################################################################
    def parse_output(self):
        # Series that we will populate
        self.series["Time"]={}
        self.series["Substep Count"]={}
        self.series["Substep Average Time"]={}
        # Open file
        fp=open(os.path.join(os.path.join(self.directory,"output.txt")))
        substep_count,substep_time_sum=0,0
        while 1:
            line=fp.readline()
            if line=="": break
            line=line.strip()
            frame_match=self.frame_regexp.match(line)
            substep_match=self.substep_regexp.match(line)
            if frame_match:
                try:
                    frame,frametime=int(frame_match.group(1)),float(frame_match.group(2))
                    if substep_count>0: average_substep_time=substep_time_sum/substep_count
                    self.series["Time"][frame]=frametime
                    self.series["Substep Count"][frame]=substep_count
                    self.series["Substep Average Time"][frame]=average_substep_time
                    substep_count,substep_time_sum=0,0 # reset substep counts
                except: pass
            elif substep_match:
                try:
                    substep_time=float(substep_match.group(2))
                    substep_count+=1
                    substep_time_sum+=substep_time
                except: pass

    ######################################################################
    # load run data
    ######################################################################
    def load(self):
        # Information
        self.info=eval(open(os.path.join(self.directory,"info.py")).read())
        # Frame memory
        self.series["Memory"]=eval(open(os.path.join(self.directory,"memory.py")).read())
        self.parse_output()
        # find min frame and max frame
        frames=map(lambda x: x.keys(),self.series.values())
        self.frame_range=(reduce(min,map(min,frames)),reduce(max,map(max,frames)))
        # compute stats for frame times and memories
        for s in self.series.keys():
            self.series_stats[s]=self.stats(self.series[s])

    ######################################################################
    # render
    ######################################################################
    def render(self):
        platform_modifier=""
        if os.environ.has_key('PLATFORM'): platform_modifier="_"+os.environ['PLATFORM']
            ray_tracing_binary=os.path.join(sys.environ["PHYSBAM"],"Projects/ray_tracing/ray_tracing%s"%platform_modifier)
                
        if not os.path.exists(os.path.join(self.directory,"render.scene")): return False
        try: os.makedirs(os.path.join(self.directory,"render"))
        except: pass

        # render undone frames
        for frame in range(self.frame_range[0],self.frame_range[1]+1):
            frame_filename=os.path.join(self.directory,"render.%05d.png"%frame)
            if not os.path.exists(frame_filename):
                print "Directory=%s Rendering Frame %d"%(self.directory,frame)
                render_cmd="cd %s;%s render.scene %d"%(self.directory,ray_tracing_binary,frame)
                #print "  %s"%render_cmd
                render_fp=os.popen(render_cmd)
                while 1:
                    line=render_fp.readline()
                    if line=="": break
        return True
        # make divx
        # fix on 64-bit
        # os.system("cd %s;mencoder mf://*.png -mf fps=30 -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=16000:vhq -o output.avi"%self.directory)

######################################################################
# Write the HTML for a table
# headers cells have the optional form (colspan, width, header_text)
######################################################################
def HTML_Table(headers,rows):
    text=""
    text+=("<table>\n")
    for h in headers:
        text+=("  <tr bgcolor=#7799ff>\n")
        for c in h:
            if type(c)==tuple:
                span,width,header=c
                text+=("    <th colspan=%d width=\"%s\">%s</th>\n"%(span,width,header))
            else:
                text+=("    <th>%s</th>\n"%c)
        text+=("  </tr>\n")
    color,alt_color="#eeeeee","#dddddd"
    for r in rows:
        color,alt_color=alt_color,color # alterante
        text+=("  <tr bgcolor=%s>\n"%color)
        for c in r:
            text+=("    <td>%s</td>\n"%c)
        text+=("  </tr>\n")
    
    text+=("</table>\n")
    return text

######################################################################
# Class REPORT_ENGINE Generate Comparisons
######################################################################
class REPORT_ENGINE:
    def Label_Extend(self,label,original_list):
        ret=[label]
        ret.extend(original_list)
        return ret

    def __init__(self,directory,directories):
        self.output_directory=directory
        self.directories=directories
        self.data={}

        # Open file for writing
        try: os.makedirs(output_directory)
        except: pass
        self.html=open(os.path.join(self.output_directory,"index.html"),"w")

        # Read Data And build data structure
        set=0 
        for directory in self.directories:
            set+=1
            self.data[directory]=RUN("Set %d"%set,directory)
        # Compute frame range
        self.frame_range=reduce(min,map(lambda x:self.data[x].frame_range[0],self.directories)),reduce(max,map(lambda x:self.data[x].frame_range[1],self.directories))

        # Generate Report Parts
        self.Basic_Info()
        self.Series_Data_Table()
        self.Series_Graphs()
        self.Thumbnails()


    ######################################################################
    # Function Basic_Info
    ######################################################################
    def Basic_Info(self):
        # Make basic info table
        self.html.write("<h1>Basic Data</h1>\n")
        self.html.write(HTML_Table([self.Label_Extend("Directory",self.directories)],
                              [self.Label_Extend("Label",map(lambda x: self.data[x].label,self.directories)),
                               self.Label_Extend("Info",map(lambda x: "<br>".join(self.data[x].info),self.directories)),
                               self.Label_Extend("Frame Range",map(lambda x: repr(self.data[x].frame_range),self.directories)),
                               self.Label_Extend("Series Stats",map(lambda x: HTML_Table(*self.data[x].series_stats_table()),self.directories))]))

    ######################################################################
    # Function Series_Data_Table
    ######################################################################
    def Series_Data_Table(self):
        # get other series'
        series=self.data[self.directories[0]].series.keys()
        directory_count,series_count=len(self.directories),len(series)
        self.html.write("<h1>Series</h1>\n")
        labels=[self.Label_Extend("Labels",[(directory_count,"%f%%"%(100/(series_count+.2)),s) for s in series]),
                self.Label_Extend("Frame",series_count*map(lambda x: self.data[x].label,self.directories))]
        rows=[]
        for frame in range(self.frame_range[0],self.frame_range[1]+1):
            row=[str(frame)]
            for s in series:
                for d in self.directories:
                    if self.data[d].series.has_key(s) and self.data[d].series[s].has_key(frame):
                        row.append("%.2f"%self.data[d].series[s][frame])
                    else:
                        row.append("");
            rows.append(row)
        self.html.write(HTML_Table(labels,rows))

    ######################################################################
    # function Series_Graphs - requires matplotlib
    ######################################################################
    def Series_Graphs(self):
        # make series graphs
        series=self.data[self.directories[0]].series.keys()
        figure_id=1
        for s in series:
            self.html.write("<h2>Series %s</h2>\n"%s)
            # collec the series for each example
            graph_data=[]
            for d in self.directories:
                if self.data[d].series.has_key(s):
                    graph_data.append(self.data[d].series[s].keys())
                    graph_data.append(self.data[d].series[s].values())
            # throw away old graph
            hold(False)
            lines=plot(*graph_data)
            # set series labels and enable legend
            for i in range(len(lines)): setp(lines[i],label=self.data[self.directories[i]].label)
            legend(loc='best')
            # axis labels
            xlabel("Frame")
            ylabel(s)
            title(s)
            # write to standard location
            savefig(os.path.join(self.output_directory,"%d.png"%figure_id))
            # add html link and increment count'
            self.html.write("<p><img src=\"%d.png\"></P>"%figure_id)
            figure_id+=1

    ######################################################################
    # thumbnails
    ######################################################################
    def Thumbnails(self):
        self.html.write("<H1>Thumbnails</h1>")
        # get all data
        labels=self.Label_Extend("Frame",map(lambda x:self.data[x].label,self.directories))
        rows=[]
        for frame in range(self.frame_range[0],self.frame_range[1]+1):
            row=[]
            row.append("Frame %d"%frame)
            for d in self.directories:
                # choose thumb directory
                thumbdir=os.path.join(self.output_directory,self.data[d].label)
                try: os.makedirs(thumbdir)
                except: pass
                self.data[d].render()
                filename="render.%05d.png"%frame
                thumbfilename="thumb.%05d.png"%frame
                srcimg=os.path.join(self.data[d].directory,filename)
                if frame%10==0:
                    if os.path.exists(srcimg):
                        shutil.copy(srcimg,os.path.join(thumbdir,filename))
                        os.system("convert -geometry 160x120 \"%s\" \"%s\""%(srcimg,os.path.join(thumbdir,thumbfilename)))
                        row.append("<a href=\"%s\"><img src=\"%s\"></a>\n"%(os.path.join(self.data[d].label,filename),os.path.join(self.data[d].label,thumbfilename)))
                    else:
                        row.append("N/A")
                    rows.append(row)

        self.html.write(HTML_Table([labels],rows))
        # copy data 

# Generate report
REPORT_ENGINE("html",sys.argv[1:])
