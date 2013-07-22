#!/usr/bin/python
import sys
import os
import time
import sys
import os
import re
import Image
import ImageDraw
import traceback

class VIEWER:
    def __init__(self,name,viewer,keys,camera_script):
        self.name,self.viewer,self.keys,self.camera_script=name,viewer,keys,camera_script
    def generate(self,result_directory,name,output_directory):
        output_movie_file=os.path.join(result_directory,name+"_"+self.name+".mov")
        strip_file=output_movie_file.replace(".mov",".png")
        view_command="%s -so %s -keys %s -camera_script %s -offscreen %s >& /dev/null"%(self.viewer,output_movie_file,self.keys,self.camera_script,output_directory)
        print view_command
        os.system(view_command)
        self.sample_filmstrip(output_movie_file,strip_file,4)
        
        return output_movie_file,strip_file
        
    def sample_filmstrip(self,movie,image,samples):
        r=re.compile("(?:.|\n)*Duration: ([0-9]+):([0-9]+):([0-9]+).([0-9])+")
        cmd="ffmpeg -i %s"%movie
        stdin,stderr=os.popen4(cmd)
        txt=stderr.read()
        m=r.match(txt)
        if m:
            h,m,s=int(m.group(1)),int(m.group(2)),float(m.group(3)+"."+m.group(4))
            sec=(h*60.+m)*60.+s
            print "%d h %d m %f s (%f seconds)"%(h,m,s,sec)
        else:
            return
        dt=sec/(samples)
        thums=[]
        resx,resy=160,120
        for i in range(samples):
            time=dt*i
            thumb_name="/tmp/img-%d.jpg"%i
            thumb_cmd="ffmpeg -y -i %s -ss %f -s %dx%d -vcodec png -vframes 1 -an -f rawvideo  %s"%(movie,time,resx,resy,thumb_name)
    
            print thumb_cmd
            os.system(thumb_cmd)
            
            thums.append(thumb_name)

        try:
            output=Image.new('RGB',(len(thums)*resx,resy),(30,30,70))
            for i in range(len(thums)):
                thum=thums[i]
                output.paste(Image.open(thum),(i*resx,0))
                os.unlink(thum)
            output.save(image)
        except:
            traceback.print_exc(file=sys.stdout)

class TEST:
    def __init__(self,result_directory,name,directory,command,output_directory,views):
        self.result_directory=result_directory
        self.name,self.directory,self.command,self.output_directory=name,directory,command,output_directory
        self.views=views
        self.video_files=[]

    def run(self):
        oldcwd=os.getcwd()
        os.chdir(os.path.join(oldcwd,self.directory))
        start_time=time.time()
        log_output=os.path.join(self.result_directory,self.name+".log")
        self.command+=" >& "+log_output
        print self.command
        self.exitcode=os.system(self.command)
        end_time=time.time()
        self.duration=end_time-start_time
        open(os.path.join(self.result_directory,self.name+".out"),"w").write("%s seconds=%d result=%d\n"%(self.name,end_time-start_time,self.exitcode))

        # run viewer
        for view in self.views:
            files=view.generate(self.result_directory,self.name,self.output_directory)
            self.video_files.append(files)

        os.chdir(oldcwd)
        
        return self.exitcode

    def report(self,fp):
        self.logfile=self.name+".log"
        fp.write("<tr><td class=entry>%s</td><td class=entry>%d s</td><td class=entry>%d</td><td class=entry><a href=\"%s\">%s</a></td></tr>\n"%(self.name,self.duration,self.exitcode,self.logfile,self.logfile))
        for video,thumbnail in self.video_files:
            fp.write("<tr><td colspan=4><a href=\"%s\"><img src=\"%s\"><br>"%(os.path.basename(video),os.path.basename(thumbnail)))
            fp.write("%s</td></tr>\n"%(os.path.basename(video)))

class VISUAL_TESTS:
    def __init__(self,name,result_directory):
        self.name=name
        self.result_directory=result_directory
        if not os.path.exists(result_directory) and not os.path.isdir(result_directory):
            os.mkdir(result_directory)
        self.tests=[]
        arch=os.popen("uname -m").read().strip()
        if arch=="x86_64": self.platform="nocona"
        else: self.platform="pentium4"

    def run(self):
        open(os.path.join(self.result_directory,"style.css"),"w").write(open(os.path.join("Scripts/test/style.css")).read())
        
        fp=open(os.path.join(self.result_directory,"index.html"),"w")
        fp.write("""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
        <html>
        <head>
        <meta http-equiv="Content-Type" content="text/html;charset=utf-8">
        <title>PhysBAM Visual Tests</title>
        <link rel="stylesheet" type="text/css" href="style.css">
        </head>
        <body>
        
        """)
        fp.write("<h1>PhysBAM Visual Tests - %s</h1>\n"%self.name)
        fp.write("<h3>%s</h3>\n"%time.ctime())
        fp.write("<table>\n")
        fp.write("<tr><th>Test Name</th><th>Seconds</th><th>Exit Code</th><th>Log</th></tr>\n")
        for test in self.tests:
            test.run()
            test.report(fp)
        fp.write("</table>\n</body>\n</html>\n")
    
