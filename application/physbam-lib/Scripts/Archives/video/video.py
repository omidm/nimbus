#!/usr/bin/python
import Image
import ImageChops
import ImageDraw
import ImageFont
import shutil
import rangeutils
import os
import shutil

def slide(filename,size,lines):
    image_size=640,480
    img=Image.new('RGBA',image_size)
    draw=ImageDraw.Draw(img)

    fontname='/home/aselle/shared/etc/VeraBd.ttf'
    font=ImageFont.truetype(fontname,size)

    # compute line length
    offsets=[]
    widths=[]
    fonts=[]
    colors=[]
    height=0
    default_color='#ffeedd'
    for line in lines:
        color=default_color
        if type(line)==tuple:
            fonts.append(ImageFont.truetype(fontname,line[0]))
            if len(line)>2:
                color=line[2]
            line=line[1]
        else:
            fonts.append(font)
            
        w,h=draw.textsize(line,font=fonts[-1])
        offsets.append(height)
        widths.append(w)
        colors.append(color)
        height+=int(h+.10*h)

    height_start=image_size[1]/2-height/2
    for index in range(len(offsets)):
        if type(lines[index])==tuple:
            line=lines[index][1]
        else:
            line=lines[index]
        draw.text((image_size[0]/2-widths[index]/2,height_start+offsets[index]),line,font=fonts[index],fill=colors[index])

    img.save(filename)

class VIDEO:
    def __init__(self,directory):
        self.directory=directory
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        self.frame=0
        self.fps=30

    def frame_file(self,input_frame):
        return os.path.join(self.directory,"video.%06d.png"%input_frame)

    def add_frame(self,image_filename,count=1):
        #print "adding frame %s with %d copies (index=%d)"%(image_filename,count,self.frame),
        if count>1: print "Adding %-90s %3.1f s"%(image_filename,count*self.fps)
        for i in range(count):
            #print "Copying %s to %s"%(image_filename,self.frame_file(self.frame))
            #os.symlink(image_filename,self.frame_file(self.frame))
            shutil.copy(image_filename,self.frame_file(self.frame))
            self.frame+=1
        #print self.frame

    def add_directory(self,directory,start=-1,end=-1,step=1,duplicate=1):
        ranges=rangeutils.FileRanges(directory)
        files=[]
        for rkey,r in ranges.ranges.items():
            pattern=r.format(False,False)
            if pattern.endswith(".png"):
                print "Adding %s/%s"%(directory,pattern)
                if start==-1 or end==-1:
                    files=map(lambda x:os.path.join(directory,x),r.all_filenames())
                else:
                    files=map(lambda x:os.path.join(directory,r.filename(x)),range(start,end+1))
        print "Adding %-90s %3.1f s"%(directory,len(files)*self.fps)
        count=0
        print step
        for file in files:
            if count%step==0:
                for i in xrange(duplicate):
                    self.add_frame(file)
            count+=1

    def composite_slide_on_current_frame(self,slide_image_name,composite_image_name,count):
        current_frame_image=Image.open(self.frame_file(self.frame-1))
        slide_image=Image.open(slide_image_name)
        mask=ImageChops.invert(slide_image)
        composite_image=Image.composite(current_frame_image,slide_image,mask)
        composite_image.save(composite_image_name)
        for i in range(count):
            shutil.copy(composite_image_name,self.frame_file(self.frame))
            self.frame+=1

    def add_blending(self,source,dest,count=9):
        source_image=Image.open(source)
        dest_image=Image.open(dest)
        for i in range(count):
            blend_image=Image.blend(source_image,dest_image,(i+1.)/(count+1.))
            blend_image.save(self.frame_file(self.frame))
            self.frame+=1
        
    def make_movie(self,filename):
        cmd="ffmpeg -y -r %d -i %s/video.%%06d.png -vcodec mjpeg -b 15000k -r %d %s_mjpeg.avi"%(self.fps,self.directory,self.fps,filename)
        print cmd
        os.system(cmd)
        #cmd="ffmpeg -y -r %d -i %s/video.%%06d.png -vcodec xvid -b 3430k -r %d %s_divx.avi"%(self.fps,self.directory,self.fps,filename)
        #os.system(cmd)
        #cmd="ffmpeg -i %s/video.%%06d.png -vcodec mpeg4 -b 3430k -r %d %s.mov"%(self.directory,self.fps,filename)
        #os.system(cmd)

    
