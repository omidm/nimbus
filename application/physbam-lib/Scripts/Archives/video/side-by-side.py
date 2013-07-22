#!/usr/bin/python
import sys, os, re, subprocess
from optparse import OptionParser

usage = "usage: %prog [options] input-dir1 input-dir2 [...] output-dir"
parser = OptionParser(usage="")
parser.add_option("-l", "--left-margin", metavar="LEFT", type="int",
                  help="Crop off LEFT pixels from the left of each "
                  "image and enough from the right to achieve the "
                  "correct size");
parser.add_option("-r", "--right-margin", metavar="RIGHT", type="int",
                  help="Crop off RIGHT pixels from the right of each "
                  "image and enough from the left to achieve the "
                  "correct size");
parser.add_option("-m", "--missing-frame", metavar="FILE",
                  help="Use this image in place of missing input frames");
parser.add_option("-M", "--missing-color", metavar="COLOR",
                  help="Use a blank image of this color in place of missing "
                  "input frames");
parser.add_option("-d", "--divider-width", metavar="PIXELS", type="int",
                  help="Separate images by PIXELS pixels");
parser.add_option("-c", "--divider-color", metavar="COLOR",
                  help="Color to use for separator");
parser.add_option("-p", "--pattern", metavar="PATTERN",
                  default="output.%05d.png",
                  help="Output filename pattern");

(options, args) = parser.parse_args()
if len(args) < 3:
    print "%s: two sources and one destination required"%sys.argv[0]
    print "Try '%s --help' for more information."%sys.argv[0]
    exit(0)

output=args[-1];
args.pop()

r = re.compile("(.*\.)([0-9]*)\.(?:png|jpg|jpeg|bmp|gif)");
files = [{} for x in range(len(args))]
for i in range(len(args)):
    d = args[i]
    for f in os.listdir(d):
        m = r.match(f)
        if m:
            files[i][int(m.groups()[1])] = d+"/"+f

first_image = files[0].itervalues().next()

ident = subprocess.Popen(["convert", "-identify", first_image, "/dev/null"],
                     stdout=subprocess.PIPE).communicate()[0]
r = re.compile("([0-9]+)x([0-9]+)")
m = r.search(ident)
[width, height] = m.groups()
width = int(width)
height = int(height)
try:
    os.makedirs(output)
except:
    pass


eff_width = width
if options.divider_width:
    eff_width -= (len(files) - 1) * options.divider_width

def add_image(s, file, dx, dy, x, y):
    s.append("(")
    s.append(file)
    s.append("-crop")
    s.append("%ix%i+%i+%i"%(dx,dy,x,y))
    s.append(")")

def add_color_image(s, color, dx, dy):
    s.append("(")
    s.append("-size")
    s.append("%ix%i"%(dx,dy))
    s.append("xc:%s"%color)
    s.append(")")

for k in files[0].keys():
    s = ["convert"]
    for a in range(len(files)):
        if a != 0 and options.divider_width:
            add_color_image(s, options.divider_color,
                            options.divider_width, height)

        dx = (a+1)*eff_width/len(files) - a*eff_width/len(files)
        dy = height
        x = (width - dx) / 2
        y = 0
        if options.left_margin:
            x = options.left_margin
        elif options.right_margin:
            x = width - options.right_margin - dx

        try:
            add_image(s, files[a][k], dx, dy, x, y)
        except:
            if options.missing_frame:
                add_image(s, options.missing_frame, dx, dy, x, y)
            elif options.missing_color:
                add_color_image(s, options.missing_color, dx, dy)
            else:
                print "Missing image %i from directory %s"%(k,files[a])
                exit(1)

    s.append("+append")
    s.append(output + "/" + (options.pattern%k))
    p = subprocess.Popen(s);
    p.wait()

