#!/usr/bin/env python

from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option("-n", "--nodefile", dest="nodefile",
                          help="input node file")
parser.add_option("-e", "--edgefile", dest="edgefile",
                          help="input edges file (src id, dst id)")
parser.add_option("-o", "--outfile", dest="outfile",
                          help="output file (input to parmets)")

(options, args) = parser.parse_args()


################################################################################
##  Parse pages
################################################################################

node_map = dict()
edges = dict()

LOG_STEP = 1000000

print("Building node map ...")

nodefile = open(options.nodefile, 'r')
node_num = 0
for line in nodefile:
    node_num += 1
    node_map[line.strip()] = node_num
    edges[node_num] = []
    if (node_num % LOG_STEP == 0):
        print "Processed " + str(node_num) + " nodes..."
nodefile.close()

print("Building graph for metis ...")

edgefile = open(options.edgefile, 'r')
edge_num = 0
for line in edgefile:
    edge_num += 1
    words = str.split(line)
    src = node_map[words[0]]
    dst = node_map[words[1]]
    edges[src].append(dst)
    edges[dst].append(src)
    if (edge_num % LOG_STEP == 0):
        print "Processed " + str(edge_num) + " edges..."
edgefile.close()

print("Saving output for metis ...")

counter = 0;
edge_num = 0
for e in edges:
    unique = set(edges[e])
    counter += len(unique)
    edge_num += len(unique)
    if (counter > LOG_STEP):
        counter -= LOG_STEP
        print "Conted " + str(edge_num) + " lines..."

outfile = open(options.outfile, 'w')
outfile.write("%d %d 000\n" % (node_num, edge_num/2))

counter = 0;
edge_num = 0
for e in edges:
    unique = set(edges[e])
    counter += len(unique)
    edge_num += len(unique)
    text = " ".join(map(str, unique))
    text += "\n"
    outfile.write(text)
    if (counter > LOG_STEP):
        counter -= LOG_STEP
        print "Wrote " + str(edge_num) + " lines..."
outfile.close()

