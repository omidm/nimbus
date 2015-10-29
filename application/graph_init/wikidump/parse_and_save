#!/usr/bin/env python

from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option("-p", "--pagefile", dest="pagefile",
                          help="trimmed pages file (input)")
parser.add_option("-l", "--linkfile", dest="linkfile",
                          help="trimmed links file (input)")
parser.add_option("-j", "--nodenumfile", dest="nodenumfile",
                          help="output file for number of nodes")
parser.add_option("-k", "--edgenumfile", dest="edgenumfile",
                          help="output file for number of edges")
parser.add_option("-n", "--outnodefile", dest="outnodefile",
                          help="output file for nodes")
parser.add_option("-e", "--outedgefile", dest="outedgefile",
                          help="output file for edges (src id, dst id)")

(options, args) = parser.parse_args()


################################################################################
##  Parse pages
################################################################################

pages = dict()
page_ids = dict()
print("Building pages dictionary ...")
pages_file = open(options.pagefile, 'r')
nodes_file = open(options.outnodefile, 'w')
num_nodes = 0
num_all_pages = 0

nodes_done = False
isentry  = False
isstring = False
idparsed = False
nsparsed = False
ttparsed = False
identry = ""
nsentry = ""
ttentry = ""
text    = ""
identry_int = 0
nsentry_int = 0
while not nodes_done:
    c = pages_file.read(1)
    if c == "":
        # done reading file
        nodes_file.write(text)
        nodes_done = True
        text = ""
    elif not isentry:
        if c == '(':
            isentry = True
        else:
            continue
    else:
        if str.isdigit(c) and not idparsed:
            # parse id here
            while str.isdigit(c):
                identry += c
                c = pages_file.read(1)
            idparsed = True
            identry_int = int(identry)
            while c != ',':
                c = pages_file.read(1)
        elif str.isdigit(c) and not nsparsed:
            assert(idparsed)
            # parse namespace here
            while str.isdigit(c):
                nsentry += c
                c = pages_file.read(1)
            nsparsed = True
            nsentry_int = int(nsentry)
            while c != ',':
                c = pages_file.read(1)
            if nsentry_int != 0:
                ttparsed = True
        elif c == "'" and not ttparsed:
            assert(idparsed and nsparsed)
            # parse title here
            c = pages_file.read(1)
            while c != "'":
                ttentry += c
                if c == "\\":
                    # string escapes
                    ttentry += c
                    c = pages_file.read(1) # skip whatever this is
                    ttentry += c
                c = pages_file.read(1)
            ttparsed = True
            while c != ',':
                c = pages_file.read(1)
        # get enclosing ')' for entry -- string entries can have ')'
        # which may be unbalanced?
        elif c == "'" and not isstring:
            assert(idparsed and nsparsed and ttparsed)
            isstring = True
        elif c == "\\" and isstring:
            # string escapes
            c = pages_file.read(1) # skip whatever this is
        elif c == "'" and isstring:
            isstring = False
        elif not isstring and c == ')':
            num_all_pages += 1
            if nsentry_int == 0:
                pages[ttentry] = identry_int
                text = text + ("%d\n" % identry_int)
                page_ids[identry_int] = True
                num_nodes += 1
                if num_nodes % pow(10,6) == 0:
                    nodes_file.write(text)
                    text = ""
                    print("Parsed %d nodes ..." % num_nodes)
            isentry  = False
            idparsed = False
            nsparsed = False
            ttparsed = False
            identry = ""
            nsentry = ""
            ttentry = ""
        else:
            continue
pages_file.close()
nodes_file.close()
print("Built pages dictionary")

nodenum_file = open(options.nodenumfile, 'w')
nodenum_file.write("%s\n" % num_nodes)
nodenum_file.close()

################################################################################
##  Parse links and print edges to file for later partitioning
################################################################################

if not (options.linkfile and options.outedgefile):
    sys.exit(0)

print("Building edges from links...")
links_file = open(options.linkfile, 'r')
edges_file = open(options.outedgefile, 'w')
num_edges = 0
num_all_links = 0

edges_done = False
isentry  = False
isstring = False
idparsed = False
nsparsed = False
ttparsed = False
identry = ""
nsentry = ""
ttentry = ""
text    = ""
identry_int = 0
nsentry_int = 0
while not edges_done:
    c = links_file.read(1)
    if c == "":
        # done reading file
        edges_file.write(text)
        text = ""
        edges_done = True
    elif not isentry:
        if c == '(':
            isentry = True
        else:
            continue
    else:
        if str.isdigit(c) and not idparsed:
            # parse id here
            while str.isdigit(c):
                identry += c
                c = links_file.read(1)
            idparsed = True
            identry_int = int(identry)
            while c != ',':
                c = links_file.read(1)
        elif str.isdigit(c) and not nsparsed:
            assert(idparsed)
            # parse namespace here
            while str.isdigit(c):
                nsentry += c
                c = links_file.read(1)
            nsparsed = True
            nsentry_int = int(nsentry)
            while c != ',':
                c = links_file.read(1)
            if nsentry_int != 0:
                ttparsed = True
        elif c == "'" and not ttparsed:
            assert(idparsed and nsparsed)
            # parse title here
            c = links_file.read(1)
            while c != "'":
                ttentry += c
                if c == "\\":
                    # string escapes
                    ttentry += c
                    c = links_file.read(1) # skip whatever this is
                    ttentry += c
                c = links_file.read(1)
            ttparsed = True
            while c != ',':
                c = links_file.read(1)
        # get enclosing ')' for entry -- string entries can have ')'
        # which may be unbalanced?
        elif c == "'" and not isstring:
            assert(idparsed and nsparsed and ttparsed)
            isstring = True
        elif c == "\\" and isstring:
            # string escapes
            c = links_file.read(1) # skip whatever this is
        elif c == "'" and isstring:
            isstring = False
        elif not isstring and c == ')':
            num_all_links += 1
            if nsentry_int == 0:
                from_id = identry_int
                if from_id in page_ids and ttentry in pages:
                    to_id   = pages[ttentry]
                    text = text +  ("%d %d\n" % (from_id, to_id))
                    num_edges += 1
                    if num_edges % pow(10,6) == 0:
                        edges_file.write(text)
                        text = ""
                        print("Parsed %d edges ..." % num_edges)
            isentry  = False
            idparsed = False
            nsparsed = False
            ttparsed = False
            identry = ""
            nsentry = ""
            ttentry = ""
        else:
            continue
links_file.close()
edges_file.close()
print("Built edges ...")


################################################################################
##  Save number of edges and nodes
################################################################################

edgenum_file = open(options.edgenumfile, 'w')
edgenum_file.write("%s\n" % num_edges)
edgenum_file.close()
