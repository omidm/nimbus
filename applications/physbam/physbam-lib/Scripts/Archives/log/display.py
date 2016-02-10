#!/usr/bin/python
######################################################################
# Copyright 2006, Andrew Selle
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# display.py - Prints summed final timing stats of simulation
######################################################################
from optparse import OptionParser
import sys
import log_parse

usage="usage: %prog [options] <log fileanme>"
parser=OptionParser(usage)
parser.add_option("-x","--summary-only",action="store_true",dest="summary_only",default=False)
parser.add_option("-1","--start_frame",action="store",type="int",dest="start_frame",default=None)
parser.add_option("-f","--stop_frame",action="store",type="int",dest="stop_frame",default=None)
parser.add_option("-p","--noprints",action="store_false",dest="prints",default=True)
parser.add_option("-s","--nostats",action="store_false",dest="stats",default=True)
parser.add_option("-t","--noscopes",action="store_false",dest="scopes",default=True)
(options,args)=parser.parse_args()
if len(args)!=1:
    parser.error("incorrect number of arguments")
    sys.exit(1)

log_parse.print_prints=options.prints and not options.summary_only
log_parse.print_stats=options.stats and not options.summary_only
log_parse.print_scopes=options.scopes and not options.summary_only
root=log_parse.Parse(args[0],options.start_frame,options.stop_frame)
print "--- SUMMARY -----------------------------------------------------------------------------------------------------------------"
root.show()
