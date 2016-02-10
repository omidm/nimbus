#!/usr/bin/python
######################################################################
# Copyright 2006, Andrew Selle
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# compare.py - parse the log file
######################################################################
import sys
import os
import pickle
import log_parse

data=[]
for file in sys.argv[1:-1]:
    print "Parsing %s..."%file,
    sys.stdout.flush()
    data.append((file,log_parse.Parse(file)))
    print "done."
if os.path.exists(sys.argv[-1]):
    print "not overwriting %s"%sys.argv[-1]
else:
    pickle.dump(data,open(sys.argv[-1],"w"))

