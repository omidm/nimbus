#!/usr/bin/python
######################################################################
# Copyright 2006, Andrew Selle
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################
# display_cached.py
######################################################################
import pickle
import sys
import log_parse

name,data=pickle.load(open(sys.argv[1]))[int(sys.argv[2])]
print "-------------------------------------------------------------------------------------------------------"
print "DATA SET %s"%name
print "-------------------------------------------------------------------------------------------------------"
data.show()

