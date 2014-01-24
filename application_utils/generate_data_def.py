#!/usr/bin/env python

import sys
import os.path
import re
from optparse import OptionParser


## Parse the command line arguments ##

parser = OptionParser()
parser.add_option("-i", "--input", dest="infile",
                  default="data_config", type="string",
                  help="input configuration file listing data configuration")
parser.add_option("-o", "--output", dest="outfile",
                  default="data_def", type="string",
                  help="output .h|cc file for generating data defintiions")
(options, args) = parser.parse_args()
data_config_file = options.infile           # input data config file
out_h_file       = options.outfile + ".h"   # output .h file for gen code
out_cc_file      = options.outfile + ".cc"  # output .cc file for gen code


## Parsing helper functions ##

def ValidateSizeTuples(sizes, num):
    if len(sizes) != 3:
        print "Expected 3 tuples - domain, number of partitions, ghost width"
        print "Got " + str(len(sizes)) + " at line " + str(num)
        return []
    params = []
    for tup in sizes:
        if len(tup) != 3:
            print "Expect tuple elements of type (int, int, int)"
            print "Got " + str(tup) + " at line " + str(num)
            return []
        param_tup = []
        for t in tup:
            try:
                param_tup.append(int(t))
            except:
                print "Expect tuple elements of type (int, int, int)"
                print "Got " + str(tup) + " at line " + str(num)
                return []
        params.append(param_tup)
    return params

def ParseLine(line, num):
    # Parsing: begin parsing
    args = re.split(':', line)
    if len(args) != 3:
        print "\nCannot parse line " + str(num) + \
                " because it contains only " + str(len(args)-1) + " ':'\n"
        sys.exit(2)
    # Parsing: C++ class
    cpp_class = args[0]
    # Parsing: nimbus types
    nimbus_types = [x for x in re.split(' |,', args[1]) if x]
    # Parsing: tuples - domain, number of partitions, ghost width
    sizes_tuple_str = [x for x in re.split('\(|\)', args[2]) if x]
    sizes_str = []
    for tup in sizes_tuple_str:
        temp = [x for x in re.split('\s|,', tup) if x]
        if (len(temp) > 0):
            sizes_str.append(temp)
    params = ValidateSizeTuples(sizes_str, num)
    if len(params) == 0:
        print "\nError parsing line " + str(num) + "\n"
        sys.exit(2)
    return cpp_class, nimbus_types, params

class Partition():
    def __init__(self, x, dx):
        self.x = x
        self.dx = dx
    def toString(self):
        return "%i, %i, %i, %i, %i, %i" % \
                (self.x[0], self.x[1], self.x[2], \
                self.dx[0], self.dx[1], self.dx[2])

def GetPartitionIds(partn_id, params, num):
    domain = params[0]
    pnum   = params[1]
    ghostw = params[2]
    ps     = [0, 0, 0]
    psl    = [0, 0, 0]
    psize  = {0:[], 1:[], 2:[]}
    pstart = {0:[], 1:[], 2:[]}
    for dim in range(0, 3):
        psl[dim] = domain[dim]/pnum[dim]
        if domain[dim]%pnum[dim] == 0:
            ps[dim] = psl[dim]
        else:
            print "Warning: Partitions for dimension " + str(dim) + \
                    " at line " + str(num) + " are not of equal size."
            ps[dim] = psl[dim]+1
    pnum   = map(lambda x : x+2, pnum)
    for dim in range(0, 3):
        pstart[dim].append(1-ghostw[dim])
        psize[dim].append(ghostw[dim])
        for i in range(1, pnum[dim]-2):
            pstart[dim].append(pstart[dim][i-1] + psize[dim][i-1])
            psize[dim].append(ps[dim])
        l = pnum[dim]-2
        pstart[dim].append(pstart[dim][l-1] + psize[dim][l-1])
        psize[dim].append(psl[dim])
        l = l + 1
        pstart[dim].append(pstart[dim][l-1] + psize[dim][l-1])
        psize[dim].append(ghostw[dim])
    pidset = set()
    for i in range(0, pnum[0]):
        for j in range(0, pnum[1]):
            for k in range(0, pnum[2]):
                x    = [pstart[0][i], pstart[1][j], pstart[2][k]]
                dx   = [psize[0][i],  psize[1][j],  psize[2][k]]
                p    = Partition(x, dx)
                pstr = p.toString()
                if pstr in partn_id:
                    pidset.add(partn_id[pstr])
                else:
                    l = len(partn_id)
                    partn_pid[pstr] = l
                    pidset.add(l)
    return pidset

## Begin parsing and building information for code generation ##

if not os.path.isfile:
    print "Could not find file " + data_config_file + "\n"
    sys.exit(1)

print "\nReading data configuration file " + data_config_file + " ...\n"
data_config = open(data_config_file, 'r')
ntypes_pid  = {} # mapping between nimbus type and partition id set
partn_pid   = {} # mapping between partition and partition id

for num, line in enumerate(data_config):
    # Parsing: empty lines
    if not line.split():
        continue
    # Parsing: comment
    if line[0] == "#":
        continue
    cpp_class, nimbus_types, params = ParseLine(line, num)
    print cpp_class + " -- " + str(nimbus_types) + " -- " + str(params)
    # Data for code generation:
    partns = GetPartitionIds(partn_pid, params, num)
    for nt in nimbus_types:
        if nt in ntypes_pid:
            print "Redefinition of " + nt + " at line " + str(num)
            print "Ignoring the new definition ..."
        else:
            ntypes_pid[nt] = partns

data_num = 0
for d in ntypes_pid:
    data_num = data_num + len(ntypes_pid[d])
part_num = len(partn_pid)

# TODO: complete beyond this
print partn_pid
print ntypes_pid


## Code generation helper functions ##

logical_id_vector_str = "std::vector<nimbus::logical_data_id_t>"
partition_id_str      = "nimbus::ID<partition_id_t>"


## Begin code generation ##

out_h       = open(out_h_file, 'w')
out_cc      = open(out_cc_file, 'w')

# Code generation - .h file

out_h.write("#include \"shared/nimbus.h\"\n")
out_h.write("#include \"shared/geometric_region.h\"\n\n")
out_h.write("namespace application {\n\n")
out_h.write("// Helper functions for defining data objects\n")
out_h.write(logical_id_vector_str + " DefineNimbusData(nimbus::Job *jb);\n\n")
out_h.write("// Constants - partition names\n")
for p in sorted(partn_pid.items(), key=lambda x: x[1]):
    decl = "const nimbus::GeometricRegion %s(%s);\n" % \
            ("kRegion%i" % (p[1]), p[0])
    out_h.write(decl)
out_h.write("\n} // namespace application") # namespace application

# Code generation - .cc file

out_cc.write("#include \"shared/nimbus.h\"\n")
out_cc.write("#include \"shared/geometric_region.h\"\n\n")
out_cc.write("namespace application {\n\n")

out_cc.write(logical_id_vector_str + " DefineNimbusData(nimbus::Job *jb) {\n")

# partitions
out_cc.write("\n")
out_cc.write("\t// Define partitions\n")
out_cc.write("\tnimbus::Parameter part_params;\n")
for p in sorted(partn_pid.items(), key=lambda x: x[1]):
    out_cc.write("\t%s pid%i(%i);\n" % (partition_id_str, p[1], p[1]))
    out_cc.write("\tjb->DefinePartition(%s, %s, %s);\n" % \
            ("pid" + str(p[1]), "kRegion" + str(p[1]), "part_params"))
# data
out_cc.write("\n")
out_cc.write("\t//Data setup\n")
out_cc.write("\tint data_num = %i;\n" % data_num)
out_cc.write("\t%s data_ids;\n" % logical_id_vector_str)

out_cc.write("}\n") # DefineNimbusData

out_cc.write("\n} // namespace application") # namespace application
