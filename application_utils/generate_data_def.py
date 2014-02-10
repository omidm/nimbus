#!/usr/bin/env python

import sys
import os
import os.path
import re
from optparse import OptionParser


## Parse the command line arguments ##

cwd = os.getcwd()
temp = cwd[cwd.find("nimbus/application/"):]
curr_path = temp[temp.find("application/"):]

parser = OptionParser()
parser.add_option("-i", "--input", dest="infile",
                  default="data_config", type="string",
                  help="input configuration file listing data configuration")
parser.add_option("-o", "--output", dest="outfile",
                  default="data_def", type="string",
                  help="output .h|cc file for generating data defintiions")
parser.add_option("-n", "--namespace", dest="namespace",
                  default="application", type="string",
                  help="namespace to use for generated code")
parser.add_option("-a", "--app_path", dest="app_path",
                  default=curr_path, type="string",
                  help="directory where application and data config file resides")
(options, args) = parser.parse_args()
data_config_file = options.infile           # input data config file
out_str          = options.outfile
out_h_file       = out_str + ".h"           # output .h file for gen code
out_cc_file      = out_str + ".cc"          # output .cc file for gen code
namespace        = options.namespace        # namespace for gen code
app_path         = options.app_path         # app path

# flags
include_boundary_str = "include_boundary"
share_boundary_str   = "share_boundary"
use_scratch_str      = "use_scratch"

# scratch types
vertex = "vertex"
edge   = "edge"
face   = "face"

## Parsing helper functions ##

def ValidateSizeTuples(sizes, num):
    if len(sizes) != 3:
        print "\nExpected 3 tuples - domain, number of partitions, ghost width"
        print "Got " + str(len(sizes)) + " at line " + str(num)
        return []
    params = []
    for tup in sizes:
        if len(tup) != 3:
            print "\nExpected tuple elements of type (int, int, int)"
            print "Got " + str(tup) + " at line " + str(num)
            return []
        param_tup = []
        for t in tup:
            try:
                param_tup.append(int(t))
            except:
                print "\nExpect tuple elements of type (int, int, int)"
                print "Got " + str(tup) + " at line " + str(num)
                return []
        params.append(param_tup)
    return params

def ParseLine(line, num):
    # Parsing: begin parsing
    args = re.split(':', line)
    if len(args) != 4:
        print "\nCannot parse line " + str(num) + \
                " because it contains " + str(len(args)-1) + " ':'\n"
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
    if len(params) != 3:
        print "\nError parsing line " + str(num) + "\n"
        sys.exit(2)
    flags = re.split('\s|,', args[3].strip())
    include_boundary = False
    share_boundary   = False
    use_scratch      = False
    if "none" in flags:
        if (len(flags) > 1):
            print "\nInvalid flags %s at line %i\n" % (str(flags), num)
            sys.exit(2)
    if include_boundary_str in flags:
        include_boundary = True
    if share_boundary_str in flags:
        share_boundary = True
    if use_scratch_str in flags:
        use_scratch = True
    return cpp_class, nimbus_types, params, \
            {  include_boundary_str : include_boundary, \
               share_boundary_str   : share_boundary, \
               use_scratch_str      : use_scratch }
            
class Partition():
    def __init__(self, x, dx):
        self.x = x
        self.dx = dx
    def toString(self):
        return "%i, %i, %i, %i, %i, %i" % \
                (self.x[0], self.x[1], self.x[2], \
                self.dx[0], self.dx[1], self.dx[2])

def GetPartitionID(x, dx, partn_id):
    p = Partition(x, dx)
    pstr = p.toString()
    if pstr in partn_id:
        return partn_id[pstr]
    else:
        l = len(partn_id)
        partn_id[pstr] = l
        return l

def GetNimbusTypePIDMapHelper(nimbus_types, \
                              pnum, pstart, psize, use_scratch, \
                              partn_id):
    nt_pid = {}
    pidset = set()
    for i in range(0, pnum[0]):
        for j in range(0, pnum[1]):
            for k in range(0, pnum[2]):
                x  = [pstart[0][i], pstart[1][j], pstart[2][k]]
                dx = [psize[0][i],  psize[1][j],  psize[2][k]]
                if 0 in dx:
                    continue
                pidset.add(GetPartitionID(x, dx, partn_id))
    for nt in nimbus_types:
        nt_pid[nt] = pidset
    if use_scratch:
        ii = [ \
                range(0, pnum[0], 3) + range(1, pnum[0], 3), \
                range(2, pnum[0], 3) ]
        jj = [ \
                range(0, pnum[1], 3) + range(1, pnum[1], 3), \
                range(2, pnum[1], 3) ]
        kk = [ \
                range(0, pnum[2], 3) + range(1, pnum[2], 3), \
                range(2, pnum[2], 3) ]
        for nt in nimbus_types:
            # vertex
            vpidset = set()
            for i in ii[0]:
                for j in jj[0]:
                    for k in kk[0]:
                        x  = [pstart[0][i], pstart[1][j], pstart[2][k]]
                        dx = [psize[0][i],  psize[1][j],  psize[2][k]]
                        if 0 in dx:
                            continue
                        vpidset.add(GetPartitionID(x, dx, partn_id))
            for t in range(1, 9):
                nt_pid["%s_%s_%i" % (nt, vertex, t)] = vpidset
            # edge
            epidset = set()
            for dim in range(0, 3):
                for i in ii[1 if dim == 0 else 0]:
                    for j in ii[1 if dim == 1 else 0]:
                        for k in kk[1 if dim == 2 else 0]:
                            x  = [pstart[0][i], pstart[1][j], pstart[2][k]]
                            dx = [psize[0][i],  psize[1][j],  psize[2][k]]
                            if 0 in dx:
                                continue
                            epidset.add(GetPartitionID(x, dx, partn_id))
            for t in range(1, 5):
                nt_pid["%s_%s_%i" % (nt, edge, t)] = epidset
            # face
            fpidset = set()
            for dim in range(0, 3):
                for i in ii[0 if dim == 0 else 1]:
                    for j in ii[0 if dim == 1 else 1]:
                        for k in kk[0 if dim == 2 else 1]:
                            x  = [pstart[0][i], pstart[1][j], pstart[2][k]]
                            dx = [psize[0][i],  psize[1][j],  psize[2][k]]
                            if 0 in dx:
                                continue
                            fpidset.add(GetPartitionID(x, dx, partn_id))
            for t in range(1, 2):
                nt_pid["%s_%s_%i" % (nt, face, t)] = fpidset
    return nt_pid

def GetNimbusTypePIDMap(nimbus_types, params, flags, num, partn_id):
    include_boundary = flags[include_boundary_str]
    share_boundary   = flags[share_boundary_str]
    use_scratch      = flags[use_scratch_str]
    domain = params[0]
    pnum   = params[1]
    ghostw = params[2]
    ps     = [0, 0, 0]
    psl    = [0, 0, 0]
    for dim in range(0, 3):
        if domain[dim]/pnum[dim] < 2*ghostw[dim]:
            print "\nError : Invalid ghost width and partition size"
            print "Partition size is smaller than 2 x ghostwidth\n"
            sys.exit(2)
        psl[dim] = domain[dim]/pnum[dim] - 2*ghostw[dim]
        if domain[dim]%pnum[dim] == 0:
            ps[dim] = psl[dim]
        else:
            print "\nWarning: Partitions for dimension " + str(dim) + \
                    " at line " + str(num) + " are not of equal size."
            ps[dim] = psl[dim]+1
    pnum = map(lambda x : 3*x+2, pnum)
    # fill in pstart and psize that define geometric regions corresponding to
    # partitions
    psize  = {0:[0]*pnum[0], 1:[0]*pnum[1], 2:[0]*pnum[2]}
    pstart = {0:[0]*pnum[0], 1:[0]*pnum[1], 2:[0]*pnum[2]}
    for dim in range(0, 3):
        for i in range(0, pnum[dim], 3) + range(1, pnum[dim], 3):
            psize[dim][i] = ghostw[dim]
        for i in range(2, pnum[dim], 3):
            psize[dim][i] = ps[dim]
        psize[dim][pnum[dim]-3] = psl[dim]
        pstart[dim][0] = 1-ghostw[dim]
        if not include_boundary:  # include_boundary = False
            pstart[dim][0] = 1
            psize[dim][0]  = 0
            psize[dim][pnum[dim]-1] = 0
        for i in range(1, pnum[dim], 1):
            pstart[dim][i] = pstart[dim][i-1] + psize[dim][i-1]
        if share_boundary:  # share_boundary = True
            psize[dim] = map(lambda x : x+1 if x > 0 else x, psize[dim])
    # obtain nimbus type - partition id mapping
    nt_pid = GetNimbusTypePIDMapHelper(nimbus_types, \
                                       pnum, pstart, psize, use_scratch, \
                                       partn_id)
    return nt_pid


## Begin parsing and building information for code generation ##

if not os.path.isfile:
    print "\nCould not find file " + data_config_file + "\n"
    sys.exit(1)

print "\nReading data configuration file " + data_config_file + " ..."

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
    cpp_class, nimbus_types, params, flags = ParseLine(line, num)
    # Data for code generation:
    nt_pid = GetNimbusTypePIDMap(nimbus_types, params, flags,  num, partn_pid)
    for nt in nt_pid:
        if nt in ntypes_pid:
            print "\nWarning: Redefinition of " + nt + " at line " + str(num)
            print "Ignoring the new definition ..."
        else:
            ntypes_pid[nt] = nt_pid[nt]

data_num = 0
for d in ntypes_pid:
    data_num = data_num + len(ntypes_pid[d])
part_num = len(partn_pid)


## Code generation helper functions and variables ##

logical_id_vector_str = "std::vector<nimbus::logical_data_id_t>"
partition_id_set_str      = "nimbus::ID<partition_id_t>"


## Begin code generation ##

print "\nGenerating code ..."

out_h       = open(out_h_file, 'w')
out_cc      = open(out_cc_file, 'w')

# Code generation - .h file

guard_str = "NIMBUS_%s_%s_H_" % \
        (curr_path.upper().replace('/', '_'), out_str.upper())
out_h.write("#ifndef %s\n" % guard_str)
out_h.write("#define %s\n\n" % guard_str)
out_h.write("#include \"shared/nimbus.h\"\n")
out_h.write("\nnamespace %s {\n\n" % namespace)
out_h.write("// Helper functions for defining data objects\n")
out_h.write(logical_id_vector_str + " DefineNimbusData(nimbus::Job *jb);\n")
out_h.write("\n} // namespace %s\n\n" % namespace) # namespace application
out_h.write("#endif // %s" % guard_str)

# Code generation - .cc file

out_cc.write("#include \"%s/%s\"\n\n" % (app_path, out_h_file))
out_cc.write("#include \"shared/geometric_region.h\"\n")
out_cc.write("#include \"shared/nimbus.h\"\n")
out_cc.write("\nnamespace %s {\n\n" % namespace)

out_cc.write(logical_id_vector_str + " DefineNimbusData(nimbus::Job *jb) {\n\n")

# geometric regions
pt_num_str = "pnum"
pt_num     = len(partn_pid)
out_cc.write("\t// Constants - geometric regions\n")
out_cc.write("\tint %s = %i;\n" % (pt_num_str, pt_num))
out_cc.write("\tnimbus::GeometricRegion kRegions[%s];\n" % pt_num_str)
for p in sorted(partn_pid.items(), key=lambda x: x[1]):
    decl = "\tkRegions[%i] = nimbus::GeometricRegion(%s);\n" % (p[1], p[0])
    out_cc.write(decl)

# partitions
partitions_str = "partitions"
out_cc.write("\n")
out_cc.write("\t// Define partitions\n")
out_cc.write("\t%s %s[%s];\n" % (partition_id_set_str, partitions_str, pt_num_str))
out_cc.write("\tnimbus::Parameter part_params;\n")
out_cc.write("\tfor (int i = 0; i < %s; i++) {\n" % pt_num_str)
out_cc.write("\t\t%s[i] = %s(i);\n" % (partitions_str, partition_id_set_str))
decl = "\t\tjb->DefinePartition(%s, %s, %s);\n" %\
        (("%s[i]" % partitions_str), "kRegions[i]", "part_params")
out_cc.write(decl)
out_cc.write("\t}\n")

# data setup
data_ids_str = "data_ids"
data_num_str = "data_num"
out_cc.write("\n")
out_cc.write("\t// Data setup\n")
out_cc.write("\tint %s = %i;\n" % (data_num_str, data_num))
out_cc.write("\t%s %s;\n" %\
        (logical_id_vector_str, data_ids_str))
out_cc.write("\tjb->GetNewLogicalDataID(&%s, %s);\n" %\
        (data_ids_str, data_num_str))
ids_used_str = "data_ids_used"
ids_used     = 0
out_cc.write("\tint %s = %i;\n" % (ids_used_str, ids_used))
ids_to_use_str = "ids_to_use"
ids_to_use     = 0
out_cc.write("\tint %s = %i;\n" % (ids_to_use_str, ids_to_use))
np_str = "neighbor_partitions"
out_cc.write("\tnimbus::IDSet<nimbus::partition_id_t> %s;\n" % np_str)
dp_str = "data_params"
out_cc.write("\tnimbus::Parameter %s;\n" % dp_str)

# define data
for d in ntypes_pid:
    out_cc.write("\n")
    out_cc.write("\t// Define data %s\n" % d)
    ids_to_use     = len(ntypes_pid[d])
    out_cc.write("\t%s = %i;\n" % (ids_to_use_str, ids_to_use))
    pset_str = "pset_" + d
    temp = [x for x in re.split('\[|\]', str(list(ntypes_pid[d]))) if x]
    decl = "size_t " + pset_str + "[%i] = {%s}" %\
            (ids_to_use, temp[0])
    out_cc.write("\t%s;\n" % decl)
    out_cc.write("\tfor (int i = 0; i < %s; i++) {\n" % ids_to_use_str)
    out_cc.write("\t\tjb->DefineData(\"%s\", %s[%s + i], %s.elem(), %s, %s);\n" % \
            (d, data_ids_str, ids_used_str, \
            ("%s[%s[i]]" % (partitions_str, pset_str)), np_str, dp_str))
    out_cc.write("\t}\n")
    out_cc.write("\t%s += %s;\n" % (ids_used_str, ids_to_use_str))

out_cc.write("\n\treturn(%s);\n" % data_ids_str)
out_cc.write("}\n") # DefineNimbusData

out_cc.write("\n} // namespace %s" % namespace) # namespace application

print ""
