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
                  default="regions_config", type="string",
                  help="input configuration file that lists regions")
parser.add_option("-o", "--output", dest="outfile",
                  default="reg_def", type="string",
                  help="output .h file for generating region defintiions")
parser.add_option("-n", "--namespace", dest="namespace",
                  default="application", type="string",
                  help="namespace to use for generated code")
parser.add_option("-a", "--app_path", dest="app_path",
                  default=curr_path, type="string",
                  help="directory where application and data config file resides")
(options, args)  = parser.parse_args()
reg_config_file  = options.infile           # input region config file
out_str          = options.outfile
out_h_file       = out_str + ".h"           # output .h file for gen code
out_cc_file      = out_str + ".cc"          # output .cc file for gen code
namespace        = options.namespace        # namespace for gen code
app_path         = options.app_path         # app path

# flags
share_boundary_str   = "share_boundary"
generate_inner_str   = "inner"
generate_outer_str   = "outer"
generate_scratch_str = "scratch"

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
    if len(args) != 3:
        print "\nCannot parse line " + str(num) + \
                " because it contains " + str(len(args)-1) + " ':'\n"
        sys.exit(2)
    # Parsing: region name
    names = [x for x in re.split(' |,', args[0]) if x]
    if len(names) != 1:
        print "\nCannot parse line " + str(num) + \
                " because it contains an invalid first argument.\n"
        sys.exit(2)
    reg_name = names[0]
    # Parsing: tuples - domain, number of partitions, ghost width
    sizes_tuple_str = [x for x in re.split('\(|\)', args[1]) if x]
    sizes_str = []
    for tup in sizes_tuple_str:
        temp = [x for x in re.split('\s|,', tup) if x]
        if (len(temp) > 0):
            sizes_str.append(temp)
    params = ValidateSizeTuples(sizes_str, num)
    if len(params) != 3:
        print "\nError parsing line " + str(num) + "\n"
        sys.exit(2)
    flags = re.split('\s|,', args[2].strip())
    share_boundary   = False
    generate_inner   = False
    generate_outer   = False
    generate_scratch = False
    if "none" in flags:
        if (len(flags) > 1):
            print "\nInvalid flags %s at line %i\n" % (str(flags), num)
            sys.exit(2)
    if share_boundary_str in flags:
        share_boundary = True
    if generate_inner_str in flags:
        generate_inner = True
    if generate_outer_str in flags:
        generate_outer = True
    if generate_scratch_str in flags:
        generate_scratch = True
    return reg_name, params, \
            {  share_boundary_str   : share_boundary, \
               generate_inner_str   : generate_inner, \
               generate_outer_str   : generate_outer, \
               generate_scratch_str : generate_scratch }

class Region():
    def __init__(self, x, dx):
        self.x = x
        self.dx = dx
    def __str__(self):
        return "%i, %i, %i, %i, %i, %i" % \
                (self.x[0], self.x[1], self.x[2], \
                self.dx[0], self.dx[1], self.dx[2])

inner   = 'Inner'
outer   = 'Outer'
scratch = 'Scratch'

def GetRegionsBB(params, flags, num, reg_type):
    share_boundary = flags[share_boundary_str]
    domain = params[0]
    rnum   = params[1]
    ghostw = params[2]
    rs     = [0, 0, 0]
    rsl    = [0, 0, 0]
    rd     = [0, 0, 0]
    rdl    = [0, 0, 0]
    for dim in range(0, 3):
        if domain[dim]/rnum[dim] < 2*ghostw[dim]:
            print "\nError : Invalid ghost width and region size"
            print "Region size is smaller than 2 x ghostwidth\n"
            sys.exit(2)
        if reg_type == inner:
            rsl[dim] = domain[dim]/rnum[dim]
            rdl[dim] = rsl[dim]
        elif reg_type == outer:
            rsl[dim] = domain[dim]/rnum[dim] + 2*ghostw[dim]
            rdl[dim] = domain[dim]/rnum[dim]
        if domain[dim]%rnum[dim] == 0:
            rs[dim] = rsl[dim]
            rd[dim] = rdl[dim]
        else:
            print "\nWarning: Regions for dimension " + str(dim) + \
                    " at line " + str(num) + " are not of equal size."
            rs[dim] = rsl[dim]+1
    rsize  = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    rdelta = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    rstart = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    for dim in range(0, 3):
        for i in range(0, rnum[dim], 1):
            rsize[dim][i]  = rs[dim]
            rdelta[dim][i] = rd[dim]
        rsize[dim][rnum[dim]-1]  = rsl[dim]
        rdelta[dim][rnum[dim]-1] = rdl[dim]
        if reg_type == inner:
            rstart[dim][0] = 1
        elif reg_type == outer:
            rstart[dim][0] = 1 - ghostw[dim]
        for i in range(1, rnum[dim], 1):
            rstart[dim][i] = rstart[dim][i-1] + rdelta[dim][i-1]
    regions = []
    if share_boundary:
        for dim in range(0, 3):
            rsize[dim] = map(lambda x : x+1, rsize[dim])
    if reg_type == outer and \
            ghostw[0] == 0 and ghostw[1] == 0 and ghostw[2] == 0:
        return regions
    for i in range(0, rnum[0]):
        for j in range(0, rnum[1]):
            for k in range(0, rnum[2]):
                x  = [rstart[0][i], rstart[1][j], rstart[2][k]]
                dx = [rsize[0][i],  rsize[1][j],  rsize[2][k]]
                if 0 in dx:
                    continue
                r  = Region(x, dx)
                regions.append(r)
    return regions

def GetRegionsScratch(params, flags, num, reg_type):
    share_boundary = flags[share_boundary_str]
    domain = params[0]
    rnum   = params[1]
    ghostw = params[2]
    rs     = [0, 0, 0]
    rsl    = [0, 0, 0]
    for dim in range(0, 3):
        if domain[dim]/rnum[dim] < 2*ghostw[dim]:
            print "\nError : Invalid ghost width and partition size"
            print "Partition size is smaller than 2 x ghostwidth\n"
            sys.exit(2)
        rsl[dim] = domain[dim]/rnum[dim] - 2*ghostw[dim]
        if domain[dim]%rnum[dim] == 0:
            rs[dim] = rsl[dim]
        else:
            print "\nWarning: Partitions for dimension " + str(dim) + \
                    " at line " + str(num) + " are not of equal size."
            rs[dim] = rsl[dim]+1
    rnum = map(lambda x : 3*x+2, rnum)
    rsize  = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    rstart = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    for dim in range(0, 3):
        for i in range(0, rnum[dim], 3) + range(1, rnum[dim], 3):
            rsize[dim][i] = ghostw[dim]
        for i in range(2, rnum[dim], 3):
            rsize[dim][i] = rs[dim]
        rsize[dim][rnum[dim]-3] = rsl[dim]
        rstart[dim][0] = 1-ghostw[dim]
        for i in range(1, rnum[dim], 1):
            rstart[dim][i] = rstart[dim][i-1] + rsize[dim][i-1]
        if share_boundary:  # share_boundary = True
            rsize[dim] = map(lambda x : x+1 if x > 0 else x, rsize[dim])
    ii = [ \
            range(0, pnum[0], 3) + range(1, pnum[0], 3), \
            range(2, pnum[0], 3) ]
    jj = [ \
            range(0, pnum[1], 3) + range(1, pnum[1], 3), \
            range(2, pnum[1], 3) ]
    kk = [ \
            range(0, pnum[2], 3) + range(1, pnum[2], 3), \
            range(2, pnum[2], 3) ]
    regions = []
    # vertex
    for i in ii[0]:
        for j in jj[0]:
            for k in kk[0]:
                x  = [pstart[0][i], pstart[1][j], pstart[2][k]]
                dx = [psize[0][i],  psize[1][j],  psize[2][k]]
                if 0 in dx:
                    continue
                regions.append(Region(x, dx))
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
                    regions.append(Region(x, dx))
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
                    regions.append(Region(x, dx))
    return regions


def GetRegions(params, flags, num, reg_type):
    if reg_type == scratch:
        return GetRegionsScratch(params, flags, num, reg_type)
    else:
        return GetRegionsBB(params, flags, num, reg_type)

## Begin parsing and building information for code generation ##

if not os.path.isfile:
    print "\nCould not find file " + reg_config_file + "\n"
    sys.exit(1)

print "\nReading region configuration file " + reg_config_file + " ..."

reg_config = open(reg_config_file, 'r')
reg_map    = {} # mapping between region name and corresponding set of regions

for num, line in enumerate(reg_config):
    # Parsing: empty lines
    if not line.split():
        continue
    # Parsing: comment
    if line[0] == "#":
        continue
    reg_name, params, flags = ParseLine(line, num)
    # Regions:
    regions = {}
    if flags[generate_inner_str]:
        regions[inner]   = GetRegions(params, flags, num, inner)
    if flags[generate_outer_str]:
        regions[outer] = GetRegions(params, flags, num, outer)
    if flags[generate_scratch_str]:
        regions[scratch] = GetRegions(params, flags, num, scratch)
    if reg_name in reg_map:
        print "\nWarning: Redefinition of " + reg_name + " at line " + str(num)
        print "Ignoring the new definition ..."
    else:
        reg_map[reg_name] = regions

region_num = 0
for rs in reg_map:
    for l in reg_map[rs]:
        if len(reg_map[rs][l]) > 0:
            region_num = region_num + len(reg_map[rs][l])

## Code generation helper functions and variables ##


## Begin code generation ##

print "\nGenerating code ..."

out_h       = open(out_h_file,  'w')
out_cc      = open(out_cc_file, 'w')

# Code generation - .h file

guard_str = "NIMBUS_%s_%s_H_" % \
        (curr_path.upper().replace('/', '_'), out_str.upper())
out_h.write("#ifndef %s\n" % guard_str)
out_h.write("#define %s\n\n" % guard_str)
out_h.write("#include \"shared/geometric_region.h\"\n")
out_h.write("#include \"shared/nimbus.h\"\n")
out_h.write("\nnamespace %s {\n\n" % namespace)
out_h.write("// Bounding boxes (geometric regions) and helper functions for " + \
        "application use\n\n")
for rs in reg_map:
    for l in reg_map[rs]:
        ls = len(reg_map[rs][l])
        if ls > 0:
            out_h.write("extern nimbus::GeometricRegion k%s%s[%i];\n" % (rs, l, ls))
out_h.write("\n")
out_h.write("void InitializeRegions();\n")
out_h.write("\n} // namespace %s\n\n" % namespace) # namespace application
out_h.write("#endif // %s" % guard_str)

# Code generation - .cc file

out_cc.write("#include \"%s/%s\"\n\n" % (app_path, out_h_file))
out_cc.write("#include \"shared/geometric_region.h\"\n")
out_cc.write("#include \"shared/nimbus.h\"\n")
out_cc.write("\nnamespace %s {\n\n" % namespace)
for rs in reg_map:
    for l in reg_map[rs]:
        ls = len(reg_map[rs][l])
        if ls > 0:
            out_cc.write("nimbus::GeometricRegion k%s%s[%i];\n" % (rs, l, ls))
out_cc.write("\n")
out_cc.write("void InitializeRegions() {\n")
for rs in reg_map:
    for l in reg_map[rs]:
        ls = len(reg_map[rs][l])
        if ls > 0:
            out_cc.write("\t// k%s%s\n" % (rs, l))
            for i in range(0, ls, 1):
                out_cc.write("\tk%s%s[%i].Rebuild(%s);\n" % \
                        (rs, l, i, str(reg_map[rs][l][i])))
out_cc.write("}\n") # InitializeRegions
out_cc.write("\n} // namespace %s" % namespace) # namespace application

print ""
