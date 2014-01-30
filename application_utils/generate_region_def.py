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

bool_map = {"true" : True, "false" : False}

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
    share_boundary = args[2].strip()
    if share_boundary not in bool_map:
        print "\nError parsing share_boundary field at line %i\n"  % num
        sys.exit(2)
    return reg_name, params, bool_map[share_boundary]

class Region():
    def __init__(self, x, dx):
        self.x = x
        self.dx = dx
    def __str__(self):
        return "%i, %i, %i, %i, %i, %i" % \
                (self.x[0], self.x[1], self.x[2], \
                self.dx[0], self.dx[1], self.dx[2])

inner   = 'Inner'
outer = 'Outer'

def GetRegions(params, hare_boundary, num, reg_type):
    domain = params[0]
    rnum   = params[1]
    ghostw = params[2]
    rs     = [0, 0, 0]
    rsl    = [0, 0, 0]
    for dim in range(0, 3):
        if domain[dim]/rnum[dim] < 2*ghostw[dim]:
            print "\nError : Invalid ghost width and region size"
            print "Region size is smaller than 2 x ghostwidth\n"
            sys.exit(2)
        if reg_type == inner:
            rsl[dim] = domain[dim]/rnum[dim]
        elif reg_type == outer:
            rsl[dim] = domain[dim]/rnum[dim] + 2*ghostw[dim]
        if domain[dim]%rnum[dim] == 0:
            rs[dim] = rsl[dim]
        else:
            print "\nWarning: Regions for dimension " + str(dim) + \
                    " at line " + str(num) + " are not of equal size."
            rs[dim] = rsl[dim]+1
    rsize  = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    rstart = {0:[0]*rnum[0], 1:[0]*rnum[1], 2:[0]*rnum[2]}
    for dim in range(0, 3):
        for i in range(0, rnum[dim], 1):
            rsize[dim][i] = rs[dim]
        rsize[dim][rnum[dim]-1] = rsl[dim]
        if reg_type == inner:
            rstart[dim][0] = 1
        elif reg_type == outer:
            rstart[dim][0] = 1 - ghostw[dim]
        for i in range(1, rnum[dim], 1):
            rstart[dim][i] = rstart[dim][i-1] + rsize[dim][i-1]
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
                if rsize[0][i] == 0 or rsize[1][j] == 0 or rsize[2][k] == 0:
                    continue
                r  = Region(x, dx)
                regions.append(r)
    return regions

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
    reg_name, params, share_boundary = ParseLine(line, num)
    # Regions:
    regions = {}
    regions[inner]   = GetRegions(params, share_boundary, num, inner)
    regions[outer] = GetRegions(params, share_boundary, num, outer)
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
            out_h.write("nimbus::GeometricRegion k%s%s [%i];\n" % (rs, l, ls))
out_h.write("\n")
out_h.write("void InitializeRegions();\n")
out_h.write("\n} // namespace %s\n\n" % namespace) # namespace application
out_h.write("#endif // %s" % guard_str)

# Code generation - .cc file

out_cc.write("#include \"%s/%s\"\n\n" % (app_path, out_h_file))
out_cc.write("#include \"shared/geometric_region.h\"\n")
out_cc.write("#include \"shared/nimbus.h\"\n")
out_cc.write("\nnamespace %s {\n\n" % namespace)
out_cc.write("void InitializeRegions() {\n")
for rs in reg_map:
    for l in reg_map[rs]:
        ls = len(reg_map[rs][l])
        if ls > 0:
            out_cc.write("\t// k%s_%s\n" % (rs, l))
            for i in range(0, ls, 1):
                out_cc.write("\tk%s%s[%i].Rebuild(%s);\n" % \
                        (rs, l, i, str(reg_map[rs][l][i])))
out_cc.write("}\n") # InitializeRegions
out_cc.write("\n} // namespace %s" % namespace) # namespace application

print ""
