#!/usr/bin/env python

import sys
import re
import numpy
from decimal import *






regexp = '.*Advect V.*time value=\"(\d+\.\d+).*'
f = open(sys.argv[1], 'r')

# f = open('log.txt', 'r')
content = f.readlines()

force = []
advect_phi = []
advect_removed_particles = []

collect = False
for line in content:
  if not collect:
    result = re.findall('.*Calculate Dt.*', line)
    if len(result) > 0:
      force_time = 0
      advect_phi_time = 0
      advect_removed_particles_time = 0
      collect = True
      continue
  else:
    result = re.findall('.*Project.*', line)
    if len(result) > 0:
      force.append(force_time)
      advect_phi.append(advect_phi_time)
      advect_removed_particles.append(advect_removed_particles_time)
      collect = False
      continue

    result = re.findall('.*Forces\".*(\d+\.\d+)', line)
    if len(result) > 0:
      force_time += float(result[0])
      continue

    result = re.findall('.*Advect Phi\".*(\d+\.\d+)', line)
    if len(result) > 0:
      advect_phi_time += float(result[0])
      continue

    result = re.findall('.*Advect Removed Particles\".*(\d+\.\d+)', line)
    if len(result) > 0:
      advect_removed_particles_time += float(result[0])
      continue


print 'Force:'
print numpy.mean(force)
print numpy.std(force)

print 'Advect Phi:'
print numpy.mean(advect_phi)
print numpy.std(advect_phi)

print 'Advect Removed Particles:'
print numpy.mean(advect_removed_particles)
print numpy.std(advect_removed_particles)




















def find_version_table(file_name, tag):
  f = open(file_name, 'r')
  content = f.readlines();
 
  count = 0;
  table = {}
  regexp =  '.*' + tag + '.*Compute:(\w+)\s*id:\s*(\d+)\s*version_hash:\s*(\d+).*'
  for line in content:
    result = re.findall(regexp, line)
    if (len(result) >= 1):
      table[(result[0][0] + "-" + result[0][1])] = result[0][2]
      count = count + 1
    # else:
    #   print "ERROR: corrupted line in " + file_name
    #   print line
  print "Read " + str(count) + " valid lines in " + file_name
  return table


def compare_version_tables(table_ref, table):
  match = True
  for key in table.keys():
    if not table_ref.has_key(key):
      print "ERROR: could not find key in refernce table: " + key
      match = False
    elif table_ref[key] != table[key]:
      print "ERROR: the value does not match for the key: " + key
      match = False
  return match




