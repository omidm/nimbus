#!/usr/bin/env python

import sys
import re
from decimal import *



def find_version_table(file_name, tag):
  f = open(file_name, 'r')
  content = f.readlines();
 
  count = 0;
  table = {}
  regexp =  '.*' + tag + '.*Compute.*id:\s*(\d+)\s*version_hash:\s*(\d+).*'
  for line in content:
    result = re.findall(regexp, line)
    if (len(result) >= 1):
      table[result[0][0]] = result[0][1]
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


table_ref = find_version_table(sys.argv[1], "version_in")
table = find_version_table(sys.argv[2], "version_in")

if compare_version_tables(table_ref, table):
  print "PASS: version tables match."
else:
  print "FAIL: version tables do not match."




