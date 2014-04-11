#!/usr/bin/env python

import sys
import re
from decimal import *



def find_hash_table(file_name, tag):
  f = open(file_name, 'r')
  content = f.readlines();
 
  count = 0;
  table = {}
  regexp =  '.*' + tag + '.*Compute:(\w+)\s*id:\s*(\d+)\s*aggregate_hash:\s*(\d+).*'
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


def compare_hash_tables(table_ref, table):
  match = True
  for key in table.keys():
    if not table_ref.has_key(key):
      print "ERROR: could not find key in refernce table: " + key
      match = False
    elif table_ref[key] != table[key]:
      print "ERROR: the value does not match for the key: " + key
      match = False
  return match


table_ref = find_hash_table(sys.argv[1], "hash_in")
table = find_hash_table(sys.argv[2], "hash_in")

if compare_hash_tables(table_ref, table):
  print "PASS: data hash tables match."
else:
  print "FAIL: data hash tables do not match."




