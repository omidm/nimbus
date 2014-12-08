#!/usr/bin/env python 
import sys

sd_check_lc = []
sd_total_lc = []
sd_check_rc = []
sd_total_rc = []
pull_data_lc = []
pull_data_rc = []
im_check_lc = []
im_total_lc = []
im_check_rc = []
im_total_rc = []
gav_wfc = []
gav_rfc = []
gav_check = []
gav_total = []
gas_wfc = []
gas_rfc = []
gas_check = []
gas_total = []

for i in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    fname = str(i) + "_parse.txt"
    f = open(fname)
    pd_lc = 0
    pd_rc = 0
    wfc = 0
    for line in f:
        if "IM LC" in line:
            im_total_lc.append(float(line.split()[-1].strip(',')))
        elif "IM RCR" in line:
            im_total_rc.append(float(line.split()[-1].strip(',')))
        elif "IM check LC" in line:
            im_check_lc.append(float(line.split()[-1].strip(',')))
        elif "IM check RCR" in line:
            im_check_rc.append(float(line.split()[-1].strip(',')))
        elif "SD LC" in line:
            sd_total_lc.append(float(line.split()[-1].strip(',')))
        elif "SD RCS" in line:
            sd_total_rc.append(float(line.split()[-1].strip(',')))
        elif "SD check LC" in line:
            sd_check_lc.append(float(line.split()[-1].strip(',')))
        elif "SD check RCS" in line:
            sd_check_rc.append(float(line.split()[-1].strip(',')))
        elif "pdata" in line and "LC" in line:
            pd_lc += float(line.split()[-1].strip(','))
        elif "pdata" in line and "RCS" in line:
            pd_rc += float(line.split()[-1].strip(',').strip('}'))
        elif "'GAV'" in line:
            gav_total.append(float(line.split()[-1].strip(',')))
        elif "GAV check" in line:
            gav_check.append(float(line.split()[-1].strip(',')))
        elif "GAV rfc" in line:
            gav_rfc.append(float(line.split()[-1].strip(',')))
        elif "GAV wfc" in line:
            wfc += float(line.split()[-1].strip(','))
        elif "'GAS'" in line:
            gas_total.append(float(line.split()[-1].strip(',')))
        elif "GAS check" in line:
            gas_check.append(float(line.split()[-1].strip(',')))
        elif "GAS rfc" in line:
            gas_rfc.append(float(line.split()[-1].strip(',')))
        elif "GAS wfc" in line:
            gas_wfc.append(float(line.split()[-1].strip(',')))
    pull_data_lc.append(pd_lc)
    pull_data_rc.append(pd_rc)
    gav_wfc.append(wfc)

print "SD CHECK LC"
print sd_check_lc
print
print "SD_TOTAL_LC"
print sd_total_lc
print
print "SD_CHECK_RC"
print sd_check_rc
print
print "SD_TOTAL_RC"
print sd_total_rc
print
print "PULL_DATA_LC"
print pull_data_lc
print 
print "PULL_DATA_RC"
print pull_data_rc
print
print "IM_CHECK_LC"
print im_check_lc
print
print "IM_TOTAL_LC"
print im_total_lc
print
print "IM_CHECK_RC"
print im_check_rc
print
print "IM_TOTAL_RC"
print im_total_rc
print
print "GAV_WFC"
print gav_wfc
print
print "GAV_RFC"
print gav_rfc
print
print "GAV_CHECK"
print gav_check
print
print "GAV_TOTAL"
print gav_total
print
print "GAS_WFC"
print gas_wfc
print
print "GAS_RFC"
print gas_rfc
print
print "GAS_CHECK"
print gas_check
print
print "GAS_TOTAL"
print gas_total
