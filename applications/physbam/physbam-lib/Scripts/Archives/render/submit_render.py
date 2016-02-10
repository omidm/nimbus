import sys
import os
import string

if len(sys.argv) < 6:
    print "<scene file> <input_data_folder> <output_folder> <start frame> <end frame> <ray-tracing binary>"
    sys.exit()

file = open(sys.argv[1])
total_lines = ""

#read in skeleton scene file for rendering
while 1:
    line = file.readline()
    total_lines = total_lines + line
    if not line:
        break
    pass

tmp_string = total_lines.replace('<input_data>', sys.argv[2])
final_string = tmp_string.replace('<output_folder>', sys.argv[3])

file.close()

#write scene file to render
os.system("mkdir tmp_files")
tmp_directory = os.getcwd() + "/tmp_files/"
tmp_filename = "tmp_output_render1234.scene"
file = open(tmp_directory + tmp_filename, "w")
file.write(final_string)
file.close()

#now generate the qsub scripts
for i in range(string.atoi(sys.argv[4]), string.atoi(sys.argv[5])):
#    os.system("rm -f tmp_submit_script" + str(i - 1) + ".sh")    
    submit_file = open(tmp_directory + "tmp_submit_script" + str(i) + ".sh", "w")
    submit_file.write("#!/bin/sh\n")
    submit_file.write("#PBS -o pbsout" + str(i) + "\n")    
    submit_file.write("#PBS -e pbserr" + str(i) + "\n")
    submit_file.write("#PBS -l cput=1000:00:00,ncpus=1\n")
    submit_file.write(sys.argv[6] + " " + tmp_directory + tmp_filename + " " + str(i) + "\n")
    submit_file.close()
    os.system("qsub -q testqueue " + tmp_directory + "tmp_submit_script" + str(i) + ".sh")
    print ("Submitting Frame: " + str(i))
    

#now do cleanup
#os.system("rm -f " + tmp_filename)

