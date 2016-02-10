#!/usr/bin/python
import sys
import os
import tests

result_directory=sys.argv[1]
visual_test=tests.VISUAL_TESTS("Fluid 2D",sys.argv[1])

camera_directory=os.getcwd()+"/Scripts/test/cameras"
def solids_3d_test(testname,number):
    test=tests.TEST(result_directory=result_directory,
     name=testname,
     directory="Projects/solids_3d",
     command="./solids_3d_%s  -o Standard_Tests/%s %d"%(visual_test.platform,testname,number),
     output_directory="Standard_Tests/%s"%testname,
     views=[tests.VIEWER(name="std",viewer="../opengl_3d/opengl_3d_%s -w 800 -h 600"%visual_test.platform,keys="\"\"",
                         camera_script="%s/solids_3d_%d.camera"%(camera_directory,number))])

    visual_test.tests.append(test)

solids_3d_test("solids_3d_std_1_sphere",1)
solids_3d_test("solids_3d_std_2_embedded_torus",2)
solids_3d_test("solids_3d_std_5_curtain_and_ball",5)

visual_test.run()
