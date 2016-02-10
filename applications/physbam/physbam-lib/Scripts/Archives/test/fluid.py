#!/usr/bin/python
import sys

import tests

result_directory=sys.argv[1]
visual_test=tests.VISUAL_TESTS("Fluid 2D",sys.argv[1])


def fluid2d_test(testname,number,resolution):
    test=tests.TEST(result_directory=result_directory,
     name=testname,
     directory="Projects/fluids_2d",
     command="./fluids_2d_%s -resolution %d -o Standard_Tests/%s %d"%(visual_test.platform,resolution,testname,number),
     output_directory="Standard_Tests/%s"%testname,
     views=[tests.VIEWER(name="std",viewer="../opengl_2d/opengl_2d_%s -w 800 -h 600"%visual_test.platform,keys="\;\;6",camera_script="water_falling_drop_2d.camera")])

    visual_test.tests.append(test)
fluid2d_test("fluid_2d_std_1_still",1,5)
fluid2d_test("fluid_2d_std_2_falling_drop",2,5)
fluid2d_test("fluid_2d_std_3_two_source",3,5)
fluid2d_test("fluid_2d_std_4_splash",4,5)
fluid2d_test("fluid_2d_std_12_viscosity",12,5)
fluid2d_test("fluid_2d_std_13_variable_viscosity",13,5)

visual_test.run()
