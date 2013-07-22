#!/usr/bin/python
import sys
import os
import tests

result_directory=sys.argv[1]
visual_test=tests.VISUAL_TESTS("Rigid Bodies",sys.argv[1])

camera_directory=os.getcwd()+"/Scripts/test/cameras"
def rigid_bodies_test(testname,number):
    test=tests.TEST(result_directory=result_directory,
     name=testname,
     directory="Projects/rigid_bodies",
     command="./rigid_bodies_%s  -o Standard_Tests/%s %d"%(visual_test.platform,testname,number),
     output_directory="Standard_Tests/%s"%testname,
     views=[tests.VIEWER(name="std",viewer="../opengl_3d/opengl_3d_%s -w 800 -h 600"%visual_test.platform,keys="\"\"",
                         camera_script="%s/rigid_bodies_%d.camera"%(camera_directory,number))])

    visual_test.tests.append(test)

rigid_bodies_test("rigid_bodies_std_1_three_spheres",1)
rigid_bodies_test("rigid_bodies_std_2_three_pelvises",2)
rigid_bodies_test("rigid_bodies_std_3_billiard_balls",3)
rigid_bodies_test("rigid_bodies_std_4_spinning_ball_drop",4)
rigid_bodies_test("rigid_bodies_std_5_two_balls_slide_spin",5)
rigid_bodies_test("rigid_bodies_std_16_ether_test",16)
rigid_bodies_test("rigid_bodies_std_18_friction_force_propagation",18)
rigid_bodies_test("rigid_bodies_std_22_spinning_top",22)
rigid_bodies_test("rigid_bodies_std_23_spinning_cylinder",23)
rigid_bodies_test("rigid_bodies_std_24_diff_restitutions_test",24)
rigid_bodies_test("rigid_bodies_std_26_cluster_test",26)
rigid_bodies_test("rigid_bodies_std_27_cluster_test_with_kinematic_body",27)

visual_test.run()
