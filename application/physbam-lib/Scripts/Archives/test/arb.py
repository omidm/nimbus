#!/usr/bin/python
import sys
import os
import tests

result_directory=sys.argv[1]
visual_test=tests.VISUAL_TESTS("Articulated Rigid Bodies",sys.argv[1])

camera_directory=os.getcwd()+"/Scripts/test/cameras"
def articulated_rigid_bodies_test(testname,number):
    test=tests.TEST(result_directory=result_directory,
     name=testname,
     directory="Projects/articulated_rigid_bodies",
     command="./articulated_rigid_bodies_%s  -o Standard_Tests/%s %d"%(visual_test.platform,testname,number),
     output_directory="Standard_Tests/%s"%testname,
     views=[tests.VIEWER(name="std",viewer="../opengl_3d/opengl_3d_%s -w 800 -h 600"%visual_test.platform,keys="\"\"",
                         camera_script="%s/articulated_rigid_bodies_%d.camera"%(camera_directory,number))])

    visual_test.tests.append(test)

articulated_rigid_bodies_test("articulated_rigid_bodies_std_1_point_joint",1)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_2_rigid_joint",2)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_3_hinge_joint",3)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_4_twist_joint",4)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_5_string_of_rigid_joints",5)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_6_closed_loop_nonconvex_objects",6)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_7_pd_planks_form_O",7)
# articulated_rigid_bodies_test("articulated_rigid_bodies_std_8_8_block_cluster",8)
# articulated_rigid_bodies_test("articulated_rigid_bodies_std_9_8_block_cluster_breaking",9)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_10_prismatic_joint",10)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_11_constrained_hinge",11)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_12_heavy_bottom_link",12)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_13_universal_joint",13)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_14_pre_stab_test",14)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_15_rigid_joint_using_point_joints",15)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_16_twist_joint_using_point_joints",16)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_20_post_stab",20)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_27_normal_joint",27)
articulated_rigid_bodies_test("articulated_rigid_bodies_std_28_angle_joint_with_kinematic_obj",28)

visual_test.run()

1-16, 20, 27, 28.
