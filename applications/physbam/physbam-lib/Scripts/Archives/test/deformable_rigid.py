#!/usr/bin/python
import sys
import os
import tests

result_directory=sys.argv[1]
visual_test=tests.VISUAL_TESTS("Fluid 2D",sys.argv[1])

camera_directory=os.getcwd()+"/Scripts/test/cameras"
def deformable_rigid_test(testname,number):
    test=tests.TEST(result_directory=result_directory,
     name=testname,
     directory="Projects/deformable_rigid_coupling",
     command="./deformable_rigid_coupling_%s  -o Standard_Tests/%s %d"%(visual_test.platform,testname,number),
     output_directory="Standard_Tests/%s"%testname,
     views=[tests.VIEWER(name="std",viewer="../opengl_3d/opengl_3d_%s -w 800 -h 600"%visual_test.platform,keys="\"\"",
                         camera_script="%s/deformable_rigid_%d.camera"%(camera_directory,number))])

    visual_test.tests.append(test)

deformable_rigid_test("deformable_rigid_std_2_block_stack",2)
deformable_rigid_test("deformable_rigid_std_5_point_joint_2_blocks_with_deformable",5)
deformable_rigid_test("deformable_rigid_std_6_point_joint_2_blocks_one_suspended_with_deformable",6)
deformable_rigid_test("deformable_rigid_std_7_torsion_spring_with_kinematic_and_arb",7)

visual_test.run()
