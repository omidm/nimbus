create_lib_makefile.sh -f Makefile.PhysBAM PhysBAM `ls -1 *.cpp */*.cpp | grep -v OpenGL | grep -v OpenGL_Components`
create_lib_makefile.sh -f Makefile.PhysBAM_OpenGL PhysBAM_OpenGL `ls -1 OpenGL/*.cpp OpenGL_Components/*.cpp`
