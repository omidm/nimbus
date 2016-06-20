# Installing shared libraries and linker issues in MacOSX 

For further reading see:
https://blogs.oracle.com/dipol/entry/dynamic_libraries_rpath_and_mac
http://stackoverflow.com/questions/33665781/dependencies-on-boost-library-dont-have-full-path
http://stackoverflow.com/questions/4513799/how-to-set-the-runtime-path-rpath-of-an-executable-with-gcc-under-mac-osx


When installing shared libraries and handling linker calls in MacOSX, there are
few differences compared to Linux. In both systems, the linker looks for the
shared libraries in some default paths, and also you can add extra runtime
paths by setting environment variables. For example, after installing the
shared libraries in this directory, you can launch executables that need this
libraries by setting the `DYLD_LIBRARY_PATH` in MacOSX (and `LD_LIBRARY_PATH` in Linux) as follows

    $ cd nimbus/nodes/nimbus_controller/
    $ DYLD_LIBRARY_PATH="../../extern/boost/lib/:../../extern/protobuf/lib/" ./nimbus_controller -p 5800 -w 1

However, this is annoying and you want to add the linker path to the binary,
for example, the `-Wl,-rpath <path-to-lib>` directives in gcc. However, this is
not just enough, in MacOSX. In MacOSX, when linking against a shared library,
it has a header file that enforces where to find that library. You can find
that using the diagnostic tool as follows:

    otool -L <executable or shared-library>

If you check for the `nimbus_controller` executable, you will see that the
`protobuf` library has absolute path, while the `boost` libraries are relative.
That is why linker complains only about not finding `boost`. To fix this issue,
you have two options you need to replace the path name of the boost libraries
after it is generayed with `install_name_tool`. There are two ways:

### First option:

You can first generate your executable/libnimbus.so and then fix the
library paths in the generated executable. Without changing the install name of
boost libraries, make nimbus library and nimbus_controller. If you want to do
the inspection, e.g.:

    $ otool -L nimbus_controller

you will find out that, protobuf libraries have absolute path while boost paths
are relative. Now you can fix them:

    # FOR INSPECTION: otool -L nimbus_controller
    $ install_name_tool -change libboost_thread.dylib @executable_path/../../extern/boost/lib/libboost_thread.dylib nimbus_controller
    $ install_name_tool -change libboost_system.dylib @executable_path/../../extern/boost/lib/libboost_system.dylib nimbus_controller
    $ install_name_tool -change libboost_program_options.dylib @executable_path/../../extern/boost/lib/libboost_program_options.dylib nimbus_controller

    # FOR INSPECTION: otool -L libnimbus.so
    $ install_name_tool -change libboost_thread.dylib @executable_path/../../extern/boost/lib/libboost_thread.dylib libnimbus.so
    $ install_name_tool -change libboost_system.dylib @executable_path/../../extern/boost/lib/libboost_system.dylib libnimbus.so
    $ install_name_tool -change libboost_program_options.dylib @executable_path/../../extern/boost/lib/libboost_program_options.dylib libnimbus.so

Also, the boost itself, needs to be updated:

    $ FOR INSPECTION: otool -L libboost_thread.dylib
    $ install_name_tool -change libboost_system.dylib @executable_path/../../extern/boost/lib/libboost_system.dylib libboost_thread.dylib

Then you can run directly:

    $ ./nimbus_controller -p 5800 -w 1



### Second option (RECOMMENDED):

Directly fix the header file of the boost libraries, and then make your
executables and other shared libraries:

    $ install_name_tool libboost_thread.dylib -id @rpath/libboost_thread.dylib
    $ install_name_tool -change libboost_system.dylib @rpath/libboost_system.dylib libboost_thread.dylib
    $ install_name_tool libboost_system.dylib -id @rpath/libboost_system.dylib
    $ install_name_tool libboost_program_options.dylib -id @rpath/libboost_program_options.dylib


Then, remake the nimbus_controller and libnimbus.so.

**IMPORTANT**: you could also replace the `@rpath` with `$(pwd)` to give the
full-absolute-path, this is a beter solution.


# Note:

It is better to fix the header file of the shared library, upon creation. This
is done by using the `-install_name @rpath/$(LIB_NAME)` option in gcc/g++. For
example see how libnimbus is generated in `nimbus/src/Makefile`. You can also
use install name as `full-absolute-path` instead of `@rpath/...`.

