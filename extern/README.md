
# Installing gcc/g++ version 4.5 form source code

Although, Nimbus is supported with the latest g++ version, PhysBAM library
compiles only with g++ version 4.5. If you want to run fluid simulation
applications in the standalone mode or against Nimbus, you need to install g++
version 4.5, manually. The older versions of g++ are not shipped with Ubuntu
distributions including and after 12.04, and even mainstream package managers
doe not have the binaries for x86 architectures.

This tutorial gives you directions on how to install g++ version 4.5 from the
source code. It has been tested on Ubuntu 12.04, 14.04, and 16.04. For more
information visit gcc website: https://gcc.gnu.org/install/.

The following instruction are scripted in `install-g++-4.5.sh` file. 

## Prerequisites

### Install GMP
    $ sudo apt-get install libmpc-dev

Or you can install form the source. For more information visit:
https://gmplib.org/#DOWNLOAD

    $ sudo apt-get install m4
    $ tar -xjf gmp-6.1.1.tar.bz2
    $ cd gmp-6.1.1
    $ ./configure
    $ make
    $ make check
    $ sudo make install

### Install MPFR
    $ sudo apt-get install libmpfr-dev

Or you can install form the source. For more information visit:
http://www.mpfr.org/mpfr-current/mpfr.html

    $ tar -xvzf mpfr.3.1.4.tar.gz
    $ cd mpfr.3.1.4
    $ ./configure
    $ make
    $ make check
    $ sudo make install

### Install MPC
    $ sudo apt-get install libmpc-dev

Or you can install form the source. For more information visit:
http://www.multiprecision.org/index.php?prog=mpc

    $ tar -xvzf mpc-1.0.3.tar.gz
    $ cd mpc-1.0.3
    $ ./configure
    $ make
    $ make check
    $ sudo make install

### Install zip
    $ sudo apt-get install zip



## Install gcc 4.5 from source

It is straight forward to make and install from the source code:

    $ tar -xvzf gcc-4.5.3.tar.gz
    $ mkdir <objdir>
    $ cd <objdir>
    $ <srcdir>/configure --prefix=/usr/ --program-suffix=-4.5 --enable-languages=c,c++
    $ make
    $ sudo make install

by default the binaries  will be installed in `/usr/local` which is searched in
`$PATH` before `/usr/`. This effectively makes gcc/g++ 4.5 the default
compiler. It is better to change the install path as used above with the
`--prefix=/usr/` option, and then switch the compiler version with alternatives
as explained in the following.


If source code is compiled as given, there are few known issues. Here is a list
you can fix as follows:

### Fix the unknown `sys/cdefs.h` file

    $ sudo apt-get install libc6-dev-i386

### Fix the known `crti.o` file

    $ sudo ln -s /usr/lib/x86_64-linux-gnu /usr/lib64

### Fix unknown field  `siginfo'

In file `gcc-4.5.3/libgcc/../gcc/config/i386/linux-unwind.h` replace `struct siginfo` with `siginfo_t`:

    $ sed -i 's/struct siginfo/siginfo_t/g' gcc-4.5.3/libgcc/../gcc/config/i386/linux-unwind.h

### Fix `libs/libgcj.so: undefined reference to '__cxa_call_unexpected'`

If youwant to also install java compiler you need to update
`gcc-4.5.3/libjava/prims.cc` file as follows:

     #ifndef DISABLE_GETENV_PROPERTIES
    +#ifdef __GLIBC__
    +/* glibc 2.15+ provides even for C++ inline optimized ::isspace etc.
    +   Unfortunately those inlines are throw (), and call a function pointer
    +   (which is throw () too, but with -fnon-call-exceptions this results
    +   in a __cxa_call_unexpected call.  This macro disables the optimized
    +   version.  */
    +#define __NO_CTYPE 1
    +#endif
     #include <ctype.h>
     #include <java-props.h>
     #define PROCESS_GCJ_PROPERTIES process_gcj_properties()


## Change the gcc/g++ version with `alternatives`

### remove existing alternatives
    $ sudo update-alternatives --remove-all gcc 
    $ sudo update-alternatives --remove-all g++

### Install Alternatives
                         --install link name path priority

    $ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.5 10
    $ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.x 20

    $ sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.5 10
    $ sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.x 20

    $ sudo update-alternatives --install /usr/bin/cc cc /usr/bin/gcc 30
    $ sudo update-alternatives --set cc /usr/bin/gcc

    $ sudo update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++ 30
    $ sudo update-alternatives --set c++ /usr/bin/g++

#### Configure Alternatives
    $ sudo update-alternatives --config gcc
    $ sudo update-alternatives --config g++





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


## Note:

It is better to fix the header file of the shared library, upon creation. This
is done by using the `-install_name @rpath/$(LIB_NAME)` option in gcc/g++. For
example see how libnimbus is generated in `nimbus/src/Makefile`. You can also
use install name as `full-absolute-path` instead of `@rpath/...`.

