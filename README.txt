
=============================================================================
Nimbus Installation Dependencies - Linux
=============================================================================

--------------------------------------------------
*** Protocol Buffers: 
--------------------------------------------------

  * Download and install protocol buffer 2.6.0 from:
    code.google.com/p/protobuf/
  
  * Installation guide:
    - Before building the source, you may need to increase the following
      default values for message size:
        kDefaultTotalBytesLimit
        kDefaultTotalBytesWarningThreshold
        (both in ./src/google/protobuf/io/coded_stream.h)
    
    - Then, build and install as follows:
        $ sudo ./configure
        $ sudo make
        $ sudo make check
        $ sudo make install


--------------------------------------------------
*** Boost:
--------------------------------------------------

  * Download and install boost >1.55 from
    www.boost.org.
  
  * Installation guide:
    - Use bootstrap with prefix option as follows:
        $ sudo ./bootstrap.sh --prefix="/usr/local"

    - Then build and install as follows:
        $ sudo ./b2 install


--------------------------------------------------
** Note:
--------------------------------------------------

  1. You have to run the commands with "sudo" permissions so that the scripts
     could put the header files and libraries in the right directory.

  2. You have to update the known libraries to the system:
      $ sudo ldconfig

  3. For PhysBAM related applications you need old compiler g++-4.5.
     However, Ubuntu 12.04's default g++-4.6 works well too. 
     
     ** EC2 base Ubuntu 12.04 AMI: ami-fa9cf1ca


--------------------------------------------------
To update the known libraries to system
--------------------------------------------------
$ sudo ldconfig 


--------------------------------------------------
To install libraries and files using apt-get
--------------------------------------------------
$ sudo apt-get install libprotoc-dev
$ sudo apt-get install protobuf-compiler
$ sudo apt-get install libboost1.48-all-dev


--------------------------------------------------
To uninstall the old files and libraries:
--------------------------------------------------
$ sudo apt-get remove --auto-remove libprotoc-dev
$ sudo apt-get remove --auto-remove protobuf-compiler
$ sudo apt-get remove --auto-remove libboost1.48-dev



=============================================================================
Nimbus Installation Dependencies - MacOSX
=============================================================================

** Install Protocol Buffer:
  download protobuf-2.5.0-zip from
  "http://code.google.com/p/protobuf/downloads/list"
  follow the README.txt in the folder. 

** Install Boost:
  sudo port install boost
  (you may need to install  MacPorts first)
  


=============================================================================
How to remove DBG in compile time:
=============================================================================

1. uncomment -D_NIMBUS_NO_DBG from Makeinclude

2. uncomment add_definitions(-D_NIMBUS_NO_DBG)
   from application/water_multiple/CMakeLists.txt

3. comment out ADD_DEFINITIONS (-Werror)
   from application/physbam-lib/Scripts/CMake/Compiler_Flags.cmake

4. make sure that you make clean-hard for ninmus library,
   and make clean for rest.



=============================================================================
How to profile code using perf:
=============================================================================

1. Compile the code using "-fno-omit-frame-pointer -ggdb -pg" flags.
   a. Uncomment the line in Makeinclude that is for profiling with perf.
   b. make clean-hard and make.

2. Run "perf record -g <executable> <args>". It dumps per.data file.

3. Run "perf report -g fractal,0.5,caller -i <dumped-file>". It gives you the
   call stacks based on the caller trace.

4. To produce FlameGraph, download the repo from github. run the following
   command in the directory that has perf.data:
   "perf script | ./stackcollapse-perf.pl | ./flamegraph.pl > result.svg"



