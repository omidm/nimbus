
Nimbus Installation Dependencies - Linux
=========================================
*** Protocol Buffers: 
Download and install protocol buffer 2.6.0 from
    code.google.com/p/protobuf/

Note:
you may need to increase the following default values for message size:
    kDefaultTotalBytesLimit
    kDefaultTotalBytesWarningThreshold
both in ./src/google/protobuf/io/coded_stream.h before building the source.

*** Boost:
Download and install boost >1.55 from
    www.boost.org.

Note:
use bootstrap with prefix option as follows:
    sudo ./bootsrap.sh --prefix="/usr/local"


--------------------------------------------------
To update the known libraries to system
--------------------------------------------------
sudo ldconfig 


--------------------------------------------------
To install libraries and files using apt-get
--------------------------------------------------
sudo apt-get install libprotoc-dev
sudo apt-get install protobuf-compiler
sudo apt-get install libboost1.48-all-dev

--------------------------------------------------
To uninstall the old files and libraries:
--------------------------------------------------
sudo apt-get remove --auto-remove libprotoc-dev
sudo apt-get remove --auto-remove protobuf-compiler
sudo apt-get remove --auto-remove libboost1.48-dev




Nimbus Installation Dependencies - MacOSX
=========================================

** Install Protocol Buffer:
  download protobuf-2.5.0-zip from
  "http://code.google.com/p/protobuf/downloads/list"
  follow the README.txt in the folder. 

** Install Boost:
  sudo port install boost
  (you may need to install  MacPorts first)
  


How to Remove DBG in compile time:

1. uncomment -D_NIMBUS_NO_DBG from Makeinclude

2. uncomment add_definitions(-D_NIMBUS_NO_DBG)
   from application/water_multiple/CMakeLists.txt

3. comment out ADD_DEFINITIONS (-Werror)
   from application/physbam-lib/Scripts/CMake/Compiler_Flags.cmake

4. make sure that you make clean-hard for ninmus library,
   and make clean for rest.
