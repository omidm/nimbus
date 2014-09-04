
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
  

