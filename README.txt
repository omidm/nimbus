
Nimbus Installation Dependencies - Linux
=========================================

sudo apt-get install libprotoc-dev
sudo apt-get install protobuf-compiler
sudo apt-get install libboost1.48-all-dev


----------------------------------------
To uninstall the old boost:
----------------------------------------
sudo apt-get autoremove libboost1.48-dev

----------------------------------------
To install the latest version:
----------------------------------------
* Download the latest copy of boost
tar --bzip2 -xf /path/to/boost_1_xx_0.tar.bz2
sudo ./bootsrap.sh --prefix="/usr/local"
sudo ./b2

----------------------------------------
To update the known libraries
----------------------------------------
sudo ldconfig


Nimbus Installation Dependencies - MacOSX
=========================================

** Install Protocol Buffer:
  download protobuf-2.5.0-zip from
  "http://code.google.com/p/protobuf/downloads/list"
  follow the README.txt in the folder. 

** Install Boost:
  sudo port install boost
  (you may need to install  MacPorts first)
  

