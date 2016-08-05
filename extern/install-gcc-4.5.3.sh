#!/bin/bash

# installs gcc version 4.5 and required dependencies. See the README.md file
# for more information.

### Install GMP
sudo apt-get install --yes libgmp-dev
# ./install-gmp-6.1.1.sh

### Install MPFR
sudo apt-get install --yes libmpfr-dev
# ./install-mpfr-3.1.4.sh

### Install MPC
sudo apt-get install --yes libmpc-dev
# ./install-mpc-1.0.3.sh

### Install zip
sudo apt-get install --yes zip

### Update search paths:
echo 'CPATH=/usr/include/x86_64-linux-gnu:$CPATH' >> ~/.profile
echo 'CPATH=/usr/include/x86_64-linux-gnu:$CPATH' >> ~/.bash_profile
source ~/.profile

echo 'LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LIBRARY_PATH' >> ~/.profile
echo 'LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LIBRARY_PATH' >> ~/.bash_profile
source ~/.profile


TAR_FILE="gcc-4.5.3.tar.gz"
SOURCE_DIR="_temp_gcc_source"
OBJ_DIR="_temp_gcc_obj"

rm -rf ${SOURCE_DIR} && mkdir -p ${SOURCE_DIR}
rm -rf ${OBJ_DIR} && mkdir -p ${OBJ_DIR}

cd ${SOURCE_DIR}
ABS_SOURCE_DIR=$(pwd)
cd -

### Install gcc 4.5 from source
### configure
tar -xvzf ${TAR_FILE} -C ${SOURCE_DIR} --strip-components=1
cd ${OBJ_DIR}
${ABS_SOURCE_DIR}/configure --prefix=/usr/ --program-suffix=-4.5 --enable-languages=c,c++

### fix minor problems with the source file
sed -i 's/struct siginfo/siginfo_t/g' ${ABS_SOURCE_DIR}/libgcc/../gcc/config/i386/linux-unwind.h
sudo apt-get install --yes libc6-dev-i386
sudo ln -s /usr/lib/x86_64-linux-gnu /usr/lib64

### build
make -j $(nproc)
sudo make install
cd -

rm -rf ${SOURCE_DIR}
rm -rf ${OBJ_DIR}

### switch gcc/g++ version
if ( ls /usr/bin/gcc-4.5 &> /dev/null ); then
  echo "SUCCESS: gcc-4.5 was installed, setting up the symlink"
  sudo rm /usr/bin/gcc
  sudo ln -s /usr/bin/gcc-4.5 /usr/bin/gcc
else
  echo "ERROR: gcc-4.5 was NOT installed!"
fi

if ( ls /usr/bin/g++-4.5 &> /dev/null ); then
  echo "SUCCESS: g++-4.5 was installed, setting up the symlink"
  sudo rm /usr/bin/g++
  sudo ln -s /usr/bin/g++-4.5 /usr/bin/g++
else
  echo "ERROR: g++-4.5 was NOT installed!"
fi

