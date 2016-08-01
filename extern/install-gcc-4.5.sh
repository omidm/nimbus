#!/bin/bash

# installs gcc version 4.5 and required dependencies. See the README.md file
# for moreinformation.

TAR_FILE="gmp-6.1.1.tar.bz2"
UNTAR_DIR="_temp_gmp_build"

rm -rf ${UNTAR_DIR} && mkdir -p ${UNTAR_DIR}

### Install GMP
sudo apt-get install --yes m4
tar -xjf ${TAR_FILE} -C ${UNTAR_DIR} --strip-components=1
cd ${UNTAR_DIR}
./configure
make -j 12
make check
sudo make install
cd -

rm -rf ${UNTAR_DIR}


### Install MPFR
sudo apt-get install --yes libmpfr-dev

### Install MPC
sudo apt-get install --yes libmpc-dev

### Install zip
sudo apt-get install --yes zip



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
make -j 12
sudo make install
cd -

### switch gcc/g++ version
sudo rm /usr/bin/gcc
sudo ln -s /usr/bin/gcc-4.5 /usr/bin/gcc
sudo rm /usr/bin/g++
sudo ln -s /usr/bin/g++-4.5 /usr/bin/g++


rm -rf ${SOURCE_DIR}
rm -rf ${OBJ_DIR}

