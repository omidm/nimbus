#!/bin/bash

# installs a version of protobuf library, header files, and executables locally.

TAR_FILE="gmp-6.1.1.tar.bz2"
UNTAR_DIR="_temp_gmp_build"

rm -rf ${UNTAR_DIR} && mkdir -p ${UNTAR_DIR}

tar -xjf ${TAR_FILE} -C ${UNTAR_DIR} --strip-components=1

sudo apt-get install --yes m4

cd ${UNTAR_DIR}
./configure
make -j 16
make check
sudo make install
cd -

rm -rf ${UNTAR_DIR}

