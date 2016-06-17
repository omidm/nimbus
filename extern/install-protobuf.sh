#!/bin/bash

# installs a version of protobuf library, header files, and executables locally.

TAR_FILE="protobuf-2.6.1.tar.gz"
UNTAR_DIR="_temp_protobuf_build"
INSTALL_DIR="protobuf"

rm -rf ${UNTAR_DIR} && mkdir -p ${UNTAR_DIR}
rm -rf ${INSTALL_DIR} && mkdir -p ${INSTALL_DIR}

cd ${INSTALL_DIR}
ABS_INSTALL_DIR=$(pwd)
cd -

tar -xf ${TAR_FILE} -C ${UNTAR_DIR} --strip-components=1

cd ${UNTAR_DIR}
./configure --prefix=${ABS_INSTALL_DIR} --exec-prefix=${ABS_INSTALL_DIR}
make -j 12
make check
make install
cd -

rm -rf ${UNTAR_DIR}

