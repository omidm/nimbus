#!/bin/bash

# installs a version of protobuf library, header files, and executables locally.

TAR_FILE="mpc-1.0.3.tar.gz"
UNTAR_DIR="_temp_mpc_build"

rm -rf ${UNTAR_DIR} && mkdir -p ${UNTAR_DIR}

tar -xvzf ${TAR_FILE} -C ${UNTAR_DIR} --strip-components=1

cd ${UNTAR_DIR}
./configure
make -j $(nproc)
make check
sudo make install
cd -

rm -rf ${UNTAR_DIR}

