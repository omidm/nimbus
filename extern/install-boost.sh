#!/bin/bash

# installs a version of boost library and header files locally, it does not
# install the entire binaries, only those required by nimbus.

TAR_FILE="boost_1_61_0.tar.bz2"
UNTAR_DIR="_temp_boost_build"
INSTALL_DIR="boost"
LIBS_TO_INSTALL="thread,system,program_options"

rm -rf ${UNTAR_DIR} && mkdir -p ${UNTAR_DIR}
rm -rf ${INSTALL_DIR} && mkdir -p ${INSTALL_DIR}

cd ${INSTALL_DIR}
ABS_INSTALL_DIR=$(pwd)
cd -

tar --bzip2 -xf ${TAR_FILE} -C ${UNTAR_DIR} --strip-components=1

cd ${UNTAR_DIR}
./bootstrap.sh --prefix=${ABS_INSTALL_DIR} --with-libraries=${LIBS_TO_INSTALL}
./b2 install
cd -

rm -rf ${UNTAR_DIR}

