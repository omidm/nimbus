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
echo 'export CPATH=/usr/include/x86_64-linux-gnu' >> ~/.bashrc
echo 'export CPATH=/usr/include/x86_64-linux-gnu' >> ~/.profile
echo 'export CPATH=/usr/include/x86_64-linux-gnu' >> ~/.bash_profile
echo 'export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu' >> ~/.bashrc
echo 'export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu' >> ~/.profile
echo 'export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu' >> ~/.bash_profile
source ~/.profile


TAR_FILE="gcc-4.5.3.tar.gz"
SOURCE_DIR="_temp_gcc_source"
OBJ_DIR="_temp_gcc_obj"
INSTALL_DIR="gcc-4.5"

rm -rf ${SOURCE_DIR} && mkdir -p ${SOURCE_DIR}
rm -rf ${OBJ_DIR} && mkdir -p ${OBJ_DIR}

# get the absolute path to the directory of the script
ABS="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

### Install gcc 4.5 from source
### configure
tar -xvzf ${TAR_FILE} -C ${SOURCE_DIR} --strip-components=1
cd ${OBJ_DIR}
${ABS}/${SOURCE_DIR}/configure --prefix=${ABS}/${INSTALL_DIR} --program-suffix=-4.5 --enable-languages=c,c++ MAKEINFO=missing

### fix minor problems with the source file
sed -i 's/struct siginfo/siginfo_t/g' ${ABS}/${SOURCE_DIR}/libgcc/../gcc/config/i386/linux-unwind.h
sudo apt-get install --yes libc6-dev-i386
sudo ln -s /usr/lib/x86_64-linux-gnu /usr/lib64

### build
make -j $(nproc)
sudo make install
cd -

rm -rf ${SOURCE_DIR}
rm -rf ${OBJ_DIR}

### give directives for switching gcc/g++ versions
if ( ls ${ABS}/${INSTALL_DIR}/bin/gcc-4.5 &> /dev/null ); then
  echo "SUCCESS: gcc-4.5 was installed at ${ABS}/${INSTALL_DIR}/bin/gcc-4.5"
  echo "*****************************************************************************"
  echo "**** NOTE: you need to either open a new shell or execute 'source ~/.profile'"
  echo "**** NOTE: use update-alternatives to switch between the compiler versions."
  echo "*****************************************************************************"
  # It is better to use the update-alternatives to switch the compilers instead of symlink.
  # sudo rm /usr/bin/gcc
  # sudo ln -s /usr/bin/gcc-4.5 /usr/bin/gcc
else
  echo "ERROR: gcc-4.5 was NOT installed!"
fi

if ( ls ${ABS}/${INSTALL_DIR}/bin/g++-4.5 &> /dev/null ); then
  echo "SUCCESS: g++-4.5 was installed at ${ABS}/${INSTALL_DIR}/bin/g++-4.5"
  echo "*****************************************************************************"
  echo "**** NOTE: you need to either open a new shell or execute 'source ~/.profile'"
  echo "**** NOTE: use update-alternatives to switch between the compiler versions."
  echo "*****************************************************************************"
  # It is better to use the update-alternatives to switch the compilers instead of symlink.
  # sudo rm /usr/bin/g++
  # sudo ln -s /usr/bin/g++-4.5 /usr/bin/g++
else
  echo "ERROR: g++-4.5 was NOT installed!"
fi

