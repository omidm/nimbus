# Building Nimbus

## Dependencies

Current version of Nimbus installs external dependencies such as `boost`,
`protobuf` and `leveldb`, as part of the build, automatically. Followings are a
small set of dependencies that you need to install manually.

#### Basic

    $ sudo apt-get install git
    $ sudo apt-get install gcc
    $ sudo apt-get install g++
    $ sudo apt-get install make

#### For plotting performance results

    $ sudo apt-get install python
    $ sudo apt-get install python-numpy
    $ sudo apt-get install python-matplotlib

#### For compiling Markdown docs

    $ sudo apt-get install pandoc


#### For creating straggler scenarios

    $ sudo apt-get install cgroup-bin

#### For deploying on EC2

    $ sudo apt install python-pip
    $ pip install -U boto

### For graphical simulations

    $ sudo apt-get install cmake-curses-gui
    $ sudo apt-get install freeglut3
    $ sudo apt-get install freeglut3-dev
    $ sudo apt-get install libqt4-opengl
    $ sudo apt-get install libqt4-opengl-dev
    $ sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
    $ sudo apt-get install libjpeg-dev
    $ sudo apt-get install libpng-dev
    $ sudo apt-get install libx11-6
    $ sudo apt-get install libx11-dev
    $ sudo apt-get install libxi-dev
    $ sudo apt-get install xdotool


## Building Nimbus Library

To build Nimbus core library, external dependencies, and default applications:

    $ make

## Building Graphical Simulations

For running graphical simulations with Nimbus you need to install PhysBAM
first. PhysBAM is officially supported with `gcc` version 4.5, and might have
problems with older versions. We provide scripts to install `gcc-4.5` on Ubuntu
12.04, 14.04, and 16.04. For more information please refer to documentation, at
[Installing GCC-4.5](installing-gcc.html). After installing `gcc-4.5`, first
clean the old Nimbus build and then rebuild everything with the new compiler.
From Nimbus root:

    $ make clean-hard
    $ make
    $ make physbam

It might take several minutes. For faster build use `-j <number-of-threads>`.
Then, run and display water and smoke simulations by issuing:

    $ make test-graphgics


## Instructions For Installing Dependencies Manually

It is highly recommended to follow the automatic Nimbus build that installs the
dependencies, but in case you want to install them manually, following
guidelines might be helpful.

### Ubuntu

#### Protocol Buffers: 

Download and install protocol buffer 2.6.0 from:
<http://code.google.com/p/protobuf/>. Installation guide:

* Before building the source, you may need to increase the following
  default values for message size:
    `kDefaultTotalBytesLimit`, and
    `kDefaultTotalBytesWarningThreshold`
    (both in `./src/google/protobuf/io/coded_stream.h`).

* Then, build and install as follows:

        $ sudo ./configure
        $ sudo make
        $ sudo make check
        $ sudo make install

#### Boost:

Download and install boost >1.55 from <http://www.boost.org>.
Installation guide:
    
* Use bootstrap with prefix option as follows:

        $ sudo ./bootstrap.sh --prefix="/usr/local"

* Then build and install as follows:

        $ sudo ./b2 install


#### Note:

1. You have to run the commands with `sudo` permissions so that the scripts
   could put the header files and libraries in the right directory.

2. If you already installed the `protobuf` through `apt-get`, there could be link
     errors. You need to remove the package:
     
        $ sudo apt-get remove --auto-remove libprotoc-dev
        $ sudo apt-get remove --auto-remove protobuf-compiler

3. You have to update the known libraries to the system:

        $ sudo ldconfig

4. To update the known libraries to system

        $ sudo ldconfig 

5. To install libraries and files using apt-get

        $ sudo apt-get install libprotoc-dev
        $ sudo apt-get install protobuf-compiler
        $ sudo apt-get install libboost1.48-all-dev

6. To uninstall the old files and libraries:

        $ sudo apt-get remove --auto-remove libprotoc-dev
        $ sudo apt-get remove --auto-remove protobuf-compiler
        $ sudo apt-get remove --auto-remove libboost1.48-dev


### MacOSX

#### Protocol Buffers: 

* Download and install Protocol Buffer from: 
  <http://code.google.com/p/protobuf/downloads/list>.
  Download `protobuf-2.5.0-zip` and follow the `README.txt` in the folder. 

#### Boost:

* Install boost with `Macports` (you may need to install  `MacPorts` first):

        $ sudo port install boost

#### Note:
1. For installing `brew` on El Captain run the command from homebrew website:

        $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

    but, `brew update` does not work directly, following is a fix:

        $ cd ~/.homebrew/Library
        $ git pull origin master
        $ brew update

