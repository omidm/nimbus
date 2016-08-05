#!/bin/bash

# installs required dependencies for building and running physbam applications.

### cmake and ccmake
sudo apt-get install --yes cmake-curses-gui

### opengl
sudo apt-get install --yes freeglut3
sudo apt-get install --yes freeglut3-dev
sudo apt-get install --yes libqt4-opengl
sudo apt-get install --yes libqt4-opengl-dev

### openmpi
sudo apt-get install --yes openmpi-bin openmpi-doc libopenmpi-dev

### jpeg, png
sudo apt-get install --yes libjpeg-dev
sudo apt-get install --yes libpng-dev
 
### x11

sudo apt-get install --yes libx11-6
sudo apt-get install --yes libx11-dev
sudo apt-get install --yes libxi-dev

