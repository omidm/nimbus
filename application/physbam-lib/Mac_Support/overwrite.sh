#!/bin/bash

cp main.cpp ../../physbam-app/opengl_2d/main.cpp
cp CMakeLists.txt ../Public_Library/CMakeLists.txt
cp OPENGL_AGL_PBUFFER.cpp ../Public_Library/PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_AGL_PBUFFER.cpp
cp PHOTON_MAP.cpp ../Public_Library/PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/PHOTON_MAP.cpp
cp THREADED_RIGIDS.h ../Public_Library/PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/THREADED_RIGIDS.h
cp Is_NaN.h ../Public_Library/PhysBAM_Tools/Math_Tools/Is_NaN.h
cp Public_Library.cmake ../Scripts/CMake/Public_Library.cmake
