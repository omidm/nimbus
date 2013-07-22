######################################################################
# Copyright 2011, Rahul Sheth.
# This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
######################################################################

#basics
#previous versions of CMake do not support ExternalProject or ProcessorCount
CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
## allow more human readable "if then else" constructs
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
## Use verbose make file
#SET (CMAKE_VERBOSE_MAKEFILE TRUE)

#You should set CC and CXX environment variables for CMake 
#to find a compiler on UNIX systems

#Set a release build type if nothing is defined
IF(NOT WIN32)
    #Set a build type based on the folder name if nothing is defined
    IF(NOT CMAKE_BUILD_TYPE)
        GET_FILENAME_COMPONENT(PHYSBAM_BUILD_FOLDER_NAME ${CMAKE_CURRENT_BINARY_DIR} NAME)
        SET(CMAKE_BUILD_TYPE ${PHYSBAM_BUILD_FOLDER_NAME} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
    ENDIF(NOT CMAKE_BUILD_TYPE)
ENDIF(NOT WIN32)

#project name & version
PROJECT(PhysBAM)

SET(SCRIPTS_PATH ${PROJECT_SOURCE_DIR}/../Scripts)
#add custom find modules to CMake
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SCRIPTS_PATH}/CMake/Modules)

SET(${PROJECT_NAME}_MAJOR_VERSION 0)
SET(${PROJECT_NAME}_MINOR_VERSION 1)
SET(${PROJECT_NAME}_PATCH_LEVEL 0)

######################################################
# Build Options and Compiler Flags                   #
######################################################
INCLUDE("${SCRIPTS_PATH}/CMake/Build_Options.cmake")
INCLUDE("${SCRIPTS_PATH}/CMake/Compiler_Flags.cmake")

######################################################
# Source Group Creation - Important for IDEs         #
######################################################
INCLUDE("${SCRIPTS_PATH}/CMake/Source_Group.cmake")

######################################################
# Dependencies                                       #
######################################################
INCLUDE(ExternalProject)
SET(EP_PREFIX ${PROJECT_SOURCE_DIR}/../External_Libraries)
# zlib
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ExternalProject_Add(zlib 
        PREFIX ${EP_PREFIX}
        PATCH_COMMAND ${CMAKE_COMMAND} -E remove <SOURCE_DIR>/zconf.h
        CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
    SET(PHYSBAM_DEPENDENCIES zlib)
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
IF (WIN32)
    SET(ZLIB_LIBRARY debug ${EP_PREFIX}/lib/zlibd.lib optimized ${EP_PREFIX}/lib/zlib.lib)
ELSE (WIN32)
    SET(ZLIB_LIBRARY debug ${EP_PREFIX}/lib/libzd.so optimized ${EP_PREFIX}/lib/libz.so)
ENDIF (WIN32)

SET(PHYSBAM_LINK_LIBRARIES ${ZLIB_LIBRARY})

#PNG
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ExternalProject_Add(libpng 
        PREFIX ${EP_PREFIX}
        CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
    ADD_DEPENDENCIES(libpng zlib)
    IF(PHYSBAM_USE_LIBPNG)
        LIST(APPEND PHYSBAM_DEPENDENCIES libpng)
    ENDIF(PHYSBAM_USE_LIBPNG)
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
IF (WIN32)
    SET(PNG_LIBRARY debug ${EP_PREFIX}/lib/libpng15d.lib optimized ${EP_PREFIX}/lib/libpng15.lib)
ELSE (WIN32)
    SET(PNG_LIBRARY debug ${EP_PREFIX}/lib/libpng15d.so optimized ${EP_PREFIX}/lib/libpng15.so)
ENDIF (WIN32)

IF(PHYSBAM_USE_LIBPNG)
    LIST(APPEND PHYSBAM_LINK_LIBRARIES ${PNG_LIBRARY})
ENDIF(PHYSBAM_USE_LIBPNG)

#JPEG
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ExternalProject_Add(libjpeg 
        PREFIX ${EP_PREFIX}
        CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
    #Artificial dependency added to prevent building libraries simultaneously
    ADD_DEPENDENCIES(libjpeg libpng)
    IF(PHYSBAM_USE_LIBJPEG)
        LIST(APPEND PHYSBAM_DEPENDENCIES libjpeg)
    ENDIF(PHYSBAM_USE_LIBJPEG)
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
IF (WIN32)
    SET(JPEG_LIBRARY debug ${EP_PREFIX}/lib/jpegd.lib optimized ${EP_PREFIX}/lib/jpeg.lib)
ELSE (WIN32)
    SET(JPEG_LIBRARY debug ${EP_PREFIX}/lib/libjpegd.a optimized ${EP_PREFIX}/lib/libjpeg.a)
ENDIF (WIN32)

IF(PHYSBAM_USE_LIBJPEG)
    LIST(APPEND PHYSBAM_LINK_LIBRARIES ${JPEG_LIBRARY})
ENDIF(PHYSBAM_USE_LIBJPEG)

#freeglut and OpenGL
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ExternalProject_Add(freeglut 
        PREFIX ${EP_PREFIX}
        CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
    #Artificial dependency added to prevent building libraries simultaneously
    ADD_DEPENDENCIES(freeglut libjpeg)
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
IF (WIN32)
    SET(GLUT_LIBRARY debug ${EP_PREFIX}/lib/freeglut.lib optimized ${EP_PREFIX}/lib/freeglut.lib)
    FILE(COPY ${EP_PREFIX}/src/glheaders/wglext.h DESTINATION ${EP_PREFIX}/include)
ELSE (WIN32)
    SET(GLUT_LIBRARY debug ${EP_PREFIX}/lib/libfreeglut.so optimized ${EP_PREFIX}/lib/libfreeglut.so)
    FILE(COPY ${EP_PREFIX}/src/glheaders/glxext.h DESTINATION ${EP_PREFIX}/include)
ENDIF (WIN32)
FIND_PACKAGE(OpenGL REQUIRED)

#pthreads
IF (WIN32)
    IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
        ExternalProject_Add(pthreads-win32 
            PREFIX ${EP_PREFIX})
        #Artificial dependency added to prevent building libraries simultaneously
        #this dependency is Windows-only, so don't use it for other cross-platform dependencies
        ADD_DEPENDENCIES(pthreads-win32 freeglut)
        IF(PHYSBAM_USE_PTHREADS)
            LIST(APPEND PHYSBAM_DEPENDENCIES pthreads-win32)
        ENDIF(PHYSBAM_USE_PTHREADS)
    ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    SET(PTHREADS_LIBRARY debug ${EP_PREFIX}/lib/pthreadVC2d.lib optimized ${EP_PREFIX}/lib/pthreadVC2.lib)
ELSE (WIN32)
    IF(PHYSBAM_USE_PTHREADS)
        FIND_PACKAGE(Threads)
        SET(PTHREADS_LIBRARY ${CMAKE_THREAD_LIBS_INIT})
    ENDIF(PHYSBAM_USE_PTHREADS)
ENDIF (WIN32)

IF(PHYSBAM_USE_PTHREADS)
    LIST(APPEND PHYSBAM_LINK_LIBRARIES ${PTHREADS_LIBRARY})
ENDIF(PHYSBAM_USE_PTHREADS)

#fltk - if you move jpeg, png, zlib, or pthreads, edit the fltk CMake files
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ExternalProject_Add(fltk 
        PREFIX ${EP_PREFIX}
        CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
    #if pthreads-win32 target doesn't exist, it will simply be ignored
    ADD_DEPENDENCIES(fltk pthreads-win32 libjpeg libpng zlib)
    #Artificial dependency added to prevent building libraries simultaneously
    ADD_DEPENDENCIES(fltk freeglut)
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
IF (WIN32)
    SET(FLTK_LIBRARY debug ${EP_PREFIX}/lib/fltkd.lib ${EP_PREFIX}/lib/fltkformsd.lib ${EP_PREFIX}/lib/fltkgld.lib ${EP_PREFIX}/lib/fltkimagesd.lib
        optimized ${EP_PREFIX}/lib/fltk.lib ${EP_PREFIX}/lib/fltkforms.lib ${EP_PREFIX}/lib/fltkgl.lib ${EP_PREFIX}/lib/fltkimages.lib)
    SET(FLTK_INCLUDE_DIR ${EP_PREFIX}/include)
ELSE(WIN32)
ENDIF(WIN32)

#StackWalker - required for Windows backtrace/debugging support
IF(WIN32)
    IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
        ExternalProject_Add(stackwalker 
            PREFIX ${EP_PREFIX}
            CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER})
        LIST(APPEND PHYSBAM_DEPENDENCIES stackwalker)
        #Artificial dependency added to prevent building libraries simultaneously 
        #this dependency is Windows-only, so don't use it for other cross-platform dependencies
        ADD_DEPENDENCIES(stackwalker fltk)
    ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    SET(STACKWALKER_LIBRARY debug ${EP_PREFIX}/lib/stackwalkerd.lib optimized ${EP_PREFIX}/lib/stackwalker.lib)
    LIST(APPEND PHYSBAM_LINK_LIBRARIES ${STACKWALKER_LIBRARY})
ENDIF (WIN32)

#OpenMPI
#build on Windows only, assume installed on Linux
IF (WIN32)
    IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
        ExternalProject_Add(openmpi 
            PREFIX ${EP_PREFIX}
            URL ${EP_PREFIX}/src/openmpi.tar.gz
            CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} 
                -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_INSTALL_PREFIX=${EP_PREFIX})
        #Artificial dependency added to prevent building libraries simultaneously
        #this dependency is Windows-only, so don't use it for other cross-platform dependencies
        ADD_DEPENDENCIES(openmpi stackwalker)
    ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    SET(OPENMPI_PAL_LIB debug ${EP_PREFIX}/lib/libopen-pald.lib optimized ${EP_PREFIX}/lib/libopen-pal.lib)
    SET(OPENMPI_RTE_LIB debug ${EP_PREFIX}/lib/libopen-rted.lib optimized ${EP_PREFIX}/lib/libopen-rte.lib)
    SET(OPENMPI_C_LIB debug ${EP_PREFIX}/lib/libmpid.lib optimized ${EP_PREFIX}/lib/libmpi.lib)
    SET(OPENMPI_CXX_LIB debug ${EP_PREFIX}/lib/libmpi_cxxd.lib optimized ${EP_PREFIX}/lib/libmpi_cxx.lib)
    SET(OPENMPI_EXTRA_LIB advapi32 ws2_32 shlwapi)
    SET(MPI_C_LIBRARIES ${EP_PREFIX}/lib/libopen-pal.lib ${EP_PREFIX}/lib/libopen-rte.lib ${EP_PREFIX}/lib/libmpi.lib advapi32 ws2_32 shlwapi)
    SET(MPI_CXX_LIBRARIES ${EP_PREFIX}/lib/libmpi_cxx.lib)
    SET(CMAKE_PREFIX_PATH ${EP_PREFIX})
    IF(PHYSBAM_USE_OPENMPI)
        LIST(APPEND PHYSBAM_DEPENDENCIES openmpi)
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${OPENMPI_PAL_LIB})
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${OPENMPI_RTE_LIB})
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${OPENMPI_C_LIB})
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${OPENMPI_CXX_LIB})
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${OPENMPI_EXTRA_LIB})
    ENDIF(PHYSBAM_USE_OPENMPI)
ENDIF (WIN32)
IF(PHYSBAM_USE_OPENMPI)
    FIND_PACKAGE(MPI)
    IF(NOT WIN32)
        INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
    ENDIF(NOT WIN32)
ENDIF(PHYSBAM_USE_OPENMPI)

#LAM-MPI - no Windows support
IF(PHYSBAM_USE_LAMMPI)
    FIND_PACKAGE(MPI)
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    LIST(APPEND PHYSBAM_LINK_LIBRARIES ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
ENDIF(PHYSBAM_USE_LAMMPI)

#CUDA
IF(PHYSBAM_USE_CUDA)
    FIND_PACKAGE(CUDA REQUIRED)
    INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
ENDIF(PHYSBAM_USE_CUDA)

#OptiX
IF(PHYSBAM_USE_OPTIX)
    FIND_PACKAGE(OptiX REQUIRED)
    INCLUDE_DIRECTORIES(${OPTIX_INCLUDE})
ENDIF(PHYSBAM_USE_OPTIX)

#Threading Building Blocks - TBB
IF(PHYSBAM_USE_TBB)
    FIND_PACKAGE(TBB REQUIRED)
    INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
ENDIF(PHYSBAM_USE_TBB)

#find external header files
INCLUDE_DIRECTORIES(${EP_PREFIX}/include)

# make it so normal headers can be found
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

#create debug libraries with -debug suffix, release/optimized libraries have no suffix
SET(CMAKE_DEBUG_POSTFIX "-debug")

IF(WIN32)
    INCLUDE(GenerateExportHeader)
ENDIF(WIN32)

#Create shared or static libs based on the variable
IF(PHYSBAM_BUILD_SHARED_LIBS)
    SET(PHYSBAM_LIB_TYPE SHARED)
ELSE(PHYSBAM_BUILD_SHARED_LIBS)
    SET(PHYSBAM_LIB_TYPE STATIC)
ENDIF(PHYSBAM_BUILD_SHARED_LIBS)

#PhysBAM_Tools
FILE(GLOB_RECURSE physbam_tools_cpp "PhysBAM_Tools/*.cpp")
Create_Source_Group("PhysBAM_Tools" ${PROJECT_SOURCE_DIR} ${physbam_tools_cpp})
FILE(GLOB_RECURSE physbam_tools_h "PhysBAM_Tools/*.h")
Create_Source_Group("PhysBAM_Tools" ${PROJECT_SOURCE_DIR} ${physbam_tools_h})
ADD_LIBRARY(PhysBAM_Tools ${PHYSBAM_LIB_TYPE} ${physbam_tools_cpp} ${physbam_tools_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Tools BASE_NAME PhysBAM_Tools EXPORT_MACRO_NAME PhysBAM_Tools_EXPORT 
        EXPORT_FILE_NAME PhysBAM_Tools_Export.h STATIC_DEFINE PhysBAM_Tools_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ADD_DEPENDENCIES(PhysBAM_Tools ${PHYSBAM_DEPENDENCIES})
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
TARGET_LINK_LIBRARIES(PhysBAM_Tools ${PHYSBAM_LINK_LIBRARIES})
SET(PHYSBAM_INSTALL_TARGETS PhysBAM_Tools)

#PhysBAM_Geometry
FILE(GLOB_RECURSE physbam_geometry_cpp "PhysBAM_Geometry/*.cpp")
Create_Source_Group("PhysBAM_Geometry" ${PROJECT_SOURCE_DIR} ${physbam_geometry_cpp})
FILE(GLOB_RECURSE physbam_geometry_h "PhysBAM_Geometry/*.h")
Create_Source_Group("PhysBAM_Geometry" ${PROJECT_SOURCE_DIR} ${physbam_geometry_h})
ADD_LIBRARY(PhysBAM_Geometry ${PHYSBAM_LIB_TYPE} ${physbam_geometry_cpp} ${physbam_geometry_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Geometry BASE_NAME PhysBAM_Geometry EXPORT_MACRO_NAME PhysBAM_Geometry_EXPORT 
        EXPORT_FILE_NAME PhysBAM_Geometry_Export.h STATIC_DEFINE PhysBAM_Geometry_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Geometry PhysBAM_Tools)
TARGET_LINK_LIBRARIES(PhysBAM_Geometry PhysBAM_Tools)
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Geometry)

#PhysBAM_Solids
FILE(GLOB_RECURSE physbam_solids_cpp "PhysBAM_Solids/*.cpp")
Create_Source_Group("PhysBAM_Solids" ${PROJECT_SOURCE_DIR} ${physbam_solids_cpp})
FILE(GLOB_RECURSE physbam_solids_h "PhysBAM_Solids/*.h")
Create_Source_Group("PhysBAM_Solids" ${PROJECT_SOURCE_DIR} ${physbam_solids_h})
ADD_LIBRARY(PhysBAM_Solids ${PHYSBAM_LIB_TYPE} ${physbam_solids_cpp} ${physbam_solids_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Solids BASE_NAME PhysBAM_Solids EXPORT_MACRO_NAME PhysBAM_Solids_EXPORT 
        EXPORT_FILE_NAME PhysBAM_Solids_Export.h STATIC_DEFINE PhysBAM_Solids_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Solids PhysBAM_Geometry)
TARGET_LINK_LIBRARIES(PhysBAM_Solids PhysBAM_Geometry)
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Solids)

#PhysBAM_Fluids
FILE(GLOB_RECURSE physbam_fluids_cpp "PhysBAM_Fluids/*.cpp")
Create_Source_Group("PhysBAM_Fluids" ${PROJECT_SOURCE_DIR} ${physbam_fluids_cpp})
FILE(GLOB_RECURSE physbam_fluids_h "PhysBAM_Fluids/*.h")
Create_Source_Group("PhysBAM_Fluids" ${PROJECT_SOURCE_DIR} ${physbam_fluids_h})
ADD_LIBRARY(PhysBAM_Fluids ${PHYSBAM_LIB_TYPE} ${physbam_fluids_cpp} ${physbam_fluids_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Fluids BASE_NAME PhysBAM_Fluids EXPORT_MACRO_NAME PhysBAM_Fluids_EXPORT 
        EXPORT_FILE_NAME PhysBAM_Fluids_Export.h STATIC_DEFINE PhysBAM_Fluids_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Fluids PhysBAM_Geometry)
TARGET_LINK_LIBRARIES(PhysBAM_Fluids PhysBAM_Dynamics PhysBAM_Geometry)
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Fluids)

#PhysBAM_Dynamics
FILE(GLOB_RECURSE physbam_dynamics_cpp "PhysBAM_Dynamics/*.cpp")
Create_Source_Group("PhysBAM_Dynamics" ${PROJECT_SOURCE_DIR} ${physbam_dynamics_cpp})
FILE(GLOB_RECURSE physbam_dynamics_h "PhysBAM_Dynamics/*.h")
Create_Source_Group("PhysBAM_Dynamics" ${PROJECT_SOURCE_DIR} ${physbam_dynamics_h})
ADD_LIBRARY(PhysBAM_Dynamics ${PHYSBAM_LIB_TYPE} ${physbam_dynamics_cpp} ${physbam_dynamics_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Dynamics BASE_NAME PhysBAM_Dynamics EXPORT_MACRO_NAME PhysBAM_Dynamics_EXPORT 
        EXPORT_FILE_NAME PhysBAM_Dynamics_Export.h STATIC_DEFINE PhysBAM_Dynamics_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Dynamics PhysBAM_Solids PhysBAM_Fluids)
TARGET_LINK_LIBRARIES(PhysBAM_Dynamics PhysBAM_Fluids PhysBAM_Solids)
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Dynamics)

#PhysBAM_Rendering/PhysBAM_OpenGL
FILE(GLOB_RECURSE physbam_rendering_opengl_cpp "PhysBAM_Rendering/PhysBAM_OpenGL*/*.cpp")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_opengl_cpp})
FILE(GLOB_RECURSE physbam_rendering_opengl_h "PhysBAM_Rendering/PhysBAM_OpenGL*/*.h")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_opengl_h})
ADD_LIBRARY(PhysBAM_Rendering-OpenGL ${PHYSBAM_LIB_TYPE} ${physbam_rendering_opengl_cpp}
            ${physbam_rendering_opengl_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Rendering-OpenGL BASE_NAME PhysBAM_Rendering-OpenGL 
        EXPORT_MACRO_NAME PhysBAM_Rendering-OpenGL_EXPORT EXPORT_FILE_NAME PhysBAM_Rendering-OpenGL_Export.h
        STATIC_DEFINE PhysBAM_Rendering-OpenGL_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Rendering-OpenGL PhysBAM_Dynamics)
IF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
    ADD_DEPENDENCIES(PhysBAM_Rendering-OpenGL freeglut)
ENDIF(PHYSBAM_BUILD_EXTERNAL_LIBRARIES)
TARGET_LINK_LIBRARIES(PhysBAM_Rendering-OpenGL PhysBAM_Dynamics ${OPENGL_LIBRARY} ${GLUT_LIBRARY})
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Rendering-OpenGL)

#PhysBAM_Rendering/PhysBAM_Ray_Tracing
FILE(GLOB_RECURSE physbam_rendering_ray_tracing_cpp "PhysBAM_Rendering/PhysBAM_Ray_Tracing/*.cpp")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_ray_tracing_cpp})
FILE(GLOB_RECURSE physbam_rendering_ray_tracing_h "PhysBAM_Rendering/PhysBAM_Ray_Tracing/*.h")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_ray_tracing_h})
ADD_LIBRARY(PhysBAM_Rendering-Ray_Tracing ${PHYSBAM_LIB_TYPE} ${physbam_rendering_ray_tracing_cpp} 
            ${physbam_rendering_ray_tracing_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Rendering-Ray_Tracing BASE_NAME PhysBAM_Rendering-Ray_Tracing 
        EXPORT_MACRO_NAME PhysBAM_Rendering-Ray_Tracing_EXPORT EXPORT_FILE_NAME PhysBAM_Rendering-Ray_Tracing_Export.h 
        STATIC_DEFINE PhysBAM_Rendering-Ray_Tracing_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Rendering-Ray_Tracing PhysBAM_Dynamics)
TARGET_LINK_LIBRARIES(PhysBAM_Rendering-Ray_Tracing PhysBAM_Dynamics)
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Rendering-Ray_Tracing)

#PhysBAM_Rendering/PhysBAM_OptiX
IF(PHYSBAM_USE_OPTIX)
FILE(GLOB_RECURSE physbam_rendering_optix_cpp "PhysBAM_Rendering/PhysBAM_OptiX/*.cpp")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_optix_cpp})
FILE(GLOB_RECURSE physbam_rendering_optix_h "PhysBAM_Rendering/PhysBAM_OptiX/*.h")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_optix_h})
FILE(GLOB_RECURSE physbam_rendering_optix_cu "PhysBAM_Rendering/PhysBAM_OptiX/*.cu")
Create_Source_Group("PhysBAM_Rendering" ${PROJECT_SOURCE_DIR} ${physbam_rendering_optix_cu})
SET(CUDA_GENERATED_OUTPUT_DIR_OLD ${CUDA_GENERATED_OUTPUT_DIR})
SET(CUDA_GENERATED_OUTPUT_DIR ${PROJECT_SOURCE_DIR}/PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Shaders/Compiled)
#CUDA 4/5 does not compile with MSVC11 yet, so you need to use MSVC10
#this should go away if/when the CUDA toolkit is updated to support MSVC11
#NOTE: You need to remove the block that defines ccbin_flags in FindCUDA.cmake
IF(MSVC11)
    GET_FILENAME_COMPONENT(PHYSBAM_MSVC10_PATH [HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\10.0\\Setup\\VC;ProductDir] ABSOLUTE)
    IF("${PHYSBAM_MSVC10_PATH}" STREQUAL "")
        MESSAGE(FATAL_ERROR "Could not find Visual Studio 2010 on this computer. You need to install both 2010 and 2012 for CUDA to work.")
    ENDIF("${PHYSBAM_MSVC10_PATH}" STREQUAL "")
    SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --compiler-bindir \"${PHYSBAM_MSVC10_PATH}/bin\"")
ENDIF(MSVC11)
CUDA_COMPILE_PTX(physbam_rendering_optix_ptx ${physbam_rendering_optix_cu})
SET(CUDA_GENERATED_OUTPUT_DIR ${CUDA_GENERATED_OUTPUT_DIR_OLD})
ADD_LIBRARY(PhysBAM_Rendering-OptiX ${PHYSBAM_LIB_TYPE} ${physbam_rendering_optix_cpp} ${physbam_rendering_optix_cu} 
            ${physbam_rendering_optix_ptx} ${physbam_rendering_optix_h})
IF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
    GENERATE_EXPORT_HEADER(PhysBAM_Rendering-OptiX BASE_NAME PhysBAM_Rendering-OptiX 
        EXPORT_MACRO_NAME PhysBAM_Rendering-OptiX_EXPORT EXPORT_FILE_NAME PhysBAM_Rendering-OptiX_Export.h 
        STATIC_DEFINE PhysBAM_Rendering-OptiX_Static)
ENDIF(WIN32 AND ${PHYSBAM_LIB_TYPE} STREQUAL SHARED)
ADD_DEPENDENCIES(PhysBAM_Rendering-OptiX PhysBAM_Rendering-OpenGL)
TARGET_LINK_LIBRARIES(PhysBAM_Rendering-OptiX PhysBAM_Rendering-OpenGL ${OPTIX_LIBRARY} ${OPTIXU_LIBRARY} ${CUDA_CUDART_LIBRARY})
LIST(APPEND PHYSBAM_INSTALL_TARGETS PhysBAM_Rendering-OptiX)
ENDIF(PHYSBAM_USE_OPTIX)

#Save a targets file that can be imported into a project for easy linking against Public_Library
EXPORT(TARGETS ${PHYSBAM_INSTALL_TARGETS} FILE Public_Library-targets.cmake)
