#This file is for including Public_Library in Projects

#Build Public_Library
INCLUDE(ExternalProject)
SET(PB_PREFIX ${PHYSBAM_ROOT}/Public_Library)
#add custom find modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PHYSBAM_ROOT}/Scripts/CMake/Modules)
IF(WIN32)
    SET(PB_PREFIX_BUILD ${PHYSBAM_ROOT}/Public_Library/Build)
ELSE(WIN32)
    SET(PB_PREFIX_BUILD ${PHYSBAM_ROOT}/Public_Library/Build/${CMAKE_BUILD_TYPE})
ENDIF(WIN32)
#make sure the CMAKE_ARGS line never goes above 8000 characters, or the build will fail in VS
ExternalProject_Add(PhysBAM 
    TMP_DIR ${PB_PREFIX_BUILD}/tmp
    STAMP_DIR ${PB_PREFIX_BUILD}/stamp
    SOURCE_DIR ${PB_PREFIX}
    BINARY_DIR ${PB_PREFIX_BUILD}
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DPHYSBAM_USE_LIBPNG=${PHYSBAM_USE_LIBPNG} -DPHYSBAM_USE_LIBJPEG=${PHYSBAM_USE_LIBJPEG} -DPHYSBAM_USE_PTHREADS=${PHYSBAM_USE_PTHREADS}
        -DPHYSBAM_USE_OPENMPI=${PHYSBAM_USE_OPENMPI} -DPHYSBAM_USE_CUDA=${PHYSBAM_USE_CUDA} -DPHYSBAM_COMPILE_WITHOUT_DYADIC_SUPPORT=${PHYSBAM_COMPILE_WITHOUT_DYADIC_SUPPORT}
        -DPHYSBAM_COMPILE_WITHOUT_RLE_SUPPORT=${PHYSBAM_COMPILE_WITHOUT_RLE_SUPPORT} -DPHYSBAM_COMPILE_WITHOUT_DOUBLE_SUPPORT=${PHYSBAM_COMPILE_WITHOUT_DOUBLE_SUPPORT}
        -DPHYSBAM_BUILD_SHARED_LIBS=${PHYSBAM_BUILD_SHARED_LIBS} -DPHYSBAM_USE_LAMMPI=${PHYSBAM_USE_LAMMPI} -DPHYSBAM_USE_OPTIX=${PHYSBAM_USE_OPTIX}
        -DPHYSBAM_COMPILE_WITHOUT_READ_WRITE_SUPPORT=${PHYSBAM_COMPILE_WITHOUT_READ_WRITE_SUPPORT} -DPHYSBAM_USE_TBB=${PHYSBAM_USE_TBB}
    INSTALL_COMMAND "")
SET(EL_PREFIX ${PHYSBAM_ROOT}/External_Libraries)
INCLUDE_DIRECTORIES(${EL_PREFIX}/include)
INCLUDE_DIRECTORIES(${PB_PREFIX})

#trigger a build for PhysBAM so that library changes are incorporated
ExternalProject_Add_Step(PhysBAM CheckRebuild
    COMMENT "Triggering a PhysBAM rebuild, in case files have changed..."
    DEPENDERS build
    DEPENDEES configure
    ALWAYS 1)

#Set paths to each dependency library
IF (WIN32)
    SET(ZLIB_LIBRARY debug ${EL_PREFIX}/lib/zlibd.lib optimized ${EL_PREFIX}/lib/zlib.lib)
ELSE (WIN32)
    SET(ZLIB_LIBRARY debug ${EL_PREFIX}/lib/libzd.so optimized ${EL_PREFIX}/lib/libz.so)
ENDIF (WIN32)

IF (WIN32)
    SET(PNG_LIBRARY debug ${EL_PREFIX}/lib/libpng15d.lib optimized ${EL_PREFIX}/lib/libpng15.lib)
ELSE (WIN32)
    SET(PNG_LIBRARY debug ${EL_PREFIX}/lib/libpng15d.so optimized ${EL_PREFIX}/lib/libpng15.so)
ENDIF (WIN32)

IF (WIN32)
    SET(JPEG_LIBRARY debug ${EL_PREFIX}/lib/jpegd.lib optimized ${EL_PREFIX}/lib/jpeg.lib)
ELSE (WIN32)
    SET(JPEG_LIBRARY debug ${EL_PREFIX}/lib/libjpegd.a optimized ${EL_PREFIX}/lib/libjpeg.a)
ENDIF (WIN32)

IF (WIN32)
    SET(PTHREADS_LIBRARY debug ${EL_PREFIX}/lib/pthreadVC2d.lib optimized ${EL_PREFIX}/lib/pthreadVC2.lib)
ELSE (WIN32)
    SET(PTHREADS_LIBRARY ${CMAKE_THREAD_LIBS_INIT})
ENDIF (WIN32)

IF (WIN32)
    SET(FLTK_LIBRARY debug ${EL_PREFIX}/lib/fltkd.lib ${EL_PREFIX}/lib/fltkformsd.lib ${EL_PREFIX}/lib/fltkgld.lib ${EL_PREFIX}/lib/fltkimagesd.lib
        optimized ${EL_PREFIX}/lib/fltk.lib ${EL_PREFIX}/lib/fltkforms.lib ${EL_PREFIX}/lib/fltkgl.lib ${EL_PREFIX}/lib/fltkimages.lib)
    SET(FLTK_INCLUDE_DIR ${EL_PREFIX}/include)
ELSE(WIN32)

ENDIF(WIN32)

IF (WIN32)
    SET(GLUT_LIBRARY debug ${EL_PREFIX}/lib/freeglut.lib optimized ${EL_PREFIX}/lib/freeglut.lib)
ELSE (WIN32)
    SET(GLUT_LIBRARY debug ${EL_PREFIX}/lib/libfreeglut.so optimized ${EL_PREFIX}/lib/libfreeglut.so)
ENDIF (WIN32)

FIND_PACKAGE(OpenGL)

IF (WIN32)
    SET(OPENMPI_PAL_LIB debug ${EL_PREFIX}/lib/libopen-pald.lib optimized ${EL_PREFIX}/lib/libopen-pal.lib)
    SET(OPENMPI_RTE_LIB debug ${EL_PREFIX}/lib/libopen-rted.lib optimized ${EL_PREFIX}/lib/libopen-rte.lib)
    SET(OPENMPI_C_LIB debug ${EL_PREFIX}/lib/libmpid.lib optimized ${EL_PREFIX}/lib/libmpi.lib)
    SET(OPENMPI_CXX_LIB debug ${EL_PREFIX}/lib/libmpi_cxxd.lib optimized ${EL_PREFIX}/lib/libmpi_cxx.lib)
    SET(OPENMPI_EXTRA_LIB advapi32 ws2_32 shlwapi)
    SET(MPI_C_LIBRARIES ${EL_PREFIX}/lib/libopen-pal.lib ${EL_PREFIX}/lib/libopen-rte.lib ${EL_PREFIX}/lib/libmpi.lib advapi32 ws2_32 shlwapi)
    SET(MPI_CXX_LIBRARIES ${EL_PREFIX}/lib/libmpi_cxx.lib)
    SET(CMAKE_PREFIX_PATH ${EL_PREFIX})
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
    IF(NOT WIN32)
        INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
        LIST(APPEND PHYSBAM_LINK_LIBRARIES ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
    ENDIF(NOT WIN32)
ENDIF(PHYSBAM_USE_OPENMPI)

IF(PHYSBAM_USE_LAMMPI)
    INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
    LIST(APPEND PHYSBAM_LINK_LIBRARIES ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
ENDIF(PHYSBAM_USE_LAMMPI)

IF(PHYSBAM_USE_CUDA)
    INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
ENDIF(PHYSBAM_USE_CUDA)

IF(PHYSBAM_USE_OPTIX)
    INCLUDE_DIRECTORIES(${OPTIX_INCLUDE})
ENDIF(PHYSBAM_USE_OPTIX)

IF(PHYSBAM_USE_TBB)
    INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
ENDIF(PHYSBAM_USE_TBB)

SET(CMAKE_DEBUG_POSTFIX "-debug")

#The following adds the sources for Public_Library to a dummy project in VS so
#that IntelliSense works properly within projects
IF(MSVC)
    FILE(GLOB_RECURSE physbam_source_cpp "${PB_PREFIX}/*.cpp")
    Create_Source_Group("" ${PB_PREFIX} ${physbam_source_cpp})
    FILE(GLOB_RECURSE physbam_source_h "${PB_PREFIX}/*.h")
    Create_Source_Group("" ${PB_PREFIX} ${physbam_source_h})
    FILE(GLOB_RECURSE physbam_source_cu "${PB_PREFIX}/*.cu")
    Create_Source_Group("" ${PB_PREFIX} ${physbam_source_cu})
    ADD_LIBRARY(PhysBAM-Source EXCLUDE_FROM_ALL ${physbam_source_cpp} ${physbam_source_h} ${physbam_source_cu})
    SET_TARGET_PROPERTIES(PhysBAM-Source PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD 1)
ENDIF(MSVC)
