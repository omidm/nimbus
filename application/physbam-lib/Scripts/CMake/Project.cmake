###########################################################################
#Build Type
###########################################################################
IF(NOT WIN32)
    #Set a build type based on the folder name if nothing is defined
    IF(NOT CMAKE_BUILD_TYPE)
        GET_FILENAME_COMPONENT(PHYSBAM_BUILD_FOLDER_NAME ${CMAKE_CURRENT_BINARY_DIR} NAME)
        SET(CMAKE_BUILD_TYPE ${PHYSBAM_BUILD_FOLDER_NAME} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
    ENDIF(NOT CMAKE_BUILD_TYPE)
ENDIF(NOT WIN32)

###########################################################################
#Public Library
###########################################################################

#tell CMake where PhysBAM root is
SET(PHYSBAM_ROOT ${PROJECT_SOURCE_DIR}/../..)
#include the build options file
INCLUDE(${PHYSBAM_ROOT}/Scripts/CMake/Build_Options.cmake)
#include the compiler flags file
INCLUDE(${PHYSBAM_ROOT}/Scripts/CMake/Compiler_Flags.cmake)
#include the source group file
INCLUDE(${PHYSBAM_ROOT}/Scripts/CMake/Source_Group.cmake)
#this include is static
INCLUDE(${PHYSBAM_ROOT}/Scripts/CMake/Public_Library.cmake)
#this include is created upon install of Public_Library - build it separately first if it doesn't exist
IF(WIN32)
    INCLUDE(${PHYSBAM_ROOT}/Public_Library/Build/Public_Library-targets.cmake)
ELSE(WIN32) #UNIX or Mac OS
    INCLUDE(${PHYSBAM_ROOT}/Public_Library/Build/${CMAKE_BUILD_TYPE}/Public_Library-targets.cmake)
ENDIF(WIN32)
#find the correct libraries (does not append to project executables)
SET(CMAKE_DEBUG_POSTFIX "-debug")
#load cached variables from Public_Library to prevent erroneous builds by variable/compiler mismatch
IF(WIN32)
    SET(PHYSBAM_CMAKE_CACHE_PATH ${PHYSBAM_ROOT}/Public_Library/Build)
ELSE(WIN32)
    SET(PHYSBAM_CMAKE_CACHE_PATH ${PHYSBAM_ROOT}/Public_Library/Build/${CMAKE_BUILD_TYPE})
ENDIF(WIN32)
LOAD_CACHE(${PHYSBAM_CMAKE_CACHE_PATH})

###########################################################################
#External Libraries for Projects
###########################################################################

#Make sure the linker can find external libraries
IF(PHYSBAM_USE_LIBPNG)
    LIST(APPEND LIBRARIES_TO_LINK ${PNG_LIBRARY})
ENDIF()

IF(PHYSBAM_USE_LIBJPEG)
    LIST(APPEND LIBRARIES_TO_LINK ${JPEG_LIBRARY})
ENDIF()

IF(PHYSBAM_USE_OPENMPI OR PHYSBAM_USE_LAMMPI)
    LIST(APPEND LIBRARIES_TO_LINK ${MPI_LIBRARY})
ENDIF(PHYSBAM_USE_OPENMPI OR PHYSBAM_USE_LAMMPI)

#Handle project-only libraries
#Compile project-only external libraries? Yes = ON, No = OFF (default)
SET(PHYSBAM_BUILD_PROJECT_SPECIFIC_EXTERNAL_LIBRARIES OFF CACHE BOOL "Build external libraries that only projects can use")

#Use CLAPACK? Yes = ON, No = OFF (default)
SET(PHYSBAM_PROJECT_USE_CLAPACK OFF CACHE BOOL "Use CLAPACK in this project")

IF(PHYSBAM_BUILD_PROJECT_SPECIFIC_EXTERNAL_LIBRARIES)
    ExternalProject_Add(clapack 
        PREFIX ${EL_PREFIX}
        CMAKE_ARGS -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} 
            -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER} -D CMAKE_INSTALL_PREFIX=${EL_PREFIX})
ENDIF(PHYSBAM_BUILD_PROJECT_SPECIFIC_EXTERNAL_LIBRARIES)

IF (WIN32)
    SET(F2C_LIBRARY debug ${EL_PREFIX}/lib/libf2cd.lib optimized ${EL_PREFIX}/lib/libf2c.lib) 
    SET(BLAS_LIBRARY debug ${EL_PREFIX}/lib/blasd.lib optimized ${EL_PREFIX}/lib/blas.lib)
    SET(CLAPACK_LIBRARY debug ${EL_PREFIX}/lib/lapackd.lib optimized ${EL_PREFIX}/lib/lapack.lib)
ELSE (WIN32)
    SET(F2C_LIBRARY debug ${EL_PREFIX}/lib/libf2cd.a optimized ${EL_PREFIX}/lib/libf2c.a)
    SET(BLAS_LIBRARY debug ${EL_PREFIX}/lib/libblasd.a optimized ${EL_PREFIX}/lib/libblas.a)
    SET(CLAPACK_LIBRARY debug ${EL_PREFIX}/lib/liblapackd.a optimized ${EL_PREFIX}/lib/liblapack.a)
ENDIF (WIN32)

IF(PHYSBAM_PROJECT_USE_CLAPACK)
    ADD_DEFINITIONS(-DUSE_CLAPACK)
    LIST(APPEND LIBRARIES_TO_LINK ${CLAPACK_LIBRARY})
    LIST(APPEND LIBRARIES_TO_LINK ${F2C_LIBRARY})
    LIST(APPEND LIBRARIES_TO_LINK ${BLAS_LIBRARY})
ENDIF(PHYSBAM_PROJECT_USE_CLAPACK)

#Use QT4? Yes = ON, No = OFF (default)
SET(PHYSBAM_PROJECT_USE_QT4 OFF CACHE BOOL "Use QT4 in this project (install it first)")

IF(PHYSBAM_PROJECT_USE_QT4)
    FIND_PACKAGE(Qt4 REQUIRED)
    SET(PHYSBAM_PROJECT_USE_QT4NETWORK OFF CACHE BOOL "Enable QT4 Network Module")
    SET(QT_USE_QTNETWORK ${PHYSBAM_PROJECT_USE_QT4NETWORK})
    SET(PHYSBAM_PROJECT_USE_QT4OPENGL OFF CACHE BOOL "Enable QT4 OpenGL Module")
    SET(QT_USE_QTOPENGL ${PHYSBAM_PROJECT_USE_QT4OPENGL})
    SET(PHYSBAM_PROJECT_USE_QT4SQL OFF CACHE BOOL "Enable QT4 SQL Module")
    SET(QT_USE_QTSQL ${PHYSBAM_PROJECT_USE_QT4SQL})
    SET(PHYSBAM_PROJECT_USE_QT4XML OFF CACHE BOOL "Enable QT4 XML Module")
    SET(QT_USE_QTXML ${PHYSBAM_PROJECT_USE_QT4XML})
    SET(PHYSBAM_PROJECT_USE_QT4SVG OFF CACHE BOOL "Enable QT4 SVG Module")
    SET(QT_USE_QTSVG ${PHYSBAM_PROJECT_USE_QT4SVG})
    SET(PHYSBAM_PROJECT_USE_QT4TEST OFF CACHE BOOL "Enable QT4 Test Module")
    SET(QT_USE_QTTEST ${PHYSBAM_PROJECT_USE_QT4TEST})
    SET(PHYSBAM_PROJECT_USE_QT4DBUS OFF CACHE BOOL "Enable QT4 DBUS Module")
    SET(QT_USE_QTDBUS ${PHYSBAM_PROJECT_USE_QT4DBUS})
    SET(PHYSBAM_PROJECT_USE_QT4SCRIPT OFF CACHE BOOL "Enable QT4 Script Module")
    SET(QT_USE_QTSCRIPT ${PHYSBAM_PROJECT_USE_QT4SCRIPT})
    SET(PHYSBAM_PROJECT_USE_QT4WEBKIT OFF CACHE BOOL "Enable QT4 Webkit Module")
    SET(QT_USE_QTWEBKIT ${PHYSBAM_PROJECT_USE_QT4WEBKIT})
    SET(PHYSBAM_PROJECT_USE_QT4XMLPATTERNS OFF CACHE BOOL "Enable QT4 XML Patterns Module")
    SET(QT_USE_QTXMLPATTERNS ${PHYSBAM_PROJECT_USE_QT4XMLPATTERNS})
    SET(PHYSBAM_PROJECT_USE_QT4PHONON OFF CACHE BOOL "Enable QT4 Phonon Module")
    SET(QT_USE_PHONON ${PHYSBAM_PROJECT_USE_QT4PHONON})
    INCLUDE(${QT_USE_FILE})
    QT4_WRAP_CPP(${PROJECT_NAME}_qt_h_proc ${${PROJECT_NAME}_qt_h})
    QT4_WRAP_UI(${PROJECT_NAME}_qt_forms_proc ${${PROJECT_NAME}_qt_forms})
    QT4_ADD_RESOURCES(${PROJECT_NAME}_qt_resources_proc ${${PROJECT_NAME}_qt_resources})
    LIST(APPEND LIBRARIES_TO_LINK ${QT_LIBRARIES})
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
ENDIF(PHYSBAM_PROJECT_USE_QT4)

###########################################################################
#Project Creation
###########################################################################

#add folder hierarchy for IDEs
SET(PROJECT_SOURCE_FILES ${${PROJECT_NAME}_h} ${${PROJECT_NAME}_cpp})
Create_Source_Group("" ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_FILES})
IF(PHYSBAM_PROJECT_USE_QT4)
    IF(NOT ${${PROJECT_NAME}_qt_h} STREQUAL "")
        Create_Source_Group("" ${PROJECT_SOURCE_DIR} ${${PROJECT_NAME}_qt_h})
        SOURCE_GROUP("\\QT moc Generated Files" FILES ${${PROJECT_NAME}_qt_h_proc})
    ENDIF(NOT ${${PROJECT_NAME}_qt_h} STREQUAL "")
    LIST(APPEND PROJECT_SOURCE_FILES ${${PROJECT_NAME}_qt_h} ${${PROJECT_NAME}_qt_h_proc}
         ${${PROJECT_NAME}_qt_forms_proc} ${${PROJECT_NAME}_qt_resources_proc})
ENDIF(PHYSBAM_PROJECT_USE_QT4)

ADD_EXECUTABLE(${PROJECT_NAME} ${PROJECT_SOURCE_FILES})
ADD_DEPENDENCIES(${PROJECT_NAME} PhysBAM)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${LIBRARIES_TO_LINK})

IF(PHYSBAM_BUILD_PROJECT_SPECIFIC_EXTERNAL_LIBRARIES)
    IF(PHYSBAM_PROJECT_USE_CLAPACK)
        ADD_DEPENDENCIES(${PROJECT_NAME} clapack)
    ENDIF(PHYSBAM_PROJECT_USE_CLAPACK)
ENDIF(PHYSBAM_BUILD_PROJECT_SPECIFIC_EXTERNAL_LIBRARIES)

SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON 
        DEBUG_OUTPUT_NAME ${PROJECT_NAME}-debug RELEASE_OUTPUT_NAME ${PROJECT_NAME})

IF(NOT SKIP_INSTALL_LIBRARIES AND NOT SKIP_INSTALL_ALL)
    INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR})
ENDIF()

IF(MSVC)
    ADD_CUSTOM_TARGET(AUTOMATIC_INSTALL ALL
        COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_CONFIG_NAME=${CMAKE_CFG_INTDIR}
        -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
    ADD_DEPENDENCIES(AUTOMATIC_INSTALL ${PROJECT_NAME})
ELSE(MSVC)
    ADD_CUSTOM_TARGET(AUTOMATIC_INSTALL ALL
        COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_CONFIG_NAME=${CMAKE_BUILD_TYPE}
        -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
    ADD_DEPENDENCIES(AUTOMATIC_INSTALL ${PROJECT_NAME})
ENDIF(MSVC)
