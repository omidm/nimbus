
#
#  Copyright (c) 2010 NVIDIA Corporation.  All rights reserved.
#
#  NVIDIA Corporation and its licensors retain all intellectual property and proprietary
#  rights in and to this software, related documentation and any modifications thereto.
#  Any use, reproduction, disclosure or distribution of this software and related
#  documentation without an express license agreement from NVIDIA Corporation is strictly
#  prohibited.
#
#  TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, THIS SOFTWARE IS PROVIDED *AS IS*
#  AND NVIDIA AND ITS SUPPLIERS DISCLAIM ALL WARRANTIES, EITHER EXPRESS OR IMPLIED,
#  INCLUDING, BUT NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  PARTICULAR PURPOSE.  IN NO EVENT SHALL NVIDIA OR ITS SUPPLIERS BE LIABLE FOR ANY
#  SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES WHATSOEVER (INCLUDING, WITHOUT
#  LIMITATION, DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION, LOSS OF
#  BUSINESS INFORMATION, OR ANY OTHER PECUNIARY LOSS) ARISING OUT OF THE USE OF OR
#  INABILITY TO USE THIS SOFTWARE, EVEN IF NVIDIA HAS BEEN ADVISED OF THE POSSIBILITY OF
#  SUCH DAMAGES
#

# Locate the OptiX distribution.  Search relative to the SDK first, then look in the system.

# Our initial guess will be within the SDK.
# PhysBAM Edit: Likely not in the SDK directory on Windows. Set the environment variable OPTIX_SDK_DIR.
set(OPTIX_INSTALL_DIR $ENV{OPTIX_SDK_DIR} CACHE PATH "Path to OptiX installed location.")

# The distribution contains both 32 and 64 bit libraries.  Adjust the library
# search path based on the bit-ness of the build.  (i.e. 64: bin64, lib64; 32:
# bin, lib).  Note that on Mac, the OptiX library is a universal binary, so we
# only need to look in lib and not lib64 for 64 bit builds.
if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND NOT APPLE)
  set(bit_dest "64")
else()
  set(bit_dest "")
endif()

macro(OPTIX_find_api_library name version)
  SET(CAP_NAME name)
  STRING(TOUPPER ${name} CAP_NAME)
  find_library(${CAP_NAME}_LIBRARY
    NAMES ${name}.${version} ${name}
    PATHS "${OPTIX_INSTALL_DIR}/lib${bit_dest}"
    NO_DEFAULT_PATH
    )
  find_library(${CAP_NAME}_LIBRARY
    NAMES ${name}.${version} ${name}
    )
  if(WIN32)
    find_file(${CAP_NAME}_DLL
      NAMES ${name}.${version}.dll
      PATHS "${OPTIX_INSTALL_DIR}/bin${bit_dest}"
      NO_DEFAULT_PATH
      )
    find_file(${CAP_NAME}_DLL
      NAMES ${name}.${version}.dll
      )
  endif()
endmacro()

OPTIX_find_api_library(optix 1)
OPTIX_find_api_library(optixu 1)

# Include
find_path(OPTIX_INCLUDE
  NAMES optix.h
  PATHS "${OPTIX_INSTALL_DIR}/include"
  NO_DEFAULT_PATH
  )
find_path(OPTIX_INCLUDE
  NAMES optix.h
  )

# Check to make sure we found what we were looking for
function(OptiX_report_error error_message)
  if(OptiX_FIND_REQUIRED)
    message(FATAL_ERROR "${error_message}")
  else(OptiX_FIND_REQUIRED)
    if(NOT OptiX_FIND_QUIETLY)
      message(STATUS "${error_message}")
    endif(NOT OptiX_FIND_QUIETLY)
  endif(OptiX_FIND_REQUIRED)
endfunction()

if(NOT OPTIX_LIBRARY)
  OptiX_report_error("optix library not found.  Please locate before proceeding.")
endif()
if(NOT OPTIX_INCLUDE)
  OptiX_report_error("OptiX headers (optix.h and friends) not found.  Please locate before proceeding.")
endif()

# Macro for setting up dummy targets
function(OptiX_add_imported_library name lib_location dll_lib dependent_libs)
  set(CMAKE_IMPORT_FILE_VERSION 1)

  # Create imported target
  add_library(${name} SHARED IMPORTED)

  # Import target "optix" for configuration "Debug"
  if(WIN32)
    set_target_properties(${name} PROPERTIES
      IMPORTED_IMPLIB "${lib_location}"
      #IMPORTED_LINK_INTERFACE_LIBRARIES "glu32;opengl32"
      IMPORTED_LOCATION "${dll_lib}"
      IMPORTED_LINK_INTERFACE_LIBRARIES "${dependent_libs}"
      )
  elseif(UNIX)
    set_target_properties(${name} PROPERTIES
      #IMPORTED_LINK_INTERFACE_LIBRARIES "glu32;opengl32"
      IMPORTED_LOCATION "${lib_location}"
      # We don't have versioned filenames for now, and it may not even matter.
      #IMPORTED_SONAME "${optix_soname}"
      IMPORTED_LINK_INTERFACE_LIBRARIES "${dependent_libs}"
      )
  else()
    # Unknown system, but at least try and provide the minimum required
    # information.
    set_target_properties(${name} PROPERTIES
      IMPORTED_LOCATION "${lib_location}"
      IMPORTED_LINK_INTERFACE_LIBRARIES "${dependent_libs}"
      )
  endif()

  # Commands beyond this point should not need to know the version.
  set(CMAKE_IMPORT_FILE_VERSION)
endfunction()

# Sets up a dummy target
OptiX_add_imported_library(optix "${OPTIX_LIBRARY}" "${OPTIX_DLL}" "${CUDA_CUDA_LIBRARY};${OPENGL_LIBRARIES}")
OptiX_add_imported_library(optixu   "${OPTIXU_LIBRARY}"   "${OPTIXU_DLL}"   "")

# Since liboptix.1.dylib is built with an install name of @rpath, we need to
# compile our samples with the rpath set to where optix exists.
if(APPLE)
  get_filename_component(_optix_path_to_optix "${OPTIX_LIBRARY}" PATH)
  if(_optix_path_to_optix)
    set( _optix_rpath "-Wl,-rpath,${_optix_path_to_optix}" )
  endif()
  get_filename_component(_optix_path_to_optixu "${OPTIXU_LIBRARY}" PATH)
  if(_optixu_path_to_optix)
    if(NOT _optixu_path_to_optix STREQUAL _optix_path_to_optixu)
      # optixu and optix are in different paths.  Make sure there isn't an optixu next to
      # the optix library.
      get_filename_component(_optix_name_of_optixu "${OPTIXU_LIBRARY}" NAME)
      if(EXISTS "${_optix_path_to_optix}/${_optix_name_of_optixu}")
        message(WARNING " optixu library found next to optix library that is not being used.  Due to the way we are usin
g rpath, the copy of optixu next to optix will be used during loading instead of the one you intended.  Consider putting the libraries in the same directory or moving ${OPTIXU_LIBRARY} out of the way.")
      endif()       
    endif()         
    set( _optixu_rpath "-Wl,-rpath,${_optixu_path_to_optix}" )
  endif()           
  set( optix_rpath ${_optixu_rpath} ${_optix_rpath} )
endif()             

