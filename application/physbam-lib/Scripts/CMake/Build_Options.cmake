#These options can be overriden by project-specific configurations

#New options needed to be passed to the external project using -D<option_name>=${<option_name>} in 
#Public_Library.cmake for project-specific configurations to be passed to the library,
#and an IF(<option_name>) block needs to be added to Compiler_Flags.cmake

#Compile with libpng support? Yes = ON (default), No = OFF
SET(PHYSBAM_USE_LIBPNG ON CACHE BOOL "Include PNG support")
#Compile with libjpeg support? Yes = ON (default), No = OFF
SET(PHYSBAM_USE_LIBJPEG ON CACHE BOOL "Include JPEG support")
#Compile with threading support? Yes = ON, No = OFF (default)
SET(PHYSBAM_USE_PTHREADS ON CACHE BOOL "Include threading support")
#Compile with OpenMPI support? Yes = ON, No = OFF (default)
SET(PHYSBAM_USE_OPENMPI OFF CACHE BOOL "Include OpenMPI support (not compatible with LAM-MPI)")
#Compile with LAM-MPI support? Yes = ON, No = OFF (default)
SET(PHYSBAM_USE_LAMMPI OFF CACHE BOOL "Include LAM-MPI support (not compatible with OpenMPI or Windows)")
#Compile with CUDA support? Yes = ON, No = OFF (default)
SET(PHYSBAM_USE_CUDA OFF CACHE BOOL "Include NVIDIA CUDA support")
#Compile with OptiX support? Yes = ON, No = OFF (default)
SET(PHYSBAM_USE_OPTIX OFF CACHE BOOL "Include NVIDIA OptiX support (requires CUDA)")
#Compile with TBB support? Yes = ON, No = OFF (default)
SET(PHYSBAM_USE_TBB OFF CACHE BOOL "Include Intel Threading Building Block support")
#Compile without dyadic support? Yes = ON (default), No = OFF
SET(PHYSBAM_COMPILE_WITHOUT_DYADIC_SUPPORT ON CACHE BOOL "Don't include dyadic support")
#Compile without run-length encoding support? Yes = ON (default), No = OFF
SET(PHYSBAM_COMPILE_WITHOUT_RLE_SUPPORT ON CACHE BOOL "Don't include RLE support")
#Compile without double precision support? Yes = ON, No = OFF (default)
SET(PHYSBAM_COMPILE_WITHOUT_DOUBLE_SUPPORT OFF CACHE BOOL "Don't include double support")
#Compile without read/write support? Yes = ON, No = OFF (default)
SET(PHYSBAM_COMPILE_WITHOUT_READ_WRITE_SUPPORT OFF CACHE BOOL "Don't include read/write support")
#Compile shared or static libraries? SHARED = ON, STATIC = OFF (default)
SET(PHYSBAM_BUILD_SHARED_LIBS OFF CACHE BOOL "Compile shared (on) or static (off) libraries")
#Compile external libraries? Yes = ON (default), No = OFF
SET(PHYSBAM_BUILD_EXTERNAL_LIBRARIES ON CACHE BOOL "Build external libraries - turn this off after the very first time")