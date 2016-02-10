#PhysBAM Edit: Custom CMake File
PROJECT(PTHREADS_W32 C)
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

SET(DLL_VER    2)

SET(PTHREADS_W32_LIB_TYPE VC)

SET(CLEANUP __CLEANUP_C)

SET(CMAKE_C_FLAGS_RELEASE "/DHAVE_PTW32_CONFIG_H /W3 /MD /nologo /O2 /Ob2")
SET(CMAKE_C_FLAGS_DEBUG "/DHAVE_PTW32_CONFIG_H /W3 /MDd /nologo /Z7")

set(PTHREADS_LIB_NAME pthread${PTHREADS_W32_LIB_TYPE}${DLL_VER})
set(CMAKE_DEBUG_POSTFIX "d")

SET(PTHREADS_W32_HDRS
    config.h implement.h semaphore.h pthread.h need_errno.h)

SET(PTHREADS_W32_SRCS
    attr.c barrier.c cancel.c cleanup.c condvar.c create.c dll.c
    autostatic.c errno.c exit.c fork.c global.c misc.c mutex.c
    nonportable.c private.c rwlock.c sched.c semaphore.c signal.c
    spin.c sync.c tsd.c)

#The following line is probably unnecessary, but VS doesn't find sched.h any other way
INCLUDE_DIRECTORIES("./")
ADD_LIBRARY(${PTHREADS_LIB_NAME} SHARED ${PTHREADS_W32_HDRS} ${PTHREADS_W32_SRCS})
SET_TARGET_PROPERTIES(${PTHREADS_LIB_NAME} PROPERTIES VERSION 2.9)

SET(CMAKE_INSTALL_PREFIX "../../")

SET(PTHREADS_W32_HDRS_INST
    config.h semaphore.h pthread.h sched.h)

IF(NOT SKIP_INSTALL_LIBRARIES AND NOT SKIP_INSTALL_ALL)
    INSTALL(TARGETS ${PTHREADS_LIB_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib )
ENDIF()
IF(NOT SKIP_INSTALL_HEADERS AND NOT SKIP_INSTALL_ALL )
    INSTALL(FILES ${PTHREADS_W32_HDRS_INST} DESTINATION include)
ENDIF()