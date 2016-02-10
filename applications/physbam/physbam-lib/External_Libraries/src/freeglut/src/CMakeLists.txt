SET(FREEGLUT_INT_HEADERS freeglut_internal.h
                         freeglut_teapot_data.h)
                         
SET(FREEGLUT_EXT_HEADERS "../include/GL/freeglut.h"
                         "../include/GL/freeglut_ext.h"
                         "../include/GL/freeglut_std.h"
                         "../include/GL/glut.h")
                         
SET(FREEGLUT_SRC freeglut_callbacks
                 freeglut_cursor
                 freeglut_display
                 freeglut_ext
                 freeglut_font
                 freeglut_font_data
                 freeglut_gamemode
                 freeglut_geometry
                 freeglut_glutfont_definitions
                 freeglut_init
                 freeglut_input_devices
                 freeglut_joystick
                 freeglut_main
                 freeglut_menu
                 freeglut_misc
                 freeglut_overlay
                 freeglut_spaceball
                 freeglut_state
                 freeglut_stroke_mono_roman
                 freeglut_stroke_roman
                 freeglut_structure
                 freeglut_teapot
                 freeglut_videoresize
                 freeglut_window)

SET(FREEGLUT_SHARED_SRC ${FREEGLUT_SRC})
SET(FREEGLUT_STATIC_SRC ${FREEGLUT_SRC})

IF(WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    LIST(APPEND FREEGLUT_SHARED_SRC freeglutdll.def freeglut.rc)
ENDIF(WIN32)

ADD_DEFINITIONS(-DHAVE_CONFIG_H)

IF(${FREEGLUT_STATIC})
	ADD_LIBRARY(freeglut STATIC ${FREEGLUT_STATIC_SRC} ${FREEGLUT_INT_HEADERS} ${FREEGLUT_EXT_HEADERS})
	ADD_DEFINITIONS(-DFREEGLUT_STATIC)
	SET_TARGET_PROPERTIES(freeglut PROPERTIES OUTPUT_NAME "freeglut_static")
ELSE()
	ADD_LIBRARY(freeglut SHARED ${FREEGLUT_SHARED_SRC} ${FREEGLUT_INT_HEADERS} ${FREEGLUT_EXT_HEADERS})
	ADD_DEFINITIONS(-DFREEGLUT_EXPORTS)
ENDIF()

#PhysBAM Edit: Need to install files to External_Libraries
IF(NOT SKIP_INSTALL_LIBRARIES AND NOT SKIP_INSTALL_ALL)
    INSTALL(TARGETS freeglut
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib )
ENDIF()
IF(NOT SKIP_INSTALL_HEADERS AND NOT SKIP_INSTALL_ALL )
    INSTALL(FILES ${FREEGLUT_EXT_HEADERS} DESTINATION "include/GL")
    IF(WIN32)
        INSTALL(FILES "../include/GL/glext.h" DESTINATION "include/GL")
    ENDIF()
ENDIF()