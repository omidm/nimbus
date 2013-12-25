
TARGET = face_simple_rw face_overlap_rw particles_simple_rw particles_overlap_rw
all: check $(TARGET)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

PHYSBAM_ROOT = $(NIMBUS_ROOT)/application/physbam-lib/Public_Library/Build/Debug

IFLAGS += -I$(LDIR) -I$(NIMBUS_ROOT)/application/physbam-lib/Public_Library
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

LFLAGS += -lpthread -lnimbus -lboost_thread-mt -lboost_system-mt -lprotobuf
#LFLAGS += -lPhysBAM_Dynamics -lPhysBAM_Fluids -lPhysBAM_Geometry
#LFLAGS += -lPhysBAM_Solids -lPhysBAM_Tools
LFLAGS += -lz

DYNAMICS = $(PHYSBAM_ROOT)/libPhysBAM_Dynamics.a
FLUIDS   = $(PHYSBAM_ROOT)/libPhysBAM_Fluids.a
GEOMETRY = $(PHYSBAM_ROOT)/libPhysBAM_Geometry.a
SOLIDS   = $(PHYSBAM_ROOT)/libPhysBAM_Solids.a
TOOLS    = $(PHYSBAM_ROOT)/libPhysBAM_Tools.a
GL       = $(PHYSBAM_ROOT)/libPhysBAM_Rendering-OpenGL.a


LDFLAGS += -L. -L/usr/local/lib
LDFLAGS += -L$(NIMBUS_ROOT)/application/physbam-lib/Public_Library/Build/Release
LDFLAGS += -L$(NIMBUS_ROOT)/application/physbam-lib/External_Libraries/lib/

# I derived this order for the linkages from looking at what PhysBAM's
# cmake files does; it is fragile and frustrating to debug, so change
# at your own risk. -pal
PHYSBAM_LIBS = $(DYNAMICS) $(FLUIDS) $(DYNAMICS) $(FLUIDS) $(SOLIDS) $(GEOMETRY) $(TOOLS)

PHYSBAM_DLIBS = $(subst .a,-debug.a,$(PHYSBAM_LIBS))

CFLAGS += -pthread

# Add the dependencies here
%.o: %.cc $(NIMBUS_ROOT)/data/physbam/translator_physbam.h
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

face_simple_rw: face_simple_rw.o
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $@.o $(PHYSBAM_DLIBS) -o $@ $(LFLAGS)

face_overlap_rw: face_overlap_rw.o
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $@.o $(PHYSBAM_DLIBS) -o $@ $(LFLAGS)

particles_simple_rw: particles_simple_rw.o
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $@.o $(PHYSBAM_DLIBS) -o $@ $(LFLAGS)

particles_overlap_rw: particles_overlap_rw.o
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $@.o $(PHYSBAM_DLIBS) -o $@ $(LFLAGS)


clean: clean-files
	\rm -f $(TARGET)

