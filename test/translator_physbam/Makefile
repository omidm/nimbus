
TARGET = simple_rw overlap_rw
all: check $(TARGET)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR) -I$(NIMBUS_ROOT)/application/physbam-lib/Public_Library
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

LFLAGS += -lpthread -lnimbus -lboost_thread-mt -lboost_system-mt -lprotobuf
LFLAGS += -lPhysBAM_Dynamics -lPhysBAM_Fluids -lPhysBAM_Geometry -lPhysBAM_Solids
LFLAGS += -lPhysBAM_Tools -lz
LDFLAGS += -L. -L/usr/local/lib
LDFLAGS += -L$(NIMBUS_ROOT)/application/physbam-lib/Public_Library/Build/Release
LDFLAGS += -L$(NIMBUS_ROOT)/application/physbam-lib/External_Libraries/lib/

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

simple_rw: $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) simple_rw.o -o $@ $(LFLAGS)

overlap_rw: $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) overlap_rw.o -o $@ $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

