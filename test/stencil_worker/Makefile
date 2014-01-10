
TARGET_SINGLE = worker_single
TARGET_MULTI_C = worker_multi_c
TARGET_MULTI = worker_multi
all: check $(TARGET_SINGLE) $(TARGET_MULTI_C) $(TARGET_MULTI)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

ADIR_SINGLE = ${NIMBUS_ROOT}/application/stencil_1d_single/
ADIR_MULTI_C = ${NIMBUS_ROOT}/application/stencil_1d_multi_c/
ADIR_MULTI = ${NIMBUS_ROOT}/application/stencil_1d_multi/

LFLAGS_SINGLE = $(LFLAGS)
LFLAGS_SINGLE += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -lstencil_1d_single -lprotobuf
LDFLAGS_SINGLE = $(LDFLAGS)
LDFLAGS_SINGLE += -L$(ADIR_SINGLE) -Wl,-rpath $(ADIR_SINGLE)

LFLAGS_MULTI_C = $(LFLAGS)
LFLAGS_MULTI_C += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -lstencil_1d_multi_c -lprotobuf
LDFLAGS_MULTI_C = $(LDFLAGS)
LDFLAGS_MULTI_C += -L$(ADIR_MULTI_C) -Wl,-rpath $(ADIR_MULTI_C)

LFLAGS_MULTI = $(LFLAGS)
LFLAGS_MULTI += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -lstencil_1d_multi -lprotobuf
LDFLAGS_MULTI = $(LDFLAGS)
LDFLAGS_MULTI += -L$(ADIR_MULTI) -Wl,-rpath $(ADIR_MULTI)

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

application_single:
	cd ${ADIR_SINGLE}; make clean; make; cd -

application_multi_c:
	cd ${ADIR_MULTI_C}; make clean; make; cd -

application_multi:
	cd ${ADIR_MULTI}; make clean; make; cd -

$(TARGET_SINGLE): application_single $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS_SINGLE) $(LFLAGS_SINGLE) $(OBJECT_FILES) -o $(TARGET_SINGLE) $(LFLAGS_SINGLE)

$(TARGET_MULTI_C): application_multi_c $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS_MULTI_C) $(LFLAGS_MULTI_C) $(OBJECT_FILES) -o $(TARGET_MULTI_C) $(LFLAGS_MULTI_C)

$(TARGET_MULTI): application_multi $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS_MULTI) $(LFLAGS_MULTI) $(OBJECT_FILES) -o $(TARGET_MULTI) $(LFLAGS_MULTI)

clean: clean-files
	\rm -f $(TARGET_SINGLE) $(TARGET_MULTI_C) $(TARGET_MULTI)

