
TARGET = worker
all: check $(TARGET)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

ADIR = ${NIMBUS_ROOT}/application/stencil_1d_single/
 
LFLAGS += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -lstencil_1d_single
LDFLAGS += -L$(ADIR) -Wl,-rpath=$(ADIR)

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

application_t:
	cd ${ADIR}; make; cd -

$(TARGET): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES) application_t
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) $(APP_OBJECT_FILES) -o $(TARGET) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

