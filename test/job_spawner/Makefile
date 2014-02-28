
TARGET = worker
APPLICATION = application

default: check $(TARGET) 
all: check $(APPLICATION) $(TARGET) 

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CPP_FILES))

ADIR = ${NIMBUS_ROOT}/application/job_spawner/

LFLAGS += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread -ljob_spawner -lprotobuf
LDFLAGS += -L$(ADIR) -Wl,-rpath $(ADIR)

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@

$(APPLICATION):
	cd $(ADIR); make clean; make; cd -

$(TARGET): $(CPP_FILES) $(H_FILES) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) -o $(TARGET) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

