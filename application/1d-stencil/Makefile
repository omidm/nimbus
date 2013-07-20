# File for 1-d stencil application

TARGET = app
all: check $(TARGET)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

LDIR = ../../lib/
IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)

OBJECT_FILES = app.o main.o
 
LFLAGS += -lnimbus -lboost_thread-mt -lboost_system-mt -lpthread
LDFLAGS += -L.

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

$(TARGET): $(CPP_FILES) $(H_FILES) $(DEPS) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) $(OBJECT_FILES) -o $(TARGET) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

#mpi: $(CPP_FILES) $(H_FILES) 
#	mpiCC -o scheduler-mpi $(CPP_FILES) $(INCL)                              
  
