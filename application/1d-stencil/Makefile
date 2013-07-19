# File for 1-d stencil application

all: check app

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

LDIR = ../../lib/
IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)

OBJECT_FILES = app.o
 
LFLAGS += -lboost_thread-mt -lboost_system-mt

# Add the dependencies here
DEPS = $(LDIR)scheduler.o $(LDIR)application.o $(LDIR)job.o $(LDIR)cluster.o $(LDIR)worker.o $(LDIR)data.o $(LDIR)parser.o $(LDIR)server.o

%.o: %.cc
	$(CPP) $(CPP_FLAGS) $(IFLAGS) -c -o $@ $< $(FLAGS)

app: $(CPP_FILES) $(H_FILES) $(DEPS) $(OBJECT_FILES)
	$(CPP) $(CPP_FLAGS) $(LDFLAGS) $(LFLAGS) $(IFLAGS) -o app $(CPP_FILES) $(DEPS)

clean:
	rm app
  
clean-all:
	rm app $(LDIR)*.o        
  
#mpi: $(CPP_FILES) $(H_FILES) 
#	mpiCC -o scheduler-mpi $(CPP_FILES) $(INCL)                              
  
