# make
NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

COMPILER = g++
FLAGS= -Wall -g
LDIR = ../../lib/
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)

IFLAGS = -I/opt/local/include 
LFLAGS = 

# Add the dependencies here
DEPS = $(LDIR)scheduler.o $(LDIR)application.o $(LDIR)job.o $(LDIR)cluster.o $(LDIR)worker.o $(LDIR)data.o

$(LDIR)%.o: $(LDIR)%.cpp $(LDIR)$(wildcard .h)
	$(COMPILER) $(IFLAGS) -c -o $@ $< $(FLAGS)

app: $(CPP_FILES) $(H_FILES) $(DEPS) 
	$(COMPILER) $(IFLAGS)  -o app -I $(LDIR) $(CPP_FILES) $(DEPS)

clean:
	rm app
  
clean-all:
	rm app $(LDIR)*.o        
  
#mpi: $(CPP_FILES) $(H_FILES) 
#	mpiCC -o scheduler-mpi $(CPP_FILES) $(INCL)                              
  
