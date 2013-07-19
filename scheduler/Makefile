# make

NIMBUS_ROOT = ..
include ../Makeinclude

COMPILER = g++
FLAGS= -Wall -g
ROOT = ~/cloud/src/nimbus/
LDIR = $(ROOT)lib/
ADIR = $(ROOT)application/1d-stencil/
LIBS = -lpthread -lboost_system -lboost_thread
LINK = -L /home/omidm/Downloads/boost_1_53_0/stage/lib
INCL = -I /home/omidm/Downloads/boost_1_53_0
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cpp)

# Add the dependencies here
DEPS = $(LDIR)scheduler.o $(LDIR)application.o $(LDIR)job.o $(LDIR)cluster.o $(LDIR)worker.o $(LDIR)data.o $(LDIR)server.o $(LDIR)client.o $(LDIR)parser.o $(ADIR)app.o

$(LDIR)%.o: $(LDIR)%.cpp $(LDIR)%.h
	$(COMPILER) -c -o $@ $< $(FLAGS) $(LIBS) $(LINK) $(INCL)

$(ADIR)%.o: $(ADIR)%.cpp $(ADIR)$(wildcard *.h)
	$(COMPILER) -c -o $@ $< $(FLAGS) -I $(LDIR) $(LIBS) $(LINK) $(INCL)

scheduler: $(CPP_FILES) $(H_FILES) $(DEPS) 
	$(COMPILER) -o scheduler -I $(LDIR) $(CPP_FILES) $(DEPS) $(LIBS) $(LINK) $(INCL)

clean:
	rm scheduler         
  
clean-all:
	rm scheduler $(LDIR)*.o $(ADIR)*.o        
  
#mpi: $(CPP_FILES) $(H_FILES) 
#	mpiCC -o scheduler-mpi $(CPP_FILES) $(INCL)                              
  
