# make
COMPILER = g++
FLAGS= -Wall -g
LDIR = ../../Library/
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cpp)

# Add the dependencies here
DEPS = $(LDIR)scheduler.o $(LDIR)application.o $(LDIR)job.o $(LDIR)cluster.o $(LDIR)worker.o $(LDIR)data.o

$(LDIR)%.o: $(LDIR)%.cpp $(LDIR)$(wildcard .h)
	$(COMPILER) -c -o $@ $< $(FLAGS)

app: $(CPP_FILES) $(H_FILES) $(DEPS) 
	$(COMPILER) -o app -I $(LDIR) $(CPP_FILES) $(DEPS)

clean:
	rm app
  
clean-all:
	rm app $(LDIR)*.o        
  
#mpi: $(CPP_FILES) $(H_FILES) 
#	mpiCC -o scheduler-mpi $(CPP_FILES) $(INCL)                              
  
