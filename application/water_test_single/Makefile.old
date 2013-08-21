# make

# nimbus root
NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

# compiler and flags
COMPILER = g++
FLAGS= -Wall -g

# nimbus libraries to link to
LDIR = $(NIMBUS_ROOT)lib/
LIBS = -lpthread -lboost_system -lboost_thread

ADIR = $(NIMBUS_ROOT)/application/water_test/
APP = water_app
H_FILES = $(wildcard *.h)
CPP_FILES = $(wildcard *.cc)

app : $(APP).o

$(APP).o: $(CPP_FILES) $(H_FILES)
	$(COMPILER) -c $(APP).cc -o $(APP).o -I $(NIMBUS_ROOT) $(LIBS)

clean:
	rm *.o
