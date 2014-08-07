# Makefile for test application

SERVER = server
CLIENT = client
all: check $(SERVER) $(CLIENT)

NIMBUS_ROOT = ../..
include ${NIMBUS_ROOT}/Makeinclude

LDIR = ../../shared/
IFLAGS += -I$(LDIR)
H_FILES = $(wildcard *.h)
CFILES = $(wildcard *.cc)
OBJECT_FILES = $(subst .cc,.o,$(CFILES))

LFLAGS += -lnimbus -lboost_thread -lboost_system -lpthread
LDFLAGS += -L.

# Add the dependencies here
%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

$(SERVER): $(CPP_FILES) $(H_FILES) $(DEPS) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) server.o -o $(SERVER) $(LFLAGS)


$(CLIENT): $(CPP_FILES) $(H_FILES) $(DEPS) $(OBJECT_FILES)
	$(CPP) $(CFLAGS) $(IFLAGS) $(LDFLAGS) $(LFLAGS) client.o -o $(CLIENT) $(LFLAGS)

clean: clean-files
	\rm -f $(TARGET)

