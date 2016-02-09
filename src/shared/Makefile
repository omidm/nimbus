all: check proto_t shared_t 

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc) 
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

PDIR_SOURCE = protobuf_source/

shared_t: proto_t $(OBJFILES) $(HFILES)

proto_t: 
	cd $(PDIR_SOURCE) && make && cd ..

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	cd $(PDIR_SOURCE) && make clean && cd ..

