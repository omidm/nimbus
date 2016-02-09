all: proto_t data_t

NIMBUS_ROOT = ../../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

PDIR_SOURCE = protobuf_source/

IFLAGS += -I$(NIMBUS_ROOT)/application/physbam-lib/Public_Library

proto_t:
	cd $(PDIR_SOURCE) && make && cd ..

data_t: proto_t $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	cd $(PDIR_SOURCE) && make clean && cd ..
