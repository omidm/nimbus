
all: check messages shared_t 

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc) 
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))
PROTOCFILES = $(wildcard protobufs/*.cc)
PROTOFILES = $(subst .cc,.o,$(PROTOCFILES))

shared_t: $(OBJFILES) $(HFILES) $(PROTOFILES)

messages: 
	mkdir -p protobufs
	protoc --cpp_out=protobufs *.proto

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	\rm protobufs/*.pb.*

