all: check worker

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS = -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

worker: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
