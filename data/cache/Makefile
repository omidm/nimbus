all: check cache_t

NIMBUS_ROOT = ../..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

cache_t: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	cd $(PDIR_SOURCE); make clean; cd ..
