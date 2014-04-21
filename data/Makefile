all: check data_t physbam cache

NIMBUS_ROOT = ..
include $(NIMBUS_ROOT)/Makeinclude

CFLAGS += -fPIC
CFILES = $(wildcard *.cc)
HFILES = $(wildcard *.h)
OBJFILES = $(subst .cc,.o,$(CFILES))

PHYSBAM_CFILES = $(wildcard physbam/*.cc)
PHYSBAM_OFILES = $(subst .cc,.o,$(PHYSBAM_CFILES))

IFLAGS += -I$(NIMBUS_ROOT)/application/physbam-lib/Public_Library

data_t: $(OBJFILES) $(HFILES)

%.o: %.cc
	$(CPP) $(CFLAGS) $(IFLAGS) -c $< -o $@ 

clean: clean-files
	echo "  cleaning physbam/"; cd physbam; make clean; cd ../cache; make clean; cd ..

.PHONY: physbam
.PHONY: cache

physbam: 
	cd physbam; make; cd ..

cache:
	cd cache; make; cd ..
